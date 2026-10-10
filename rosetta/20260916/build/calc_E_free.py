#!/usr/bin/env python3
"""
calc_E_free.py — 计算每个 PROTAC 分子的自由态 MMFF94 能量 (E_free)

用途:
  为 summary.csv 中的每行计算 ΔE = Best_Energy_kJ - E_free，消除不同 linker(完整PROTAC)
  化学基线差异，使能量特征跨分子可比。

使用方式 (与 batch_pipeline.sh 环境一致):
  cd build
  mamba run -n openfe_env python3 calc_E_free.py \\
      --mae-dir ../linkers \\
      --schrodinger /opt/schrodinger \\
      --output-dir ../outputs

输出:
  {output_dir}/E_free.csv           — 每个 MAE 的 E_free 值
  {output_dir}/enhanced_summary.csv — summary.csv + 新增 ΔE 列

依赖:
  - Schrödinger: structconvert (MAE → SDF)
  - RDKit: 读取 SDF + MMFF94 力场最小化
  - 运行环境: openfe_env (同 pipeline.py)
"""

import argparse, csv, os, shutil, subprocess, sys, tempfile, time


# ═══════════════════════════════════════════════════════════════
#  工具函数 (复用 pipeline.py 的风格)
# ═══════════════════════════════════════════════════════════════

def find_schrodinger(schro_path=None):
    """定位 Schrödinger 安装路径 — 同 pipeline.py 的 find_schrodinger"""
    for p in (schro_path, os.environ.get('SCHRODINGER', ''),
              '/opt/schrodinger', '/home/zhennan/schrodinger'):
        if p and os.path.isdir(p):
            return p
    print('  ERROR: Schrödinger not found. Set $SCHRODINGER or use --schrodinger')
    sys.exit(1)


def run_cmd(cmd, desc='', log=None):
    """运行外部命令，失败时退出"""
    print(f'  [run] {desc or cmd[0]}')
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f'  ERROR: {desc or cmd[0]} failed (code={result.returncode})')
        print(f'  stderr: {result.stderr[:500]}')
        sys.exit(1)
    if log:
        with open(log, 'a') as f:
            f.write(f'# {desc}\n# CMD: {" ".join(cmd)}\n')
            f.write(result.stdout + '\n')
    return result.stdout


# ═══════════════════════════════════════════════════════════════
#  核心: 从 MAE 文件计算 E_free
# ═══════════════════════════════════════════════════════════════

def calc_E_free_from_mae(mae_path, schro, tmpdir, n_confs=100, seed=42):
    """
    对单个 MAE 文件:
      1. structconvert MAE → SDF
      2. RDKit 读入、加氢
      3. MMFF94 全自由最小化 (无任何固定点)
      4. 多次随机扭转采样 + 最小化, 取全局最低 → E_free

    返回: (mae_basename, E_free_value or None)
    """
    import random
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolTransforms

    random.seed(seed)
    mae_base = os.path.splitext(os.path.basename(mae_path))[0]
    structconv = f'{schro}/utilities/structconvert'

    # ── 1. MAE → SDF (同 pipeline.py line 177-179) ──
    sdf_out = os.path.join(tmpdir, f'{mae_base}.sdf')
    if not os.path.exists(sdf_out):
        run_cmd([structconv, mae_path, sdf_out],
                desc=f'structconvert {mae_base}', log=None)
    else:
        print(f'  [skip] {sdf_out} 已存在')

    # ── 2. RDKit 读取 ──
    mol = Chem.SDMolSupplier(sdf_out, removeHs=False)[0]
    if mol is None:
        print(f'  ERROR: 无法读取 {sdf_out}')
        return mae_base, None

    n_heavy = mol.GetNumAtoms()
    print(f'  {mae_base}: {n_heavy} atoms (from SDF)')

    # ── 3. 加氢 ──
    mol_h = Chem.AddHs(mol, addCoords=True)
    n_total = mol_h.GetNumAtoms()
    print(f'           {n_total} atoms (after AddHs, {n_total - n_heavy} H)')

    # ── 4. 获取 MMFF94 力场属性 ──
    props = AllChem.MMFFGetMoleculeProperties(mol_h, mmffVariant='MMFF94')
    if props is None:
        print(f'  WARN: MMFF94 参数化失败, 尝试 UFF 预处理后重试')
        ff_uff = AllChem.UFFGetMoleculeForceField(mol_h)
        if ff_uff is None:
            print(f'  ERROR: UFF 也无法参数化, 跳过 {mae_base}')
            return mae_base, None
        ff_uff.Minimize(maxIts=500)
        props = AllChem.MMFFGetMoleculeProperties(mol_h, mmffVariant='MMFF94')
        if props is None:
            print(f'  ERROR: 重试后 MMFF94 仍无法参数化, 跳过 {mae_base}')
            return mae_base, None

    # ── 5. 从起始结构直接最小化 (同 pipeline.py line 362, 2000步) ──
    ff = AllChem.MMFFGetMoleculeForceField(mol_h, props)
    ff.Minimize(maxIts=2000)
    e_initial = ff.CalcEnergy()
    best_e = e_initial
    print(f'    [起始结构]      E = {e_initial:.3f}')

    # ── 6. 随机扭转采样 + 最小化 (同 pipeline.py lines 379-409, 但无固定点) ──
    #    识别可旋转键 (非环、非末端、非H的单键)
    rot_bonds = []
    for b in mol_h.GetBonds():
        if b.IsInRing() or b.GetBondType() != Chem.BondType.SINGLE:
            continue
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if mol_h.GetAtomWithIdx(i).GetAtomicNum() == 1 or \
           mol_h.GetAtomWithIdx(j).GetAtomicNum() == 1:
            continue
        # 排除末端原子 (邻居数==1)
        if mol_h.GetAtomWithIdx(i).GetDegree() == 1 or \
           mol_h.GetAtomWithIdx(j).GetDegree() == 1:
            continue
        rot_bonds.append((i, j))

    print(f'    Rotatable bonds: {len(rot_bonds)}')

    if len(rot_bonds) == 0:
        print(f'    → E_free = {best_e:.3f}  (无可旋转键)')
        return mae_base, best_e

    batch_size = max(10, n_confs // 10)
    n_ok = 0

    for it in range(n_confs):
        mc = Chem.Mol(mol_h)
        cf = mc.GetConformer()

        # 随机化所有可旋转键二面角
        for i, j in rot_bonds:
            try:
                inb = [n.GetIdx() for n in mc.GetAtomWithIdx(i).GetNeighbors()
                       if n.GetIdx() != j]
                jnb = [n.GetIdx() for n in mc.GetAtomWithIdx(j).GetNeighbors()
                       if n.GetIdx() != i]
                if not inb or not jnb:
                    continue
                # 优先用重原子作为参考
                ip = inb[0]
                for n in inb:
                    if mc.GetAtomWithIdx(n).GetAtomicNum() != 1:
                        ip = n; break
                jp = jnb[0]
                for n in jnb:
                    if mc.GetAtomWithIdx(n).GetAtomicNum() != 1:
                        jp = n; break
                rdMolTransforms.SetDihedralDeg(cf, ip, i, j, jp,
                                                random.uniform(0, 360))
            except Exception:
                continue

        # MMFF94 最小化 (无固定点! 全自由)
        props2 = AllChem.MMFFGetMoleculeProperties(mc, mmffVariant='MMFF94')
        if props2 is None:
            continue
        ff2 = AllChem.MMFFGetMoleculeForceField(mc, props2)
        if ff2 is None:
            continue
        ff2.Minimize(maxIts=500)
        e = ff2.CalcEnergy()

        n_ok += 1
        if e < best_e:
            best_e = e

        if (it + 1) % batch_size == 0 or it == n_confs - 1:
            print(f'    [{it+1:>4d}/{n_confs}] ok={n_ok}, best={best_e:.3f}', end='\r')

    print()
    print(f'    → E_free({mae_base}) = {best_e:.3f}  '
          f'(来自 {n_ok}/{n_confs} 个有效构象)')
    return mae_base, best_e


# ═══════════════════════════════════════════════════════════════
#  汇总输出: E_free.csv + enhanced_summary.csv (含 ΔE)
# ═══════════════════════════════════════════════════════════════

def load_E_free(e_free_csv):
    """读取 E_free.csv, 返回 {MAE: E_free} 字典"""
    result = {}
    with open(e_free_csv, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            mae = row['MAE']
            try:
                result[mae] = float(row['E_free'])
            except (ValueError, KeyError):
                result[mae] = None
    return result


def write_E_free_csv(e_free_records, output_path):
    """写入 E_free.csv"""
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['MAE', 'E_free'])
        for mae, val in e_free_records:
            if val is not None:
                writer.writerow([mae, f'{val:.3f}'])
            else:
                writer.writerow([mae, ''])

    n_ok = sum(1 for _, v in e_free_records if v is not None)
    n_total = len(e_free_records)
    print(f'\n  ✅ E_free.csv 已生成: {output_path}  ({n_ok}/{n_total} 成功)')


def generate_enhanced_summary(summary_csv, e_free_csv, output_csv):
    """
    读取 summary.csv 和 E_free.csv,
    添加 ΔE 和 E_free 列后输出 enhanced_summary.csv
    """
    e_free_map = load_E_free(e_free_csv)

    with open(summary_csv, 'r') as f_in, \
         open(output_csv, 'w', newline='') as f_out:

        reader = csv.DictReader(f_in)
        fieldnames = reader.fieldnames[:]
        for col in ('E_free', 'ΔE'):
            if col not in fieldnames:
                fieldnames.append(col)

        writer = csv.DictWriter(f_out, fieldnames=fieldnames)
        writer.writeheader()

        n_rows = 0
        n_with_delta = 0
        for row in reader:
            mae = row.get('MAE', '')
            try:
                best = float(row.get('Best_Energy_kJ', ''))
            except (ValueError, TypeError):
                best = None

            e_free = e_free_map.get(mae)

            if best is not None and e_free is not None:
                delta = best - e_free
                row['ΔE'] = f'{delta:.3f}'
                row['E_free'] = f'{e_free:.3f}'
                n_with_delta += 1
            else:
                row['ΔE'] = ''
                row['E_free'] = ''

            writer.writerow(row)
            n_rows += 1

    print(f'  ✅ enhanced_summary.csv 已生成: {output_csv}')
    print(f'     {n_rows} 行, 其中 {n_with_delta} 行有 ΔE 值')


# ═══════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description='计算 PROTAC 自由态 MMFF94 能量 (E_free) 及 ΔE')
    parser.add_argument('--mae-dir', required=True,
                        help='MAE 文件目录 (例: ../linkers)')
    parser.add_argument('--schrodinger',
                        help='Schrödinger 安装路径 (默认: $SCHRODINGER 或自动查找)')
    parser.add_argument('--output-dir', default='.',
                        help='输出目录 (默认: 当前目录)')
    parser.add_argument('--n-confs', type=int, default=100,
                        help='随机构象采样数 (默认: 100)')
    parser.add_argument('--seed', type=int, default=42,
                        help='随机种子 (默认: 42, 同 pipeline.py)')
    parser.add_argument('--summary', default=None,
                        help='summary.csv 路径 (默认: {output_dir}/summary.csv)')
    parser.add_argument('--keep-temp', action='store_true',
                        help='保留临时目录')
    args = parser.parse_args()

    t_start = time.time()

    # ── 定位 Schrödinger ──
    schro = find_schrodinger(args.schrodinger)
    print(f'  Schrödinger: {schro}')

    # ── 解析 MAE 目录 ──
    mae_dir = os.path.abspath(args.mae_dir)
    if not os.path.isdir(mae_dir):
        print(f'  ERROR: MAE 目录不存在: {mae_dir}')
        sys.exit(1)

    mae_files = sorted([
        os.path.join(mae_dir, f)
        for f in os.listdir(mae_dir)
        if f.endswith('.mae')
    ])
    if not mae_files:
        print(f'  ERROR: {mae_dir} 中没有 .mae 文件')
        sys.exit(1)
    print(f'  找到 {len(mae_files)} 个 MAE 文件: '
          f'{", ".join(os.path.basename(m) for m in mae_files)}')

    # ── 输出目录 ──
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # ── 逐 MAE 计算 E_free ──
    tmpdir = os.path.abspath(tempfile.mkdtemp(prefix='E_free_', dir='.'))
    print(f'  Temp dir: {tmpdir}\n')

    e_free_records = []
    for i, mae_path in enumerate(mae_files):
        mae_name = os.path.splitext(os.path.basename(mae_path))[0]
        print(f'\n[{i+1}/{len(mae_files)}] {mae_name}')
        print('─' * 50)
        base, val = calc_E_free_from_mae(
            mae_path, schro, tmpdir,
            n_confs=args.n_confs, seed=args.seed + i
        )
        e_free_records.append((base, val))

    # ── 输出 E_free.csv ──
    e_free_csv = os.path.join(output_dir, 'E_free.csv')
    write_E_free_csv(e_free_records, e_free_csv)

    # ── 输出 enhanced_summary.csv ──
    summary_csv = args.summary or os.path.join(output_dir, 'summary.csv')
    if os.path.exists(summary_csv):
        enhanced_csv = os.path.join(output_dir, 'enhanced_summary.csv')
        generate_enhanced_summary(summary_csv, e_free_csv, enhanced_csv)
    else:
        print(f'\n  ⚠ 未找到 summary.csv: {summary_csv}')
        print('    跳过 enhanced_summary.csv 生成')

    # ── 清理 ──
    if not args.keep_temp:
        shutil.rmtree(tmpdir, ignore_errors=True)
    else:
        print(f'  (temp: {tmpdir})')

    elapsed = time.time() - t_start
    print(f'\n  ⏱ 总耗时: {elapsed:.1f}s')
    print(f'  📄 E_free.csv         → {e_free_csv}')
    if os.path.exists(os.path.join(output_dir, 'enhanced_summary.csv')):
        print(f'  📄 enhanced_summary.csv → '
              f'{os.path.join(output_dir, "enhanced_summary.csv")}')
    print()


if __name__ == '__main__':
    main()