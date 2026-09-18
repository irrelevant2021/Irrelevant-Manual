#!/usr/bin/env python3
"""
PROTAC Ternary — 全自动流水线 (Schrödinger + RDKit)

合并 prepwizard_fix.py / align_and_extract.py / hybrid_assembly.py /
    optimize_linker_rdkit.py 为一条指令。

用法:
  mamba run -n molscribe python3 pipeline.py \\
      --pdb  ../ppd/TER_0001.pdb \\
      --mae  ../linkers/DY-0417.mae \\
      -o     TER_0001_refined.pdb

依赖:
  - Schrödinger: prepwizard, structconvert, align_ligands
  - RDKit (conda env: molscribe)
"""

import argparse, os, re, shutil, subprocess, sys, tempfile, time


# ═══════════════════════════════════════════════════════════════
#  工具函数
# ═══════════════════════════════════════════════════════════════

def find_schrodinger(schro_path=None):
    for p in (schro_path, os.environ.get('SCHRODINGER', ''),
              '/opt/schrodinger', '/home/zhennan/schrodinger'):
        if p and os.path.isdir(p):
            return p
    print('  ERROR: Schrödinger not found. Set $SCHRODINGER or use --schrodinger')
    sys.exit(1)


def run_cmd(cmd, desc='', logfile=None):
    """Run subprocess, log output."""
    print(f'  [run] {desc or cmd[0]}')
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f'  ERROR: {desc or cmd[0]} failed (code={result.returncode})')
        print(f'  stderr: {result.stderr[:500]}')
        sys.exit(1)
    if logfile:
        with open(logfile, 'a') as f:
            f.write(f'# {desc}\n# CMD: {" ".join(cmd)}\n')
            f.write(result.stdout + '\n')
            if result.stderr:
                f.write('# STDERR:\n' + result.stderr + '\n')
    return result.stdout


def wait_for_job(job_prefix, timeout=120, logfile=None):
    """Wait for Schrödinger background job to complete (poll .log file)."""
    t0 = time.time()
    log_path = f'{job_prefix}.log'
    while time.time() - t0 < timeout:
        if os.path.exists(log_path):
            with open(log_path) as f:
                content = f.read()
            if 'successfully completed' in content:
                if logfile:
                    with open(logfile, 'a') as f:
                        f.write(f'# Job {job_prefix}: completed ({time.time()-t0:.1f}s)\n')
                return True
            if 'failed' in content or 'RuntimeError' in content:
                lines = content.strip().split('\n')
                print(f'  ERROR: job {job_prefix} failed: {lines[-3] if len(lines)>=3 else lines[-1]}')
                return False
        time.sleep(2)
    print(f'  ERROR: job {job_prefix} timed out ({timeout}s)')
    return False
# ═══════════════════════════════════════════════════════════════
#  Step 1: 提取 CBN/POI + prepwizard 修复 + SDF 转换
# ═══════════════════════════════════════════════════════════════

def step1_extract_and_fix(pdb_path, tmpdir, schro, logfile):
    print('\n' + '═' * 60)
    print('  Step 1/4: 提取 CBN/POI + prepwizard 修复')
    print('═' * 60 + '\n')
    prepwiz = f'{schro}/utilities/prepwizard'
    structconv = f'{schro}/utilities/structconvert'
    hetatm = {'CBN': [], 'POI': []}
    conect_all = []
    with open(pdb_path) as f:
        for line in f:
            rec = line[:6].strip()
            if rec == 'HETATM':
                rname = line[17:20].strip()
                if rname in hetatm:
                    hetatm[rname].append(line)
            elif rec == 'CONECT':
                conect_all.append(line)

    # prepwizard requires relative paths; work inside tmpdir
    orig_cwd = os.getcwd()
    os.chdir(tmpdir)

    out_sdfs = {}
    try:
        for lig in ('CBN', 'POI'):
            atoms = hetatm[lig]
            n_atoms = len(atoms)
            if n_atoms == 0:
                print(f'  WARNING: no {lig} found'); continue
            # Write raw PDB (relative path inside tmpdir)
            raw_pdb = f'{lig}_raw.pdb'
            serials = set()
            for ln in atoms:
                try: serials.add(int(ln[6:11]))
                except ValueError: pass
            with open(raw_pdb, 'w') as f:
                f.write(f'HEADER    EXTRACTED FROM {pdb_path}\n')
                f.write(f'HETNAM    {lig.ljust(4)} 1  {lig}\n')
                for ln in atoms: f.write(ln)
                f.write('TER\n')
                for cln in conect_all:
                    parts = cln.split()
                    if len(parts) >= 2:
                        try:
                            if int(parts[1]) in serials: f.write(cln)
                        except ValueError: pass
                f.write('END\n')
            print(f'  {lig}: {n_atoms} atoms -> {raw_pdb}')

            # prepwizard fix (in-place, relative paths)
            fixed_pdb = f'{lig}_fixed.pdb'
            shutil.copy(raw_pdb, fixed_pdb)
            run_cmd([prepwiz, '-WAIT', '-noepik', '-noprotassign',
                      '-nometaltreat', fixed_pdb, fixed_pdb],
                    desc=f'prepwizard {lig}', logfile=logfile)

            # structconvert: PDB -> SDF
            sdf_out = f'{lig}.sdf'
            run_cmd([structconv, fixed_pdb, sdf_out],
                    desc=f'structconvert {lig}', logfile=logfile)
            print(f'  -> {sdf_out}')
            out_sdfs[lig] = os.path.abspath(sdf_out)
    finally:
        os.chdir(orig_cwd)
    return out_sdfs


# ═══════════════════════════════════════════════════════════════
#  Step 2: linker 格式转换 + 双端对齐
# ═══════════════════════════════════════════════════════════════

def parse_sd_props(text):
    """Parse SD tag key-value pairs from SDF text."""
    props = {}
    lines = text.split('\n')
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith('> <') and line.endswith('>'):
            key = line[3:-1]
            i += 1
            vals = []
            while i < len(lines) and lines[i].strip() != '':
                vals.append(lines[i].strip()); i += 1
            props[key] = ' '.join(vals)
        i += 1
    return props


def step2_align(linker_mae, cbn_sdf, poi_sdf, tmpdir, schro,
                core1, core2, logfile):
    print('\n' + '═' * 60)
    print('  Step 2/4: Linker 双端对齐')
    print('═' * 60 + '\n')
    structconv = f'{schro}/utilities/structconvert'
    align = f'{schro}/utilities/align_ligands'

    # Work inside tmpdir for Schrödinger tools (they prefer relative paths)
    # Copy all required SDFs into tmpdir, then use just filenames
    linker_sdf = os.path.join(tmpdir, 'linker.sdf')
    if not os.path.exists(linker_sdf):
        if linker_mae.endswith(('.mae', '.maegz')):
            run_cmd([structconv, linker_mae, linker_sdf],
                    desc='structconvert MAE→SDF', logfile=logfile)
        elif os.path.exists(linker_mae):
            shutil.copy(linker_mae, linker_sdf)
        else:
            print(f'  ERROR: {linker_mae} not found'); sys.exit(1)

    # Reference SDFs should already be in tmpdir from step1
    results = {}
    for label, core in [('CBN', core1), ('POI', core2)]:
        ref_sdf = os.path.join(tmpdir, f'{label}.sdf')
        if not os.path.exists(ref_sdf):
            print(f'  SKIP: {label} ref not available'); continue

        # Schrödinger tools want relative paths - run from tmpdir
        orig_cwd = os.getcwd()
        os.chdir(tmpdir)
        try:
            merged = f'merged_{label}.sdf'
            with open(merged, 'w') as out:
                for fn in ('linker.sdf', f'{label}.sdf'):
                    with open(fn) as f: out.write(f.read())
            aligned = f'aligned_{label}.sdf'
            run_cmd([align, merged, '-ref', '2', '-core', core, '-o', aligned],
                    desc=f'align_ligands -> {label}', logfile=logfile)
            print('    Waiting for job...')
            if not wait_for_job(f'merged_{label}', timeout=120, logfile=logfile):
                sys.exit(1)
        finally:
            os.chdir(orig_cwd)

        # Read result (absolute path)
        aligned_abs = os.path.join(tmpdir, aligned)
        if not os.path.exists(aligned_abs):
            print(f'  ERROR: {aligned} not generated by align_ligands'); sys.exit(1)
        to_out = os.path.join(tmpdir, f'DY-0417_to_{label}.sdf')
        with open(aligned_abs) as f:
            content = f.read()
        idx = content.find('$$$$')
        if idx == -1:
            print(f'  ERROR: no $$$$ in {aligned_abs}'); sys.exit(1)
        with open(to_out, 'w') as f:
            f.write(content[:idx + 4])
        props = parse_sd_props(content[:idx + 4])
        with open(logfile, 'a') as f:
            f.write(f'\n{"="*70}\n  Alignment: DY-0417 -> {label}\n{"="*70}\n')
            for k in ('s_phase_Core_SMARTS', 'r_phase_Similarity_to_Reference',
                      's_phase_Alignment_Method', 's_f3d_frozen_atoms',
                      's_user_Torsion_Dihedral_Atoms'):
                if k in props: f.write(f'  {k}: {props[k]}\n')
        print(f'  -> {to_out}')
        results[label] = to_out
    return results
# ═══════════════════════════════════════════════════════════════
#  Step 3: BFS 域划分 + 坐标杂交 (RDKit)
# ═══════════════════════════════════════════════════════════════

def step3_hybrid_assembly(to_cbn_sdf, to_poi_sdf, tmpdir, logfile):
    print('\n' + '═' * 60)
    print('  Step 3/4: BFS 域划分 + 坐标杂交')
    print('═' * 60 + '\n')
    from rdkit import Chem
    from collections import deque

    def load_mol(p):
        m = Chem.SDMolSupplier(p, removeHs=False)[0]
        if m is None: print(f'  ERROR: cannot read {p}'); sys.exit(1)
        return m
    def parse_frozen(mol, tag='s_f3d_frozen_atoms'):
        raw = mol.GetPropsAsDict().get(tag, '')
        return {int(x)-1 for x in raw.split()} if raw else set()

    mol_cbn = load_mol(to_cbn_sdf)
    mol_poi = load_mol(to_poi_sdf)
    n = mol_cbn.GetNumAtoms()
    print(f'  Atoms: {n}, Bonds: {mol_cbn.GetNumBonds()}')
    for i in range(n):
        if mol_cbn.GetAtomWithIdx(i).GetAtomicNum() != mol_poi.GetAtomWithIdx(i).GetAtomicNum():
            print(f'  ERROR: atom {i} mismatch'); sys.exit(1)
    print('  Atom consistency: OK')

    cbn_core = parse_frozen(mol_cbn)
    poi_core = parse_frozen(mol_poi)
    print(f'  CBN core: {len(cbn_core)}, POI core: {len(poi_core)}')

    # BFS classification
    adj = [[] for _ in range(n)]
    for b in mol_cbn.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        adj[i].append(j); adj[j].append(i)

    def bfs(srcs):
        d = [-1]*n
        q = deque(srcs)
        for s in srcs: d[s] = 0
        while q:
            u = q.popleft()
            for v in adj[u]:
                if d[v]==-1: d[v]=d[u]+1; q.append(v)
        return d
    dc = bfs(cbn_core); dp = bfs(poi_core)
    L = min(dp[c] for c in cbn_core)
    is_spine = [dc[i]+dp[i]==L for i in range(n)]
    is_linker = is_spine[:]
    for i in range(n):
        if is_spine[i]:
            for v in adj[i]:
                if v not in cbn_core and v not in poi_core: is_linker[v]=True
    domain = ['']*n
    for i in cbn_core: domain[i]='CBN'
    for i in poi_core: domain[i]='POI'
    for i in range(n):
        if domain[i]=='':
            domain[i]='Linker' if is_linker[i] else ('CBN' if dc[i]<dp[i] else 'POI')
    cnt = {d: domain.count(d) for d in ('CBN','Linker','POI')}
    for d,c in cnt.items(): print(f'  {d}: {c}')

    # Coordinate hybridization
    mol = Chem.Mol(mol_poi)
    cc = mol_cbn.GetConformer(); cp = mol_poi.GetConformer()
    cn = mol.GetConformer()
    src = {'CBN':cc,'Linker':cc,'POI':cp}
    for i in range(n):
        cn.SetAtomPosition(i, src[domain[i]].GetAtomPosition(i))

    # Strip Hs + write
    old2new = {}
    ni = 0
    for oi in range(mol.GetNumAtoms()):
        if mol.GetAtomWithIdx(oi).GetAtomicNum()!=1:
            old2new[oi]=ni; ni+=1
    mol_noh = Chem.RemoveHs(mol)
    code = {'CBN':'C','Linker':'L','POI':'P'}
    ds = ''.join(code[d] for d in domain)
    nds = ''.join(ds[oi] for oi in sorted(old2new))
    mol_noh.SetProp('s_user_Protac_Domain', nds)
    mol_noh.SetProp('s_user_Protac_Domain_Counts',
                    f'CBN={nds.count("C")} Linker={nds.count("L")} POI={nds.count("P")}')
    hybrid = os.path.join(tmpdir, 'hybrid.sdf')
    Chem.SDWriter(hybrid).write(mol_noh)
    print(f'  -> {hybrid}  ({mol_noh.GetNumAtoms()} heavy atoms)')
    return hybrid, nds
# ═══════════════════════════════════════════════════════════════
#  Step 4: 扭转扰动 + MMFF + PDB 嵌入 (RDKit)
# ═══════════════════════════════════════════════════════════════

def step4_optimize_and_embed(hybrid_sdf, domain_str, orig_pdb, output,
                              n_confs, seed, resname, chain, resseq,
                              tmpdir, logfile):
    print('\n' + '═' * 60)
    print('  Step 4/4: Linker 优化 + PDB 嵌入')
    print('═' * 60 + '\n')
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolTransforms
    from rdkit.Geometry import Point3D
    import math, random, time
    random.seed(seed)
    C, P, L = 'C', 'P', 'L'

    mol = Chem.SDMolSupplier(hybrid_sdf, removeHs=False)[0]
    print(f'  Hybrid: {mol.GetNumAtoms()} heavy atoms')
    mol_h = Chem.AddHs(mol, addCoords=True)
    print(f'  With Hs: {mol_h.GetNumAtoms()}')

    heavy_to_all = {}
    hvy = 0
    for i in range(mol_h.GetNumAtoms()):
        if mol_h.GetAtomWithIdx(i).GetAtomicNum()!=1:
            heavy_to_all[hvy]=i; hvy+=1

    fixed = set()
    for oi,ch in enumerate(domain_str):
        if ch in (C,P):
            ai = heavy_to_all[oi]; fixed.add(ai)
            for nb in mol_h.GetAtomWithIdx(ai).GetNeighbors():
                if nb.GetAtomicNum()==1: fixed.add(nb.GetIdx())
    fixed = sorted(fixed)
    print(f'  Fixed atoms: {len(fixed)}')
    conf0 = mol_h.GetConformer()
    props = AllChem.MMFFGetMoleculeProperties(mol_h, mmffVariant='MMFF94')
    ff = AllChem.MMFFGetMoleculeForceField(mol_h, props) if props else \
         AllChem.UFFGetMoleculeForceField(mol_h)
    fn = 'MMFF94' if props else 'UFF'
    for idx in fixed: ff.AddFixedPoint(idx)
    t0=time.time(); ff.Minimize(maxIts=2000); e0=ff.CalcEnergy()
    print(f'  Repair ({fn}): energy={e0:.1f} ({time.time()-t0:.1f}s)')

    dlook = {}
    for hv,ai in heavy_to_all.items():
        if hv<len(domain_str): dlook[ai]=domain_str[hv]
    rb = []
    for b in mol_h.GetBonds():
        i,j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if b.IsInRing() or b.GetBondType()!=Chem.BondType.SINGLE: continue
        if mol_h.GetAtomWithIdx(i).GetAtomicNum()==1 or mol_h.GetAtomWithIdx(j).GetAtomicNum()==1: continue
        if dlook.get(i)==L and dlook.get(j)==L:
            if i not in fixed and j not in fixed: rb.append((i,j))
    print(f'  Rotatable bonds: {len(rb)}')

    all_e = []; best_mol=None; best_e=float('inf')
    bs = max(50, n_confs//10)
    for it in range(n_confs):
        mc = Chem.Mol(mol_h); cf = mc.GetConformer()
        of = {idx: Point3D(*cf.GetAtomPosition(idx)) for idx in fixed}
        for i,j in rb:
            ai=mc.GetAtomWithIdx(i); aj=mc.GetAtomWithIdx(j)
            inb=[n.GetIdx() for n in ai.GetNeighbors() if n.GetIdx()!=j]
            jnb=[n.GetIdx() for n in aj.GetNeighbors() if n.GetIdx()!=i]
            if not inb or not jnb: continue
            ip=inb[0]
            for n in inb:
                if mc.GetAtomWithIdx(n).GetAtomicNum()!=1: ip=n; break
            jp=jnb[0]
            for n in jnb:
                if mc.GetAtomWithIdx(n).GetAtomicNum()!=1: jp=n; break
            try: rdMolTransforms.SetDihedralDeg(cf,ip,i,j,jp,random.uniform(0,360))
            except: pass
        props=AllChem.MMFFGetMoleculeProperties(mc,mmffVariant='MMFF94')
        ff=AllChem.MMFFGetMoleculeForceField(mc,props) if props else AllChem.UFFGetMoleculeForceField(mc)
        if ff is None: continue
        for idx in fixed: ff.AddFixedPoint(idx)
        ff.Minimize(maxIts=300)
        cf=mc.GetConformer()
        for idx in fixed: cf.SetAtomPosition(idx,of[idx])
        props2=AllChem.MMFFGetMoleculeProperties(mc,mmffVariant='MMFF94')
        ff2=AllChem.MMFFGetMoleculeForceField(mc,props2) if props2 else AllChem.UFFGetMoleculeForceField(mc)
        if ff2 is None: continue
        for idx in fixed: ff2.AddFixedPoint(idx)
        ff2.Minimize(maxIts=200)
        e=ff2.CalcEnergy()
        all_e.append(e)
        if e<best_e: best_e=e; best_mol=mc
        if (it+1)%bs==0 or it==n_confs-1:
            print(f'    [{it+1:>4d}/{n_confs}] {len(all_e)} ok, best={best_e:.1f}', end='\r')
    print()
    if best_mol is None: print('  ERROR: all failed'); sys.exit(1)
    all_e.sort()
    print(f'  Energy: range={all_e[0]:.1f}..{all_e[-1]:.1f}, best={best_e:.1f}, '
          f'mean={sum(all_e)/len(all_e):.1f}, median={all_e[len(all_e)//2]:.1f}')
    elog = os.path.splitext(output)[0]+'.energies.txt'
    with open(elog,'w') as f:
        f.write(f'# Conformer energies\n# Total: {len(all_e)}\n# Range: {all_e[0]:.1f}-{all_e[-1]:.1f}\n'
                f'# Best: {best_e:.1f}\n# Format: rank energy(kJ/mol)\n')
        for rk,ev in enumerate(all_e,1): f.write(f'{rk:>6d}  {ev:>12.3f}\n')
    print(f'  Energy log: {elog} ({len(all_e)} entries)')

    mh = Chem.RemoveHs(best_mol)
    hl = []
    cf = mh.GetConformer()
    for ai in range(mh.GetNumAtoms()):
        at = mh.GetAtomWithIdx(ai); el = at.GetSymbol(); p = cf.GetAtomPosition(ai)
        hl.append(f'HETATM{ai+1:5d} {f"{el}{ai+1}".ljust(4)[:4]} '
                  f'{resname:3s} {chain}{resseq:4d}    '
                  f'{p.x:8.3f}{p.y:8.3f}{p.z:8.3f}{1.00:6.2f}{0.00:6.2f}          {el:>2s}')
    cn = {}
    for b in mh.GetBonds():
        i=b.GetBeginAtomIdx()+1; j=b.GetEndAtomIdx()+1
        cn.setdefault(i,[]).append(j); cn.setdefault(j,[]).append(i)
    cl = []
    for s in sorted(cn):
        nbrs=sorted(cn[s])
        for ck in range(0,len(nbrs),4):
            ll=f'CONECT{s:5d}'
            for nb in nbrs[ck:ck+4]: ll+=f'{nb:5d}'
            cl.append(ll)

    with open(orig_pdb) as f: lines=f.readlines()
    li=[i for i,l in enumerate(lines) if l.startswith('HETATM') and l[17:20].strip() in ('CBN','POI')]
    if not li: print('  ERROR: no CBN/POI in PDB'); sys.exit(1)
    first,last=li[0],li[-1]
    lines=[l for l in lines if not l.startswith('CONECT')]
    ip=len(lines)
    for i in range(len(lines)-1,-1,-1):
        if lines[i].startswith('TER'): ip=i+1; break
    out = (lines[:first]+['\n'.join(hl)+'\n']+lines[last+1:ip]
           +['\n'.join(cl)+'\n']+lines[ip:])
    with open(output,'w') as f: f.writelines(out)
    print(f'  -> {output} ({len(hl)} HETATM, {len(cl)} CONECT)')
# ═══════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description='PROTAC 全自动流水线 (Schrödinger + RDKit)')
    parser.add_argument('--pdb', required=True, help='Rosetta 三元 PDB')
    parser.add_argument('--mae', required=True, help='PROTAC linker MAE 文件')
    parser.add_argument('-o','--output', required=True, help='输出 PDB')
    parser.add_argument('--cbn-core', default='MCS', help='CBN 对齐核心')
    parser.add_argument('--poi-core', default='N#Cc1ccccc1', help='POI 对齐核心')
    parser.add_argument('--n-confs', type=int, default=200, help='扭转采样数')
    parser.add_argument('--seed', type=int, default=42, help='随机种子')
    parser.add_argument('--schrodinger', help='Schrödinger 路径')
    parser.add_argument('--keep-temp', action='store_true', help='保留中间文件')
    parser.add_argument('--resname', default='PRT', help='PDB 残基名')
    parser.add_argument('--chain', default='X', help='PDB 链')
    parser.add_argument('--resseq', type=int, default=1, help='PDB 残基序号')
    args = parser.parse_args()

    t_start = time.time()
    schro = find_schrodinger(args.schrodinger)
    tmpdir = os.path.abspath(tempfile.mkdtemp(prefix='pipeline_', dir='.'))
    logfile = os.path.join(tmpdir, 'pipeline.log')
    print(f'  Temp dir: {tmpdir}\n  Schrödinger: {schro}')
    print(f'  Input: {args.pdb}\n  Output: {args.output}\n')

    # Convert I/O paths to absolute (steps chdir into tmpdir)
    pdb_abs = os.path.abspath(args.pdb)
    mae_abs = os.path.abspath(args.mae)
    out_abs = os.path.abspath(args.output)

    try:
        sdfs = step1_extract_and_fix(pdb_abs, tmpdir, schro, logfile)
        aligned = step2_align(mae_abs, sdfs.get('CBN'), sdfs.get('POI'),
                              tmpdir, schro, args.cbn_core, args.poi_core, logfile)
        hybrid_sdf, ds = step3_hybrid_assembly(
            aligned.get('CBN'), aligned.get('POI'), tmpdir, logfile)
        step4_optimize_and_embed(hybrid_sdf, ds, pdb_abs, out_abs,
                                  args.n_confs, args.seed, args.resname,
                                  args.chain, args.resseq, tmpdir, logfile)
        print(f'\n{"="*60}\n  Pipeline 完成! {time.time()-t_start:.1f}s\n'
              f'  输出: {out_abs}\n{"="*60}')
    finally:
        if not args.keep_temp:
            shutil.rmtree(tmpdir, ignore_errors=True)
        else:
            print(f'  (temp: {tmpdir})')

if __name__ == '__main__':
    main()