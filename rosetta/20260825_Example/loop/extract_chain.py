import pymol
from pymol import cmd, stored
import subprocess
import os
import glob

def init_pymol():
    """初始化 PyMOL（只调用一次）"""
    pymol.finish_launching(['pymol', '-cq'])

def quit_pymol():
    """退出 PyMOL"""
    cmd.delete('all')
    cmd.quit()

def extract_and_renumber(pdb_file, chain_id, output_file):
    """提取指定链并重新编号 (pymol)"""
    
    # 加载 PDB
    cmd.load(pdb_file, 'protein')
    
    # 提取指定链
    selection_name = f'chain_{chain_id}'
    cmd.extract(selection_name, f'chain {chain_id}')
    
    # 重新编号从1开始
    stored.residues = []
    cmd.iterate(selection_name, 'stored.residues.append(resi)')
    
    # 创建重编号映射
    residue_map = {}
    for i, old_resi in enumerate(sorted(set(stored.residues), key=int)):
        residue_map[old_resi] = str(i + 1)
    
    # 应用新编号
    for old_resi, new_resi in residue_map.items():
        cmd.alter(f'{selection_name} and resi {old_resi}', f'resi="{new_resi}"')
    
    # 排序
    cmd.sort(selection_name)
    
    # 保存
    cmd.save(output_file, selection_name)
    
    print(f"提取的链 {chain_id} 已保存到 {output_file}")
    print(f"共 {len(residue_map)} 个残基，已从 1 开始重新编号")
    
    # 清理
    cmd.delete('all')

def fix_missing_residues(input_pdb, output_pdb, missing_dict):
    """使用 pdbfixer 补全缺失残基 (openfe_env)"""
    # 生成临时 python 脚本内容
    script_content = f"""
from openmm.app import PDBFile
from pdbfixer import PDBFixer

fixer = PDBFixer(filename='{input_pdb}')
fixer.findMissingResidues()
fixer.missingResidues = {missing_dict}
fixer.findMissingAtoms()
fixer.addMissingAtoms()
PDBFile.writeFile(fixer.topology, fixer.positions, open('{output_pdb}', 'w'))
print(f"补全缺失残基完成，输出到 {output_pdb}")
"""
    # 用 openfe_env 的 python 执行
    openfe_python = '/home/xfusion/miniforge3/envs/openfe_env/bin/python'
    result = subprocess.run(
        [openfe_python, '-c', script_content],
        capture_output=True, text=True
    )
    print(result.stdout)
    if result.stderr:
        print(result.stderr)

def merge_chain_back(original_pdb, new_chain_pdb, chain_id, output_pdb):
    """将补全后的链放回原结构，替换原来的链"""
    
    # 加载原始结构
    cmd.load(original_pdb, 'original')
    
    # 加载补全后的链（pdbfixer输出链名为A，需要改回C）
    cmd.load(new_chain_pdb, 'new_chain')
    
    # 将新链的链标识符改为原来的链标识符
    cmd.alter('new_chain', f'chain="{chain_id}"')
    
    # 删除原始结构中的旧链
    cmd.remove(f'original and chain {chain_id}')
    
    # 将新链合并到original
    cmd.create('merged', 'original or new_chain')
    
    # 保存
    cmd.save(output_pdb, 'merged')
    
    print(f"已将补全的 {chain_id} 链放回原结构")
    print(f"输出到 {output_pdb}")
    
    # 清理
    cmd.delete('all')

# 主流程
build_dir = '/data/AR/PROTAC_ternary-master/20260825/build'
loop_dir = '/data/AR/PROTAC_ternary-master/20260825/loop'

# 查找所有 *mod*.pdb 文件
mod_files = sorted(glob.glob(os.path.join(build_dir, '*mod*.pdb')))
print(f"找到 {len(mod_files)} 个文件需要处理")

# 初始化 PyMOL（只调用一次）
init_pymol()

for mod_pdb in mod_files:
    # 提取文件名（不含路径和扩展名）
    basename = os.path.splitext(os.path.basename(mod_pdb))[0]  # e.g., ternary0_mod
    
    chainC_pdb = os.path.join(loop_dir, f'{basename}_chainC.pdb')
    loop_pdb = os.path.join(loop_dir, f'{basename}_chainC_loop.pdb')
    mod_loop_pdb = os.path.join(loop_dir, f'{basename}_loop.pdb')
    
    print(f"\n{'='*50}")
    print(f"处理: {basename}")
    print(f"{'='*50}")
    
    # 步骤1: 提取链并重编号
    extract_and_renumber(mod_pdb, 'C', chainC_pdb)
    
    # 步骤2: 补全缺失残基
    missing_residues = {(0, 210): ['SER', 'VAL', 'ASP', 'PHE', 'PRO', 'GLU']}
    fix_missing_residues(chainC_pdb, loop_pdb, missing_residues)
    
    # 步骤3: 将补全的C链放回原结构
    merge_chain_back(mod_pdb, loop_pdb, 'C', mod_loop_pdb)

# 退出 PyMOL
quit_pymol()

print(f"\n{'='*50}")
print(f"全部完成！共处理 {len(mod_files)} 个文件")
print(f"输出目录: {loop_dir}")
