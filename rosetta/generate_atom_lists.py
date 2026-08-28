#!/usr/bin/env python3
"""
批量自动生成 PROTAC 三元复合物建模所需的原子列表文件

使用场景：处理大量 linker PDB 文件，自动化生成原子列表
"""

import argparse
import sys
from pathlib import Path
from rdkit import Chem
import glob


def extract_residue_heavy_atoms(mol, resname):
    """提取指定残基的所有重原子"""
    atoms = []
    for i, atom in enumerate(mol.GetAtoms()):
        info = atom.GetPDBResidueInfo()
        if info and info.GetResidueName().strip() == resname and atom.GetSymbol() != 'H':
            aname = info.GetName().strip()
            atoms.append((resname, aname))
    return atoms


def find_substructure_in_residue(mol, resname, smiles):
    """在指定残基中搜索子结构，返回匹配的原子名列表
    
    Args:
        mol: RDKit Mol 对象
        resname: 残基名称
        smiles: 要搜索的子结构 SMILES
    
    Returns:
        list: 匹配的原子名列表 [(resname, atomname), ...]
    """
    # 提取该残基的所有原子索引
    residue_atom_indices = []
    for i, atom in enumerate(mol.GetAtoms()):
        info = atom.GetPDBResidueInfo()
        if info and info.GetResidueName().strip() == resname:
            residue_atom_indices.append(i)
    
    if not residue_atom_indices:
        return []
    
    # 创建子结构 SMILES 分子
    submol = Chem.MolFromSmiles(smiles)
    if submol is None:
        return []
    
    # 在残基原子中搜索匹配
    matched_atoms = []
    matches = mol.GetSubstructMatches(submol)
    
    for match in matches:
        # 检查匹配的原子是否都在该残基中
        if all(idx in residue_atom_indices for idx in match):
            for idx in match:
                atom = mol.GetAtomWithIdx(idx)
                info = atom.GetPDBResidueInfo()
                if info:
                    aname = info.GetName().strip()
                    matched_atoms.append((resname, aname))
    
    # 去重
    return list(set(matched_atoms))


def find_matching_atoms(linker_mol, ligand_smiles):
    """在 linker 中搜索与配体匹配的子结构
    
    Returns:
        tuple: (match_count, matched_atom_names)
            - match_count: 匹配次数
            - matched_atom_names: 匹配的原子名列表
    """
    ligand = Chem.MolFromSmiles(ligand_smiles)
    if ligand is None:
        return 0, []
    
    # 尝试匹配
    matches = linker_mol.GetSubstructMatches(ligand)
    
    if not matches:
        # 尝试 sanitize 后匹配
        try:
            linker_sanitized = Chem.Mol(linker_mol)
            Chem.SanitizeMol(linker_sanitized)
            matches = linker_sanitized.GetSubstructMatches(ligand)
        except:
            return 0, []
    
    # 收集匹配的原子名
    matched_atom_names = set()
    for match in matches:
        for idx in match:
            atom = linker_mol.GetAtomWithIdx(idx)
            info = atom.GetPDBResidueInfo()
            if info:
                aname = info.GetName().strip()
                matched_atom_names.add(aname)
    
    return len(matches), sorted(list(matched_atom_names))


def process_single_linker(decoy_mol, linker_mol, lig1_name, lig2_name, 
                         lig1_smiles, lig2_smiles, linker_file, base_dir):
    """处理单个 linker 文件
    
    Args:
        decoy_mol: decoy 分子对象
        linker_mol: linker 分子对象
        lig1_name: 配体1名称
        lig2_name: 配体2名称
        lig1_smiles: 配体1 SMILES
        lig2_smiles: 配体2 SMILES
        linker_file: linker PDB 文件路径
        base_dir: 基础输出目录
    
    Returns:
        tuple: (success, message)
    """
    
    # 1. 在 decoy 的配体中搜索 SMILES 匹配（找到重叠区域）
    lig1_in_decoy = find_substructure_in_residue(decoy_mol, lig1_name, lig1_smiles)
    lig2_in_decoy = find_substructure_in_residue(decoy_mol, lig2_name, lig2_smiles)
    
    if not lig1_in_decoy or not lig2_in_decoy:
        return False, "Decoy 配体中未找到 SMILES 匹配的子结构"
    
    # 2. 在 linker 中搜索匹配的原子
    lig1_count, lig1_in_linker = find_matching_atoms(linker_mol, lig1_smiles)
    lig2_count, lig2_in_linker = find_matching_atoms(linker_mol, lig2_smiles)
    
    # 检查匹配唯一性
    if lig1_count > 1:
        return False, f"配体1 SMILES 匹配到 {lig1_count} 个结构（不唯一）"
    if lig2_count > 1:
        return False, f"配体2 SMILES 匹配到 {lig2_count} 个结构（不唯一）"
    
    if not lig1_in_linker and not lig2_in_linker:
        return False, "在 linker 中未找到匹配的原子"
    
    # 3. 生成对齐原子列表（匹配到的原子用于刚体对齐）
    decoy_align = lig1_in_decoy + lig2_in_decoy
    linker_align = [('UNL', aname) for aname in lig1_in_linker] + \
                   [('UNL', aname) for aname in lig2_in_linker]
    
    # 4. 生成删除原子列表
    # Decoy 侧：匹配到的原子全部删除（因为 linker 已包含完整结构）
    decoy_delete = lig1_in_decoy + lig2_in_decoy
    
    # Linker 侧：匹配到的原子全部保留（不删除）
    linker_delete = []
    
    # 5. 创建子目录并写入文件
    linker_name = Path(linker_file).stem
    output_dir = base_dir / linker_name
    output_dir.mkdir(parents=True, exist_ok=True)
    
    with open(output_dir / 'decoy_atom_list.txt', 'w') as f:
        for resname, aname in decoy_align:
            f.write(f"{resname} {aname}\n")
    
    with open(output_dir / 'linker_atom_list.txt', 'w') as f:
        for resname, aname in linker_align:
            f.write(f"{resname} {aname}\n")
    
    with open(output_dir / 'decoy_atom_list_delete.txt', 'w') as f:
        for resname, aname in decoy_delete:
            f.write(f"{resname} {aname}\n")
    
    with open(output_dir / 'linker_atom_list_delete.txt', 'w') as f:
        for resname, aname in linker_delete:
            f.write(f"{resname} {aname}\n")
    
    return True, f"成功处理 {linker_name}"


def print_ligand_smiles(pdb_file, resnames):
    """根据残基名提取小分子并生成 SMILES
    
    Args:
        pdb_file: PDB 文件路径
        resnames: 残基名列表，如 ['KIN', 'CBN']
    """
    # 使用 proximityBonding=True 以获取键信息
    mol = Chem.MolFromPDBFile(
        pdb_file, sanitize=False, removeHs=True, flavor=0, proximityBonding=True)
    
    if mol is None:
        print("无法读取 PDB 文件以生成 SMILES")
        return
    
    print("\nDecoy 中的小分子 SMILES:")
    print("-" * 60)
    
    for resname in resnames:
        # 找到目标残基的原子索引
        target_indices = set()
        for i, atom in enumerate(mol.GetAtoms()):
            info = atom.GetPDBResidueInfo()
            if info and info.GetResidueName().strip() == resname:
                target_indices.add(i)
        
        if not target_indices:
            print(f"  {resname}: 未找到该残基")
            continue
        
        # 用 RWMol 删除非目标原子
        rw = Chem.RWMol(mol)
        non_target = [i for i in range(mol.GetNumAtoms()) if i not in target_indices]
        for idx in sorted(non_target, reverse=True):
            rw.RemoveAtom(idx)
        
        fragment = rw.GetMol()
        try:
            Chem.SanitizeMol(fragment)
        except:
            pass
        smiles = Chem.MolToSmiles(fragment)
        print(f"  {resname}: {smiles}")
    
    print("-" * 60)


def main():
    parser = argparse.ArgumentParser(
        description="批量自动生成 PROTAC 三元复合物建模所需的原子列表文件")
    
    parser.add_argument('--decoy', required=True, help='Decoy PDB 文件路径')
    parser.add_argument('--linkers', required=True, help='Linker PDB 文件路径或目录')
    parser.add_argument('--lig1-smiles', required=True, help='配体1的 SMILES')
    parser.add_argument('--lig1-name', required=True, help='配体1的残基名')
    parser.add_argument('--lig2-smiles', required=True, help='配体2的 SMILES')
    parser.add_argument('--lig2-name', required=True, help='配体2的残基名')
    parser.add_argument('--output-dir', default='.', help='输出基础目录（默认为当前目录）')
    
    args = parser.parse_args()
    
    # 设置基础输出目录
    base_dir = Path(args.output_dir)
    base_dir.mkdir(parents=True, exist_ok=True)
    
    # 读取 decoy
    print(f"读取 decoy: {args.decoy}")
    decoy_mol = Chem.rdmolfiles.MolFromPDBFile(
        args.decoy, sanitize=False, removeHs=False, flavor=0, proximityBonding=False)
    if decoy_mol is None:
        print("错误：无法读取 decoy PDB")
        sys.exit(1)
    
    # 打印 decoy 中小分子的 SMILES
    print_ligand_smiles(args.decoy)
    
    # 获取 linker 文件列表
    linker_path = Path(args.linkers)
    if linker_path.is_dir():
        linker_files = sorted(glob.glob(str(linker_path / '*.pdb')))
    else:
        linker_files = [args.linkers]
    
    print(f"找到 {len(linker_files)} 个 linker 文件")
    print(f"配体1: {args.lig1_name} ({args.lig1_smiles})")
    print(f"配体2: {args.lig2_name} ({args.lig2_smiles})")
    print(f"输出目录: {base_dir}\n")
    
    # 批量处理
    success_count = 0
    fail_count = 0
    
    for linker_file in linker_files:
        print(f"处理: {Path(linker_file).name}", end=' ')
        
        linker_mol = Chem.rdmolfiles.MolFromPDBFile(
            linker_file, sanitize=False, removeHs=False, flavor=0, proximityBonding=False)
        
        if linker_mol is None:
            print("❌ 无法读取")
            fail_count += 1
            continue
        
        success, msg = process_single_linker(
            decoy_mol, linker_mol, 
            args.lig1_name, args.lig2_name,
            args.lig1_smiles, args.lig2_smiles,
            linker_file, base_dir)
        
        if success:
            print("✓")
            success_count += 1
        else:
            print(f"❌ {msg}")
            fail_count += 1
    
    print(f"\n处理完成: 成功 {success_count}, 失败 {fail_count}")


if __name__ == '__main__':
    main()
