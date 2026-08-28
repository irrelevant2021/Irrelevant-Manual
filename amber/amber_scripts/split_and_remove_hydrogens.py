#!/usr/bin/env python
"""
PDB 拆分与去氢脚本
功能:
  将三元复合物 PDB 拆分为:
    - protein.pdb  (蛋白，已去除所有氢原子)
    - ligand.pdb   (配体，非蛋白非水分子)

使用方法: python split_and_remove_hydrogens.py input.pdb
"""

import sys
import MDAnalysis as mda


def main():
    if len(sys.argv) != 2:
        print("用法: python split_and_remove_hydrogens.py input.pdb")
        print("示例: python split_and_remove_hydrogens.py mini_ternary20_mod_loop.pdb")
        sys.exit(1)

    input_pdb = sys.argv[1]
    print(f"加载 PDB 文件: {input_pdb}")
    u = mda.Universe(input_pdb)

    # 蛋白：标准氨基酸，去除氢原子
    protein = u.select_atoms('protein and not element H')
    # 配体：非蛋白、非水
    ligand = u.select_atoms('not protein and not resname WAT HOH')

    print(f"蛋白重原子数: {protein.n_atoms}")
    print(f"配体原子数:   {ligand.n_atoms}")
    print(f"配体残基:     {sorted(set(ligand.residues.resnames))}")

    protein.write('protein_noH.pdb')
    ligand.write('ligand.pdb')

    print("已保存: protein_noH.pdb (去氢)")
    print("已保存: ligand.pdb")


if __name__ == "__main__":
    main()
