# 实现总结

## 项目结构

```
PROTAC_ternary-master/
├── generate_atom_lists.py              # 主脚本（项目根目录）
├── GENERATE_ATOM_LISTS_README.md       # 使用说明（项目根目录）
├── ternary_model_prediction.py         # 原有的三元复合物预测脚本
├── PROTAC_ternary/                     # 原有的核心模块
└── Example1_ternary_model_prediction/  # 示例工作目录
    ├── docking_decoy.pdb               # Decoy 文件
    ├── linkers/                        # Linker PDB 文件目录
    │   ├── linker_001.pdb
    │   ├── linker_002.pdb
    │   └── linker_003.pdb
    ├── linker_001/                     # 生成的输出目录
    │   ├── decoy_atom_list.txt
    │   ├── linker_atom_list.txt
    │   ├── decoy_atom_list_delete.txt
    │   └── linker_atom_list_delete.txt
    ├── linker_002/
    │   └── ... (4 个文件)
    └── linker_003/
        └── ... (4 个文件)
```

## 功能特性

### ✅ 已实现
1. **批量处理**：支持处理整个目录下的所有 linker PDB 文件
2. **自动匹配**：根据 SMILES 在 linker 中自动搜索配体子结构
3. **独立子目录**：每个 linker 在自己的子目录中生成文件
4. **错误处理**：无法读取或匹配失败时给出明确提示
5. **统计输出**：显示成功/失败的处理数量

### 📁 输出结构
- 脚本文件位于项目根目录
- 输出文件位于执行脚本的工作目录
- 每个 linker PDB 对应一个独立的子目录
- 每个子目录包含 4 个标准命名的 txt 文件

## 使用示例

### 1. 准备工作目录
```bash
cd Example1_ternary_model_prediction/

# 确保有 linker PDB 文件
mkdir -p linkers/
cp linker_conformer.pdb linkers/linker_001.pdb
cp linker_conformer.pdb linkers/linker_002.pdb
cp linker_conformer.pdb linkers/linker_003.pdb
```

### 2. 运行脚本
```bash
python ../generate_atom_lists.py \
    --decoy docking_decoy.pdb \
    --linkers linkers/ \
    --lig1-smiles "CCO" \
    --lig1-name KIN \
    --lig2-smiles "c1ccccc1" \
    --lig2-name CBN
```

### 3. 查看输出
```bash
ls -d linker_*/
# 输出: linker_001/  linker_002/  linker_003/

ls linker_001/
# 输出: decoy_atom_list.txt  linker_atom_list.txt
#       decoy_atom_list_delete.txt  linker_atom_list_delete.txt
```

### 4. 用于三元复合物建模
```bash
cd linker_001/

python ../../ternary_model_prediction.py \
    -d ../docking_decoy.pdb \
    -l ../linkers/linker_001.pdb \
    -da decoy_atom_list.txt \
    -la linker_atom_list.txt \
    -wd decoy_atom_list_delete.txt \
    -ld linker_atom_list_delete.txt \
    -t default \
    -r rmsd.txt
```

## 测试结果

```
读取 decoy: docking_decoy.pdb
找到 3 个 linker 文件
配体1: KIN (CCO)
配体2: CBN (c1ccccc1)
输出目录: .

处理: linker_001.pdb ✓
处理: linker_002.pdb ✓
处理: linker_003.pdb ✓

处理完成: 成功 3, 失败 0
```

生成的目录结构：
- linker_001/: 4 个文件
- linker_002/: 4 个文件
- linker_003/: 4 个文件

## 技术细节

### 核心算法
1. **子结构匹配**：使用 RDKit 的 `GetSubstructMatches()` 搜索配体在 linker 中的位置
2. **对齐原子选择**：
   - 配体1：前 5 个重原子
   - 配体2：前 3 个重原子
3. **删除原子选择**：
   - Decoy 侧：两个配体的所有重原子
   - Linker 侧：匹配到的原子及其相连的 H 原子

### 依赖
- Python 3.x
- RDKit

## 文件清单

### 项目根目录
- `generate_atom_lists.py` (7.6 KB) - 主脚本
- `GENERATE_ATOM_LISTS_README.md` (3.3 KB) - 使用说明

### 工作目录（示例）
- `Example1_ternary_model_prediction/linker_001/` - 输出示例
- `Example1_ternary_model_prediction/linkers/` - Linker 输入文件

## 下一步

1. 获取实际配体的 SMILES 进行真实数据测试
2. 根据实际需求调整对齐原子数量
3. 考虑添加日志文件记录详细处理信息
