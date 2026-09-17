#!/usr/bin/env bash
# ============================================================
# single_extract_chain.sh
# 单文件版 extract_chain 逻辑
#
# 用法:
#   ./single_extract_chain.sh <input_pdb> <output_dir> [chain_id]
#
# 环境变量:
#   MISSING_RESIDUES - 缺失残基信息 (Python dict 格式)
#     默认: "{(0, 210): ['SER', 'VAL', 'ASP', 'PHE', 'PRO', 'GLU']}"
#
# 输出:
#   <output_dir>/{base}_mod_loop.pdb
#
# 环境依赖:
#   - base (PyMOL)
#   - openfe_env (pdbfixer)
# ============================================================
set -e

INPUT_PDB="$1"
OUTPUT_DIR="$2"
CHAIN_ID="${3:-C}"

if [ -z "$INPUT_PDB" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "用法: $0 <input_pdb> <output_dir> [chain_id]"
    exit 1
fi

if [ ! -f "$INPUT_PDB" ]; then
    echo "错误: 输入文件不存在: $INPUT_PDB"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# 将 OUTPUT_DIR 转换为绝对路径 (在 cd 到工作目录之前)
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"

# 获取基础名称 (去掉 _refine_mod.pdb 后缀)
BASE_NAME=$(basename "$INPUT_PDB" _refine_mod.pdb)
if [ "$BASE_NAME" = "$(basename "$INPUT_PDB")" ]; then
    BASE_NAME=$(basename "$INPUT_PDB" .pdb)
fi

# 工作目录
WORK_DIR=$(mktemp -d)
trap "rm -rf $WORK_DIR" EXIT

# 复制输入文件到工作目录
cp "$INPUT_PDB" "$WORK_DIR/input.pdb"
cd "$WORK_DIR"

echo "  [阶段3] 处理: $BASE_NAME (链: $CHAIN_ID)"

# 步骤1: 使用 PyMOL 提取链并重编号
CHAIN_C_PDB="${BASE_NAME}_chainC.pdb"
mamba run -n base pymol -cq -d "
cmd.load('input.pdb', 'protein')
cmd.extract('chain_${CHAIN_ID}', 'chain ${CHAIN_ID}')
stored.residues = []
cmd.iterate('chain_${CHAIN_ID}', 'stored.residues.append(resi)')
residue_map = {}
for i, old_resi in enumerate(sorted(set(stored.residues), key=lambda x: int(x) if x.lstrip('-').isdigit() else 0)):
    residue_map[old_resi] = str(i + 1)
for old_resi, new_resi in residue_map.items():
    cmd.alter(f'chain_${CHAIN_ID} and resi {old_resi}', f'resi=\"{new_resi}\"')
cmd.sort('chain_${CHAIN_ID}')
cmd.save('${CHAIN_C_PDB}', 'chain_${CHAIN_ID}')
cmd.quit()
" 2>/dev/null

if [ ! -f "$CHAIN_C_PDB" ]; then
    echo "  错误: PyMOL 提取链失败"
    exit 1
fi

echo "  [阶段3] 链提取完成: $CHAIN_C_PDB"

# 步骤2: 使用 pdbfixer 补全缺失残基
LOOP_PDB="${BASE_NAME}_chainC_loop.pdb"

# 设置缺失残基信息 (从环境变量读取，默认值与原始脚本一致)
if [ -z "$MISSING_RESIDUES" ]; then
    MISSING_RESIDUES="{(0, 217): ['SER', 'VAL', 'ASP', 'PHE', 'PRO', 'GLU']}"
fi

mamba run -n openfe_env python << PYTHON_SCRIPT
from openmm.app import PDBFile
from pdbfixer import PDBFixer

fixer = PDBFixer(filename='${CHAIN_C_PDB}')
fixer.findMissingResidues()
# 强制设置缺失残基 (从环境变量读取)
missing_residues = $MISSING_RESIDUES
fixer.missingResidues = missing_residues
fixer.findMissingAtoms()
fixer.addMissingAtoms()
PDBFile.writeFile(fixer.topology, fixer.positions, open('${LOOP_PDB}', 'w'))
print("  补全缺失残基完成")
PYTHON_SCRIPT

if [ ! -f "$LOOP_PDB" ]; then
    echo "  错误: pdbfixer 补全失败"
    exit 1
fi

# 步骤3: 使用 PyMOL 将补全的链放回原结构
MOD_LOOP_PDB="${BASE_NAME}_mod_loop.pdb"
mamba run -n base pymol -cq -d "
cmd.load('input.pdb', 'original')
cmd.load('${LOOP_PDB}', 'new_chain')
# 将新链的链标识符改为原来的链标识符
cmd.alter('new_chain', 'chain=\"${CHAIN_ID}\"')
# 删除原始结构中的旧链
cmd.remove('original and chain ${CHAIN_ID}')
# 将新链合并到original
cmd.create('merged', 'original or new_chain')
# 保存
cmd.save('${MOD_LOOP_PDB}', 'merged')
cmd.quit()
" 2>/dev/null

if [ ! -f "$MOD_LOOP_PDB" ]; then
    echo "  错误: PyMOL 合并失败"
    exit 1
fi

# 复制输出文件
OUTPUT_FILE="$OUTPUT_DIR/${BASE_NAME}_mod_loop.pdb"
cp "$MOD_LOOP_PDB" "$OUTPUT_FILE"

echo "  [阶段3] 完成: $OUTPUT_FILE"
