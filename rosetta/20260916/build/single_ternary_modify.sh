#!/usr/bin/env bash
# ============================================================
# single_ternary_modify.sh
# 单文件版 ternary_modify 逻辑
#
# 用法:
#   ./single_ternary_modify.sh <input_pdb> <output_dir>
#
# 输入:
#   <input_pdb>   - *_refined.pdb 文件
#   <output_dir>  - 输出目录
#
# 输出:
#   <output_dir>/{base}_refine_mod.pdb
#   <output_dir>/{base}_refined.params
#   <output_dir>/{base}_refine_lig.mol2
#   <output_dir>/{base}_refine_lig.pdb
#
# 环境依赖:
#   - BABEL (环境变量)
#   - MOL2PARAMS (环境变量)
#   - py26 conda 环境
# ============================================================
set -e

INPUT_PDB="$1"
OUTPUT_DIR="$2"

if [ -z "$INPUT_PDB" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "用法: $0 <input_pdb> <output_dir>"
    exit 1
fi

if [ ! -f "$INPUT_PDB" ]; then
    echo "错误: 输入文件不存在: $INPUT_PDB"
    exit 1
fi

if [ -z "$BABEL" ]; then
    echo "错误: 环境变量 BABEL 未设置"
    exit 1
fi

if [ -z "$MOL2PARAMS" ]; then
    echo "错误: 环境变量 MOL2PARAMS 未设置"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# 将 OUTPUT_DIR 转换为绝对路径 (在 cd 到工作目录之前)
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"

# 获取基础名称 (去掉 _refined.pdb 后缀)
BASE_NAME=$(basename "$INPUT_PDB" _refined.pdb)
# 如果文件名不包含 _refined.pdb，则去掉 .pdb
if [ "$BASE_NAME" = "$(basename "$INPUT_PDB")" ]; then
    BASE_NAME=$(basename "$INPUT_PDB" .pdb)
fi

# 工作目录
WORK_DIR=$(mktemp -d)
trap "rm -rf $WORK_DIR" EXIT

# 复制输入文件到工作目录
cp "$INPUT_PDB" "$WORK_DIR/ternary.pdb"
cd "$WORK_DIR"

echo "  [阶段2] 处理: $BASE_NAME"

# 步骤1: 提取 HETATM 到单独的 PDB
python3 << 'PYTHON_SCRIPT'
import collections
import os

def lig_to_dictionary(ter_file):
    lig_dic = collections.OrderedDict()
    with open(ter_file, "r") as f:
        counter = 0
        for line in f:
            if line.startswith("HETATM"):
                key = counter
                val = line.rstrip('\n')
                lig_dic[key] = val
                counter += 1
    return lig_dic

def lig_pdb_maker(ter_file, lig_dic):
    name = ter_file.strip(".pdb") + "_lig.pdb"
    with open(name, "w") as f:
        for value in lig_dic.values():
            f.write(value + '\n')
    print(f"  {name} 创建成功")
    return name

# 主流程
orig_lig = lig_to_dictionary("ternary.pdb")
lig_pdb = lig_pdb_maker("ternary.pdb", orig_lig)
print("  配体提取完成")
PYTHON_SCRIPT

# 步骤2: PDB -> MOL2 (使用 babel)
echo "  [阶段2] 转换 PDB -> MOL2"
LIG_PDB="ternary_lig.pdb"
LIG_MOL2="ternary_lig.mol2"
"$BABEL" -h -ipdb "$LIG_PDB" -omol2 -O "$LIG_MOL2" 2>/dev/null || {
    echo "  警告: babel 转换失败，尝试不使用 -h 参数"
    "$BABEL" -ipdb "$LIG_PDB" -omol2 -O "$LIG_MOL2"
}

# 步骤3: MOL2 -> params (使用 py26 环境)
echo "  [阶段2] 生成 params 文件"
mamba run -n py26 python2.6 "$MOL2PARAMS" "$LIG_MOL2" -n TRN --clobber -p ternary 2>/dev/null || {
    echo "  警告: py26 环境失败，尝试使用默认 python"
    python "$MOL2PARAMS" "$LIG_MOL2" -n TRN --clobber -p ternary
}

# 步骤4: 插入新配体坐标
echo "  [阶段2] 插入新配体坐标"
python3 << 'PYTHON_SCRIPT'
import collections
import os
import sys

def insert_lig(pdb2_lig, ter_file):
    lig2_dic = collections.OrderedDict()
    with open(pdb2_lig, "r") as f:
        counter = 0
        for line in f:
            lig2_dic[counter] = line
            counter += 1
    # 去掉最后一行 (TER)
    if counter > 0:
        lig2_dic.pop(counter - 1)
    
    mod_pdb = ter_file.strip(".pdb") + "_mod.pdb"
    
    # 写入蛋白质部分
    with open(mod_pdb, 'w') as ff:
        with open(ter_file, "r") as f:
            judgement = False
            for line in f:
                if line.startswith("ATOM"):
                    judgement = True
                if line.startswith("HETATM"):
                    judgement = False
                if line.startswith("TER"):
                    judgement = False
                elif judgement:
                    ff.write(line)
    
    # 追加新的配体坐标
    with open(mod_pdb, 'a') as fff:
        for key in lig2_dic:
            fff.write(lig2_dic[key])
    
    # 验证文件是否生成
    if os.path.exists(mod_pdb):
        print(f"  最终 PDB: {mod_pdb} (大小: {os.path.getsize(mod_pdb)} bytes)")
        return mod_pdb
    else:
        print(f"  错误: {mod_pdb} 未生成")
        sys.exit(1)

# 读取新生成的配体 PDB (从 params 生成)
new_lig_pdb = "ternary_0001.pdb"
if not os.path.exists(new_lig_pdb):
    print(f"  错误: {new_lig_pdb} 不存在")
    sys.exit(1)

insert_lig(new_lig_pdb, "ternary.pdb")
PYTHON_SCRIPT

# 步骤5: 复制输出文件
OUTPUT_BASE="$OUTPUT_DIR/$BASE_NAME"
cp "ternary_mod.pdb" "${OUTPUT_BASE}_refine_mod.pdb" 2>/dev/null || {
    echo "  错误: ternary_mod.pdb 未生成"
    exit 1
}

# 复制 params 文件 (如果存在)
if [ -f "ternary.params" ]; then
    cp "ternary.params" "${OUTPUT_BASE}_refined.params"
fi

# 复制 mol2 文件
cp "$LIG_MOL2" "${OUTPUT_BASE}_refine_lig.mol2" 2>/dev/null || true
cp "$LIG_PDB" "${OUTPUT_BASE}_refine_lig.pdb" 2>/dev/null || true

echo "  [阶段2] 完成: ${OUTPUT_BASE}_refine_mod.pdb"
