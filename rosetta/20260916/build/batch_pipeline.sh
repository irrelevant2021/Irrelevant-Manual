#!/usr/bin/env bash
# ============================================================
# batch_pipeline.sh
# PROTAC Ternary 批量调度脚本 — 完整 4 阶段流水线
#
# 用法:
#   ./batch_pipeline.sh                       # 默认 64 并发
#   ./batch_pipeline.sh 32                    # 32 并发
#   ./batch_pipeline.sh 64 ../ppd ../linkers  # 指定目录
#   ./batch_pipeline.sh 64 ../ppd ../linkers ../outputs
#
# 完整流水线 (4 阶段):
#   1. pipeline.py      - Schrödinger + RDKit 配体优化
#   2. ternary_modify   - 配体坐标替换 + params 生成
#   3. extract_chain    - PyMOL 链提取 + pdbfixer 补全
#   4. minimize_ppi     - Rosetta MPI 能量最小化
#
# 输出结构:
#   ../outputs/{MAE_NAME}/
#     ├── {PDB}_refined.pdb          # 阶段1
#     ├── {PDB}_refine_mod.pdb       # 阶段2
#     ├── {PDB}_refined.params       # 阶段2
#     ├── loop/
#     │   └── {PDB}_mod_loop.pdb     # 阶段3
#     └── mini/
#         ├── mini_{PDB}_mod_loop.pdb # 阶段4
#         └── {PDB}_score.sc         # 阶段4
#
# 特性:
#   - 恒定 N 并发 (后台进程池)
#   - 输出已存在时自动跳过 (可恢复)
#   - 单任务/单阶段失败不影响其他
#   - 实时进度 (每轮/每秒)
#   - 最终汇总表
# ============================================================
set -o pipefail

# ─── 参数 ───────────────────────────────────────────────
MAX_JOBS="${1:-64}"
PPD_DIR="${2:-../ppd}"
MAE_DIR="${3:-../linkers}"
OUTPUT_ROOT="${4:-../outputs}"

# ─── 环境变量配置 ───────────────────────────────────────
export BABEL="${BABEL:-/home/xfusion/miniforge3/envs/maize/bin/obabel}"
export MOL2PARAMS="${MOL2PARAMS:-/home/xfusion/rosetta/rosetta.source.release-408/main/source/scripts/python/public/molfile_to_params.py}"
export ROSETTA="${ROSETTA:-/home/xfusion/rosetta/rosetta.source.release-408}"
export ROSETTA_BIN="$ROSETTA/main/source/bin"

# conda 环境名
ENV_OPENFE="openfe_env"
ENV_PY26="py26"
ENV_BASE="base"

# ─── 路径固化 ───────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PPD_DIR="$(cd "$PPD_DIR" 2>/dev/null && pwd)" || {
    echo "ERROR: PPD_DIR ($PPD_DIR) 不存在"; exit 1
}
MAE_DIR="$(cd "$MAE_DIR" 2>/dev/null && pwd)" || {
    echo "ERROR: MAE_DIR ($MAE_DIR) 不存在"; exit 1
}
OUTPUT_ROOT="$(realpath -m "$OUTPUT_ROOT")"
PIPELINE="$SCRIPT_DIR/pipeline.py"
TERNARY_MODIFY="$SCRIPT_DIR/single_ternary_modify.sh"
EXTRACT_CHAIN="$SCRIPT_DIR/single_extract_chain.sh"
LOG_DIR="$OUTPUT_ROOT/.batch_logs"
STATUS_DIR="$SCRIPT_DIR/.batch_status_$$"
START_TIME=$(date +%s)

mkdir -p "$LOG_DIR" "$OUTPUT_ROOT" "$STATUS_DIR"

# ─── 清理 ───────────────────────────────────────────────
cleanup() {
    echo ""
    echo "  [WARN] 收到中断，正在终止子进程..."
    local pids
    pids=$(jobs -rp 2>/dev/null)
    [ -n "$pids" ] && kill $pids 2>/dev/null || true
    wait 2>/dev/null || true
    echo "  [WARN] 全部终止"
    summarize
    rm -rf "$STATUS_DIR"
    exit 1
}
trap cleanup SIGINT SIGTERM
trap 'rm -rf "$STATUS_DIR"' EXIT

# ─── 工具函数 ───────────────────────────────────────────
log() { echo "  [$(date '+%H:%M:%S')] $*"; }

summarize() {
    local ok fail skip
    ok=$(find "$STATUS_DIR" -name '*.OK'   2>/dev/null | wc -l)
    fail=$(find "$STATUS_DIR" -name '*.FAIL' 2>/dev/null | wc -l)
    skip=$(find "$STATUS_DIR" -name '*.SKIP' 2>/dev/null | wc -l)
    
    # 各阶段失败统计
    local s1_fail s2_fail s3_fail s4_fail s4_skip
    s1_fail=$(find "$STATUS_DIR" -name '*.STAGE1_FAIL' 2>/dev/null | wc -l)
    s2_fail=$(find "$STATUS_DIR" -name '*.STAGE2_FAIL' 2>/dev/null | wc -l)
    s3_fail=$(find "$STATUS_DIR" -name '*.STAGE3_FAIL' 2>/dev/null | wc -l)
    s4_fail=$(find "$STATUS_DIR" -name '*.STAGE4_FAIL' 2>/dev/null | wc -l)
    s4_skip=$(find "$STATUS_DIR" -name '*.STAGE4_SKIP' 2>/dev/null | wc -l)
    
    local elapsed=$(( $(date +%s) - START_TIME ))
    local total=$(( ok + fail + skip + s1_fail + s2_fail + s3_fail + s4_fail + s4_skip ))
    echo ""
    echo "  ═══════════════════════════════════════════"
    echo "   汇总:  成功 ${ok}   |   跳过 ${skip}   |   失败 ${fail}"
    echo "  ───────────────────────────────────────────"
    echo "   阶段1失败: ${s1_fail}   |   阶段2失败: ${s2_fail}"
    echo "   阶段3失败: ${s3_fail}   |   阶段4失败: ${s4_fail}"
    echo "   阶段4跳过: ${s4_skip} (无params)"
    echo "  ───────────────────────────────────────────"
    echo "   总计:  ${total} 个任务  |  耗时: ${elapsed}s"
    echo "  ═══════════════════════════════════════════"
    echo ""
    
    # 显示各阶段失败详情
    local has_failures=0
    if [ "$s1_fail" -gt 0 ]; then
        echo "  ❌ 阶段1失败任务:"
        find "$STATUS_DIR" -name '*.STAGE1_FAIL' -printf '    - %f\n' | sed 's/\.STAGE1_FAIL$//'
        has_failures=1
    fi
    if [ "$s2_fail" -gt 0 ]; then
        echo "  ❌ 阶段2失败任务:"
        find "$STATUS_DIR" -name '*.STAGE2_FAIL' -printf '    - %f\n' | sed 's/\.STAGE2_FAIL$//'
        has_failures=1
    fi
    if [ "$s3_fail" -gt 0 ]; then
        echo "  ❌ 阶段3失败任务:"
        find "$STATUS_DIR" -name '*.STAGE3_FAIL' -printf '    - %f\n' | sed 's/\.STAGE3_FAIL$//'
        has_failures=1
    fi
    if [ "$s4_fail" -gt 0 ]; then
        echo "  ❌ 阶段4失败任务:"
        find "$STATUS_DIR" -name '*.STAGE4_FAIL' -printf '    - %f\n' | sed 's/\.STAGE4_FAIL$//'
        has_failures=1
    fi
    if [ "$has_failures" -eq 1 ]; then
        echo ""
    fi
}

# ─── 提取关键参数到 CSV ─────────────────────────────────
extract_summary() {
    local csv_file="$OUTPUT_ROOT/summary.csv"
    local count=0
    
    echo ""
    echo "  📊 提取关键参数到 CSV..."
    
    # 写入 CSV 头部
    echo "MAE,PDB,Final_Score,Backbone_RMS,Ligand_RMSD,Best_Energy_kJ" > "$csv_file"
    
    # 遍历所有 MAE 目录
    for mae_dir in "$OUTPUT_ROOT"/*/; do
        [ -d "$mae_dir" ] || continue
        [ "$(basename "$mae_dir")" = ".batch_logs" ] && continue
        
        local mae_name
        mae_name=$(basename "$mae_dir")
        
        # 遍历所有 mini 目录中的 PDB 文件
        for mini_pdb in "$mae_dir"mini/mini_*_mod_loop.pdb; do
            [ -f "$mini_pdb" ] || continue
            
            local pdb_name
            pdb_name=$(basename "$mini_pdb" _mod_loop.pdb)
            pdb_name=${pdb_name#mini_}  # 去掉 mini_ 前缀
            
            local log_file="$LOG_DIR/$mae_name/${pdb_name}.log"
            local energies_file="$mae_dir/${pdb_name}_refined.energies.txt"
            
            # 提取参数
            local final_score="" backbone_rms="" ligand_rmsd="" best_energy=""
            
            if [ -f "$log_file" ]; then
                # Final score
                final_score=$(grep -oP "Final score: \K[-0-9.]+" "$log_file" 2>/dev/null | tail -1)
                
                # Backbone RMS to reference
                backbone_rms=$(grep -oP "Backbone RMS to reference: \K[0-9.]+" "$log_file" 2>/dev/null | tail -1)
                
                # Ligand RMSD
                ligand_rmsd=$(grep -oP "ligand RMSD for.*:\K[0-9.]+" "$log_file" 2>/dev/null | tail -1)
            fi
            
            if [ -f "$energies_file" ]; then
                # Best energy
                best_energy=$(grep -oP "# Best: \K[-0-9.]+" "$energies_file" 2>/dev/null)
            fi
            
            # 写入 CSV (使用空字符串填充缺失值)
            echo "${mae_name},${pdb_name},${final_score:-},${backbone_rms:-},${ligand_rmsd:-},${best_energy:-}" >> "$csv_file"
            count=$((count + 1))
        done
    done
    
    echo "  ✅ 已生成: $csv_file ($count 条记录)"
    echo ""
}

run_one() {
    local pdb="$1"
    local mae="$2"
    local out_base="$3"  # 输出基础路径 (不含扩展名)
    local pdb_base; pdb_base="$(basename "$pdb" .pdb)"
    local mae_base; mae_base="$(basename "$mae" .mae)"
    local task_label="${pdb_base}__${mae_base}"
    local task_log="$LOG_DIR/${mae_base}/${pdb_base}.log"

    mkdir -p "$(dirname "$task_log")"
    mkdir -p "$out_base/loop" "$out_base/mini"

    # 定义各阶段输出文件
    local stage1_out="${out_base}/${pdb_base}_refined.pdb"
    local stage2_out="${out_base}/${pdb_base}_refine_mod.pdb"
    local stage2_params="${out_base}/${pdb_base}_refined.params"
    local stage3_out="${out_base}/loop/${pdb_base}_mod_loop.pdb"
    local stage4_out="${out_base}/mini/mini_${pdb_base}_mod_loop.pdb"
    local stage4_score="${out_base}/mini/${pdb_base}_score.sc"

    # 检查最终输出是否已存在 (可恢复)
    if [ -s "$stage4_out" ]; then
        touch "$STATUS_DIR/${task_label}.SKIP"
        return 0
    fi

    # ═══════════════════════════════════════════════════════════════
    # 阶段 1: pipeline.py (Schrödinger + RDKit)
    # ═══════════════════════════════════════════════════════════════
    if [ ! -s "$stage1_out" ]; then
        echo "[$(date '+%H:%M:%S')] [阶段1] 开始: $task_label" >> "$task_log"
        if ! mamba run -n "$ENV_OPENFE" python3 "$PIPELINE" \
            --pdb "$pdb" \
            --mae "$mae" \
            -o "$stage1_out" \
            >> "$task_log" 2>&1
        then
            echo "[$(date '+%H:%M:%S')] [阶段1] 失败: $task_label" >> "$task_log"
            touch "$STATUS_DIR/${task_label}.STAGE1_FAIL"
            echo "    [STAGE1_FAIL] ${task_label}"
            return 1
        fi
        echo "[$(date '+%H:%M:%S')] [阶段1] 完成: $task_label" >> "$task_log"
    else
        echo "[$(date '+%H:%M:%S')] [阶段1] 跳过 (已存在): $task_label" >> "$task_log"
    fi

    # ═══════════════════════════════════════════════════════════════
    # 阶段 2: ternary_modify (配体坐标替换 + params)
    # ═══════════════════════════════════════════════════════════════
    if [ ! -s "$stage2_out" ]; then
        echo "[$(date '+%H:%M:%S')] [阶段2] 开始: $task_label" >> "$task_log"
        if ! bash "$TERNARY_MODIFY" "$stage1_out" "$out_base" \
            >> "$task_log" 2>&1
        then
            echo "[$(date '+%H:%M:%S')] [阶段2] 失败: $task_label" >> "$task_log"
            touch "$STATUS_DIR/${task_label}.STAGE2_FAIL"
            echo "    [STAGE2_FAIL] ${task_label}"
            return 1
        fi
        echo "[$(date '+%H:%M:%S')] [阶段2] 完成: $task_label" >> "$task_log"
    else
        echo "[$(date '+%H:%M:%S')] [阶段2] 跳过 (已存在): $task_label" >> "$task_log"
    fi

    # 检查 params 文件
    if [ ! -s "$stage2_params" ]; then
        echo "[$(date '+%H:%M:%S')] [阶段2] 警告: params 文件未生成: $stage2_params" >> "$task_log"
    fi

    # ═══════════════════════════════════════════════════════════════
    # 阶段 3: extract_chain (PyMOL 链提取 + pdbfixer)
    # ═══════════════════════════════════════════════════════════════
    if [ ! -s "$stage3_out" ]; then
        echo "[$(date '+%H:%M:%S')] [阶段3] 开始: $task_label" >> "$task_log"
        
        # 设置缺失残基信息 (从环境变量读取，默认值与原始脚本一致)
        export MISSING_RESIDUES="${MISSING_RESIDUES:-{(0, 217): ['SER', 'VAL', 'ASP', 'PHE', 'PRO', 'GLU']}}"
        
        if ! bash "$EXTRACT_CHAIN" "$stage2_out" "${out_base}/loop" "C" \
            >> "$task_log" 2>&1
        then
            echo "[$(date '+%H:%M:%S')] [阶段3] 失败: $task_label" >> "$task_log"
            touch "$STATUS_DIR/${task_label}.STAGE3_FAIL"
            echo "    [STAGE3_FAIL] ${task_label}"
            return 1
        fi
        echo "[$(date '+%H:%M:%S')] [阶段3] 完成: $task_label" >> "$task_log"
    else
        echo "[$(date '+%H:%M:%S')] [阶段3] 跳过 (已存在): $task_label" >> "$task_log"
    fi

    # ═══════════════════════════════════════════════════════════════
    # 阶段 4: minimize_ppi (Rosetta MPI)
    # ═══════════════════════════════════════════════════════════════
    if [ ! -s "$stage4_out" ]; then
        echo "[$(date '+%H:%M:%S')] [阶段4] 开始: $task_label" >> "$task_log"
        
        # 检查 params 文件是否存在
        if [ ! -s "$stage2_params" ]; then
            echo "[$(date '+%H:%M:%S')] [阶段4] 跳过: params 文件不存在" >> "$task_log"
            touch "$STATUS_DIR/${task_label}.STAGE4_SKIP"
            echo "    [STAGE4_SKIP] ${task_label} (无 params)"
            return 1
        fi
        
        # 在 mini 目录中运行 Rosetta (输出到 CWD)
        # minimize_ppi 输出: mini_{input_basename}.pdb
        if ! (
            cd "${out_base}/mini"
            "$ROSETTA_BIN/minimize_ppi.mpi.linuxgccrelease" \
                -s "../loop/${pdb_base}_mod_loop.pdb" \
                -extra_res_fa "../${pdb_base}_refined.params" \
                -jump_all \
                -out:file:scorefile "${pdb_base}_score.sc"
        ) >> "$task_log" 2>&1
        then
            echo "[$(date '+%H:%M:%S')] [阶段4] 失败: $task_label" >> "$task_log"
            touch "$STATUS_DIR/${task_label}.STAGE4_FAIL"
            echo "    [STAGE4_FAIL] ${task_label}"
            return 1
        fi
        
        # 检查输出文件 (minimize_ppi 输出: mini_{input_basename}.pdb)
        if [ ! -s "$stage4_out" ]; then
            echo "[$(date '+%H:%M:%S')] [阶段4] 警告: 输出文件未生成" >> "$task_log"
            touch "$STATUS_DIR/${task_label}.STAGE4_FAIL"
            echo "    [STAGE4_FAIL] ${task_label} (输出未生成)"
            return 1
        fi
        
        echo "[$(date '+%H:%M:%S')] [阶段4] 完成: $task_label" >> "$task_log"
    else
        echo "[$(date '+%H:%M:%S')] [阶段4] 跳过 (已存在): $task_label" >> "$task_log"
    fi

    # 全部完成
    touch "$STATUS_DIR/${task_label}.OK"
    return 0
}

# ─── 批量调度 (恒定 N 并发) ────────────────────────────
run_batches() {
    local task_count=0

    for mae in "$MAE_DIR"/*.mae; do
        [ -f "$mae" ] || continue
        mae_base=$(basename "$mae" .mae)

        for pdb in "$PPD_DIR"/*.pdb; do
            [ -f "$pdb" ] || continue
            pdb_base=$(basename "$pdb" .pdb)
            out_base="$OUTPUT_ROOT/$mae_base"
            mkdir -p "$out_base/loop" "$out_base/mini"

            task_count=$((task_count + 1))

            # 等待直到有空闲槽位
            while [ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]; do
                wait -n 2>/dev/null || true
            done

            # 启动后台任务
            (
                run_one "$pdb" "$mae" "$out_base"
            ) &

            # 每启动 50 个输出一次进度
            if [ $((task_count % 50)) -eq 0 ]; then
                local running
                running=$(jobs -rp | wc -l)
                log "已提交 ${task_count} 个, 运行中 ${running}"
            fi
        done
    done

    # 等待所有剩余任务完成
    log "全部 ${task_count} 个任务已提交, 等待剩余完成..."
    wait
    echo ""
}

# ─── 主流程 ─────────────────────────────────────────────
main() {
    echo ""
    echo "  ╔═══════════════════════════════════════════╗"
    echo "  ║   PROTAC Ternary 完整流水线               ║"
    echo "  ║   4 阶段: pipeline → modify → loop → mini ║"
    echo "  ╚═══════════════════════════════════════════╝"
    echo ""
    echo "  并发:     ${MAX_JOBS}"
    echo "  PDB目录:  ${PPD_DIR}"
    echo "  MAE目录:  ${MAE_DIR}"
    echo "  输出目录: ${OUTPUT_ROOT}"
    echo "  日志目录: ${LOG_DIR}"
    echo ""
    echo "  环境配置:"
    echo "    BABEL:       ${BABEL}"
    echo "    ROSETTA_BIN: ${ROSETTA_BIN}"
    echo ""

    # 统计任务数
    local total_pdb total_mae
    total_pdb=$(find "$PPD_DIR" -maxdepth 1 -name '*.pdb' | wc -l)
    total_mae=$(find "$MAE_DIR" -maxdepth 1 -name '*.mae' | wc -l)
    echo "  PDB: ${total_pdb}   MAE: ${total_mae}   → 共 $((total_pdb * total_mae)) 个任务"
    echo ""

    local t0=$SECONDS
    run_batches
    summarize
    extract_summary

    local elapsed=$(( SECONDS - t0 ))
    echo "  总耗时: ${elapsed}s"
    echo ""
}

# ─── 入口 ──────────────────────────────────────────────
main
