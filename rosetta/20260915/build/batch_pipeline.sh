#!/usr/bin/env bash
# ============================================================
# batch_pipeline.sh
# PROTAC Ternary 批量调度脚本 — 全组合 (所有 PDB × 所有 MAE)
#
# 用法:
#   ./batch_pipeline.sh                              # 全部默认
#   ./batch_pipeline.sh 32                           # 32 并发
#   ./batch_pipeline.sh 64 ../ppd ../linkers         # 指定目录
#   ./batch_pipeline.sh 64 ../ppd ../linkers ../outputs
#
#   # 自定义对齐核心
#   ./batch_pipeline.sh --cbn-core MCS --poi-core MCS
#   ./batch_pipeline.sh --cbn-core MCS --poi-core MCS 32 ../ppd ../linkers
#
# 输出结构:
#   ../outputs/{MAE_NAME}/{PDB_NAME}_refined.pdb
#
# 特性:
#   - 恒定 N 并发 (后台进程池)
#   - 输出已存在时自动跳过 (可恢复)
#   - 单任务失败不影响其他
#   - 实时进度 + 最终汇总表
# ============================================================
set -o pipefail

# ─── 参数 (含默认值) ────────────────────────────────────
MAX_JOBS=64
PPD_DIR="../ppd"
MAE_DIR="../linkers"
OUTPUT_ROOT="../outputs"
CBN_CORE="MCS"
POI_CORE="N#Cc1ccccc1"

# 解析命名参数 (必须出现在位置参数之前)
while [[ $# -gt 0 ]]; do
    case "$1" in
        --cbn-core) CBN_CORE="$2"; shift 2 ;;
        --poi-core) POI_CORE="$2"; shift 2 ;;
        --help|-h)
            echo "用法: $0 [--cbn-core MCS] [--poi-core SMARTS] [并发数] [pdb目录] [mae目录] [输出目录]"
            exit 0 ;;
        --) shift; break ;;
        *) break ;;
    esac
done

# 剩余位置参数
[ $# -ge 1 ] && MAX_JOBS="$1"
[ $# -ge 2 ] && PPD_DIR="$2"
[ $# -ge 3 ] && MAE_DIR="$3"
[ $# -ge 4 ] && OUTPUT_ROOT="$4"

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
ENV_NAME="molscribe"
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
    local elapsed=$(( $(date +%s) - START_TIME ))
    local total=$(( ok + fail + skip ))
    echo ""
    echo "  ═══════════════════════════════════════════"
    echo "   汇总:  成功 ${ok}   |   跳过 ${skip}   |   失败 ${fail}"
    echo "   总计:  ${total} 个任务  |  耗时: ${elapsed}s"
    echo "  ═══════════════════════════════════════════"
    echo ""
    if [ "$fail" -gt 0 ]; then
        echo "  ❌ 失败任务列表:"
        find "$STATUS_DIR" -name '*.FAIL' -printf '    - %f\n' | sed 's/\.FAIL$//'
        echo "    日志路径: $LOG_DIR/{mae}/{pdb}.log"
        echo ""
    fi
}

run_one() {
    local pdb="$1"
    local mae="$2"
    local output="$3"
    local pdb_base; pdb_base="$(basename "$pdb" .pdb)"
    local mae_base; mae_base="$(basename "$mae" .mae)"
    local task_label="${pdb_base}__${mae_base}"
    local task_log="$LOG_DIR/${mae_base}/${pdb_base}.log"

    mkdir -p "$(dirname "$task_log")"

    # 跳过已完成的
    if [ -s "$output" ]; then
        touch "$STATUS_DIR/${task_label}.SKIP"
        return 0
    fi

    # 运行 pipeline
    if mamba run -n "$ENV_NAME" python3 "$PIPELINE" \
        --pdb  "$pdb" \
        --mae  "$mae" \
        -o     "$output" \
        --cbn-core "$CBN_CORE" \
        --poi-core "$POI_CORE" \
        > "$task_log" 2>&1
    then
        touch "$STATUS_DIR/${task_label}.OK"
        return 0
    else
        touch "$STATUS_DIR/${task_label}.FAIL"
        echo "    [FAIL] ${task_label}  —  tail -5 $task_log"
        tail -5 "$task_log" | sed 's/^/           /'
        return 1
    fi
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
            out_dir="$OUTPUT_ROOT/$mae_base"
            mkdir -p "$out_dir"
            output="$out_dir/${pdb_base}_refined.pdb"

            task_count=$((task_count + 1))

            # 等待直到有空闲槽位
            while [ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]; do
                wait -n 2>/dev/null || true
            done

            # 启动后台任务
            (
                run_one "$pdb" "$mae" "$output"
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
    echo "  ║   PROTAC Ternary 批量流水线               ║"
    echo "  ║   全组合 (PDB × MAE)                      ║"
    echo "  ╚═══════════════════════════════════════════╝"
    echo ""

    # ── 任务规模 ──
    local total_pdb total_mae
    total_pdb=$(find "$PPD_DIR" -maxdepth 1 -name '*.pdb' | wc -l)
    total_mae=$(find "$MAE_DIR" -maxdepth 1 -name '*.mae' | wc -l)
    echo "  规模:   PDB ${total_pdb}  ×  MAE ${total_mae}  =  $((total_pdb * total_mae)) 个任务"
    echo "  并发:   ${MAX_JOBS}"
    echo "  目录:   PDB=${PPD_DIR}"
    echo "          MAE=${MAE_DIR}"
    echo "          输出=${OUTPUT_ROOT}"
    echo "          日志=${LOG_DIR}"
    echo ""

    # ── pipeline.py 参数表 ──
    echo "  ┌─────────────────────────────────────────────────────┐"
    echo "  │ pipeline.py 参数                                     │"
    echo "  ├──────────────────────────────┬──────────────────────┤"
    printf "  │  %-28s │  %-20s │\n" "--cbn-core (可设)"      "$CBN_CORE"
    printf "  │  %-28s │  %-20s │\n" "--poi-core (可设)"      "$POI_CORE"
    printf "  │  %-28s │  %-20s │\n" "--n-confs (固定)"       "200"
    printf "  │  %-28s │  %-20s │\n" "--seed (固定)"          "42"
    printf "  │  %-28s │  %-20s │\n" "--resname (固定)"       "PRT"
    printf "  │  %-28s │  %-20s │\n" "--chain (固定)"         "X"
    printf "  │  %-28s │  %-20s │\n" "--resseq (固定)"        "1"
    printf "  │  %-28s │  %-20s │\n" "--schrodinger (自动)"   '$SCHRODINGER'
    echo "  └──────────────────────────────┴──────────────────────┘"
    echo ""

    local t0=$SECONDS
    run_batches
    summarize

    local elapsed=$(( SECONDS - t0 ))
    echo "  总耗时: ${elapsed}s"
    echo ""
}

# ─── 入口 ──────────────────────────────────────────────
main
