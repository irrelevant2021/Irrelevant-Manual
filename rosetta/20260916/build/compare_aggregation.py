#!/usr/bin/env python3
"""
compare_aggregation.py — 比较多种多构象聚合方法对 ΔE cutoff 分析的影响

用途:
  在同一套数据上，比较四种聚合策略（winner / top_k / energy_window / boltzmann）
  在不同 ΔE cutoff 下，pIC50 与 Final_Score 回归的 R² 表现。

使用方式 (与其他 build/ 脚本风格一致):
  cd /path/to/project/build
  mamba run -n openfe_env python3 compare_aggregation.py \\
      --input  ../output/enhanced_summary.csv \\
      --output ../output/aggregation_comparison

输出:
  {output}/
    ├── comparison_summary.csv     — 所有方法 × 所有 cutoff 汇总
    ├── fig1_R2_vs_cutoff.png     — 四种方法 R² 曲线叠加
    ├── fig2_best_scatter.png     — 最佳方法 (全局 R² 最高) 散点矩阵
    │
    ├── winner/                    — 各方法子目录（含详细输出）
    ├── top_k/
    ├── energy_window/
    └── boltzmann/
"""

import argparse
import os
import sys
import warnings
import random
import copy

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

warnings.filterwarnings('ignore')


def parse_args():
    parser = argparse.ArgumentParser(
        description='Compare multi-pose aggregation methods for ΔE cutoff analysis')
    parser.add_argument('--input', required=True,
                        help='输入 CSV 路径 (enhanced_summary.csv)')
    parser.add_argument('--output', required=True,
                        help='输出根目录 (自动创建，各方法存入子目录)')
    parser.add_argument('--pIC50-col', default='pIC50_Vcap',
                        help='pIC50 列名 (默认: pIC50_Vcap)')
    parser.add_argument('--methods', type=str, default='winner,top_k,energy_window,boltzmann,top_pct',
                        help='聚合方法列表 (默认: winner,top_k,energy_window,boltzmann,top_pct)')
    parser.add_argument('--top-k', type=int, default=3,
                        help='Top-K 法的 K 值 (默认: 3)')
    parser.add_argument('--top-pct', type=float, default=2.0,
                        help='Top-N%% 法的百分比 (默认: 2.0, 即前 2%%)')
    parser.add_argument('--energy-window', type=float, default=10.0,
                        help='能量窗口法窗口大小 (REU, 默认: 10.0)')
    parser.add_argument('--boltzmann-kt', type=float, default=600.0,
                        help='Boltzmann 加权 kT (REU 单位, 默认: 600.0)')
    parser.add_argument('--backbone-rms-max', type=float, default=2.0,
                        help='Backbone_RMS 过滤阈值 (默认: 2.0)')
    parser.add_argument('--ligand-rmsd-max', type=float, default=1.0,
                        help='Ligand_RMSD 过滤阈值 (默认: 1.0)')
    parser.add_argument('--final-score-max', type=float, default=None,
                        help='Final_Score 过滤上限 (默认: None, 即不进行过滤)')
    parser.add_argument('--seed', type=int, default=42,
                        help='随机种子 (默认: 42)')
    parser.add_argument('--cutoffs', type=str, default=None,
                        help='ΔE cutoff 列表, 逗号分隔 (默认: 0,20,40,...,200,9999)')
    return parser.parse_args()
def aggregate_method(group, method_name, args):
    """
    对单个化合物 (group, 含多行 PDB) 执行指定聚合方法。
    返回聚合后的一行数据 (Series), 包含 Final_Score, ΔE 等。
    """
    if method_name == 'winner':
        return group.loc[group['Final_Score'].idxmin()]

    elif method_name == 'top_k':
        k = min(args.top_k, len(group))
        top = group.nsmallest(k, 'Final_Score')
        row = group.iloc[0].copy()
        for col in ['Final_Score', 'ΔE', 'Best_Energy_kJ', 'E_free']:
            if col in group.columns:
                row[col] = top[col].mean()
        return row

    elif method_name == 'energy_window':
        min_final = group['Final_Score'].min()
        window = group[group['Final_Score'] <= min_final + args.energy_window]
        if window.empty:
            window = group.nsmallest(1, 'Final_Score')
        row = group.iloc[0].copy()
        for col in ['Final_Score', 'ΔE', 'Best_Energy_kJ', 'E_free']:
            if col in group.columns:
                row[col] = window[col].mean()
        return row

    elif method_name == 'boltzmann':
        kt = args.boltzmann_kt
        scores = group['Final_Score'].values
        scores_shifted = scores - scores.min()
        weights = np.exp(-scores_shifted / kt)
        if weights.sum() == 0:
            weights = np.ones_like(weights) / len(weights)
        else:
            weights = weights / weights.sum()
        row = group.iloc[0].copy()
        for col in ['Final_Score', 'ΔE', 'Best_Energy_kJ', 'E_free']:
            if col in group.columns:
                row[col] = np.average(group[col].values, weights=weights)
        return row

    elif method_name == 'top_pct':
        pct = args.top_pct / 100.0
        k = max(1, int(len(group) * pct))
        top = group.nsmallest(k, 'Final_Score')
        row = group.iloc[0].copy()
        for col in ['Final_Score', 'ΔE', 'Best_Energy_kJ', 'E_free']:
            if col in group.columns:
                row[col] = top[col].mean()
        return row

    else:
        raise ValueError(f"Unknown method: {method_name}")

def main():
    args = parse_args()
    os.makedirs(args.output, exist_ok=True)

    random.seed(args.seed)
    np.random.seed(args.seed)

    if args.cutoffs:
        CUTOFFS = sorted(set(int(c) for c in args.cutoffs.split(',') if c.strip()))
    else:
        CUTOFFS = sorted(set([0, 20, 40, 60, 80, 100] + list(range(120, 251, 10)) + [9999]))

    selected_methods = [m.strip() for m in args.methods.split(',') if m.strip()]
    print(f"  Methods: {selected_methods}")

    # ─── 1. Load & QC ───────────────────────────────────────────────────
    print("=" * 60)
    print("Step 1: Loading data and applying QC filter")
    print("=" * 60)

    pIC50_col = args.pIC50_col
    df = pd.read_csv(args.input)
    print(f"Total rows loaded: {len(df)}")

    if 'Code' in df.columns and 'Compound' not in df.columns:
        df.rename(columns={'Code': 'Compound'}, inplace=True)
    elif 'MAE' in df.columns and 'Compound' not in df.columns:
        df.rename(columns={'MAE': 'Compound'}, inplace=True)

    before = len(df)
    df = df.dropna(subset=[pIC50_col])
    print(f"Dropped {before - len(df)} rows with missing pIC50")

    before_qc = len(df)
    df = df[(df['Backbone_RMS'] < args.backbone_rms_max) &
            (df['Ligand_RMSD']  < args.ligand_rmsd_max)].copy()
    print(f"QC filter (Backbone_RMS<{args.backbone_rms_max}, "
          f"Ligand_RMSD<{args.ligand_rmsd_max}): "
          f"{before_qc} -> {len(df)} rows")

    # Final_Score filter (< args.final_score_max, 仅当设置了有效值)
    if args.final_score_max is not None:
        before_fs = len(df)
        df = df[df['Final_Score'] < args.final_score_max].copy()
        print(f"Final_Score filter (< {args.final_score_max}): "
              f"{before_fs} -> {len(df)} rows")
    else:
        print("Final_Score filter: disabled (--final-score-max is None)")

    df['ΔE'] = pd.to_numeric(df['ΔE'], errors='coerce')
    df = df.dropna(subset=['ΔE'])
    df = df.sort_values(['Compound', 'PDB']).reset_index(drop=True)

    unique_compounds = df['Compound'].unique()
    print(f"Unique compounds after QC: {len(unique_compounds)}")
    print(f"pIC50 range: {df[pIC50_col].min():.3f} - {df[pIC50_col].max():.3f}")
    print(f"ΔE  range:   {df['ΔE'].min():.3f} - {df['ΔE'].max():.3f}")
# ─── 2. Multi-method ΔE cutoff scanning ───────────────────────────
    all_results = []

    for method_name in selected_methods:
        print(f"\n  Method: {method_name}")
        print(f"  {'─' * 50}")

        method_dir = os.path.join(args.output, method_name)
        os.makedirs(method_dir, exist_ok=True)

        method_results = []

        for cutoff in CUTOFFS:
            if cutoff == 9999:
                sub = df.copy()
            else:
                sub = df[df['ΔE'] < cutoff].copy()

            if sub.empty:
                print(f"    cutoff={cutoff:5d}: no data -> skip")
                continue

            groups = sub.groupby('Compound', sort=False)
            agg_rows = []
            for _, g in groups:
                try:
                    agg_rows.append(aggregate_method(g, method_name, args))
                except Exception:
                    continue
            if not agg_rows:
                continue
            agg = pd.DataFrame(agg_rows).reset_index(drop=True)

            N = len(agg)
            if N < 3:
                print(f"    cutoff={cutoff:5d}: N={N} (<3) -> skip")
                continue

            X = agg['Final_Score'].values.reshape(-1, 1)
            y = agg[pIC50_col].values
            reg = LinearRegression()
            reg.fit(X, y)
            r2 = r2_score(y, reg.predict(X))
            _, _, r_val, p_val, _ = stats.linregress(X.flatten(), y)

            method_results.append({
                'cutoff':      cutoff,
                'N_compounds': N,
                'R²':          round(r2, 4),
                'Pearson_r':   round(r_val, 4),
                'p_value':     f"{p_val:.4e}",
                'method':      method_name,
            })
            print(f"    cutoff={cutoff:5d}: N={N:2d}, R²={r2:.4f}, "
                  f"r={r_val:.4f}, p={p_val:.4e}")

        if method_results:
            pd.DataFrame(method_results).to_csv(
                os.path.join(method_dir, 'cutoff_summary.csv'), index=False)
            all_results.extend(method_results)
    if not all_results:
        print("\n  ❌ No results generated. Check your data.")
        sys.exit(1)

    all_df = pd.DataFrame(all_results)
    all_df.to_csv(os.path.join(args.output, 'comparison_summary.csv'), index=False)
    print(f"\n  ✅ comparison_summary.csv saved")

    print(f"\n  Best R² by method:")
    for method in selected_methods:
        mrows = all_df[all_df['method'] == method]
        if not mrows.empty:
            best = mrows.loc[mrows['R²'].idxmax()]
            print(f"    {method:20s}: R²={best['R²']:.4f} at cutoff={best['cutoff']:>5d}, "
                  f"N={best['N_compounds']}, p={best['p_value']}")
# ─── 3. Figures ─────────────────────────────────────────────────────
    print(f"\n{'=' * 60}")
    print("Generating figures")
    print("=" * 60)

    METHOD_COLORS = {
        'winner':         '#1f77b4',
        'top_k':          '#d62728',
        'top_pct':        '#9467bd',
        'energy_window':  '#2ca02c',
        'boltzmann':      '#ff7f0e',
    }
    METHOD_LABELS = {
        'winner':         'Winner (min Final)',
        'top_k':          f'Top-{args.top_k} Avg',
        'top_pct':        f'Top {args.top_pct}% Avg',
        'energy_window':  f'Window (±{args.energy_window})',
        'boltzmann':      f'Boltzmann (kT={args.boltzmann_kt})',
    }

    # ── Figure 1: R² vs ΔE cutoff (overlay) ──────────────────────────
    fig1, ax1 = plt.subplots(figsize=(12, 7))
    for method_name in selected_methods:
        mdata = all_df[all_df['method'] == method_name]
        if mdata.empty:
            continue
        plot_data = mdata[mdata['cutoff'] != 9999].sort_values('cutoff')
        if plot_data.empty:
            continue
        color = METHOD_COLORS.get(method_name, '#333')
        label = METHOD_LABELS.get(method_name, method_name)

        ax1.plot(plot_data['cutoff'], plot_data['R²'],
                 marker='o', linestyle='-', linewidth=2, markersize=7,
                 color=color, label=label)
        for _, row in plot_data.iterrows():
            ax1.annotate(f"{int(row['N_compounds'])}",
                         (row['cutoff'], row['R²']),
                         textcoords="offset points", xytext=(0, 12),
                         ha='center', fontsize=7, color=color)

        baseline = mdata[mdata['cutoff'] == 9999]
        if not baseline.empty:
            ax1.axhline(y=baseline['R²'].values[0], linestyle='--',
                        color=color, alpha=0.4, linewidth=1)
            ax1.text(plot_data['cutoff'].max() + 5,
                     baseline['R²'].values[0],
                     f"{label} baseln", fontsize=8, color=color, alpha=0.6,
                     va='center')

    ax1.set_xlabel('ΔE Cutoff (kJ/mol)', fontsize=12)
    ax1.set_ylabel('R² (pIC50 ~ Final_Score)', fontsize=12)
    ax1.set_title('Figure 1: R² vs ΔE Cutoff — Aggregation Method Comparison',
                  fontsize=14)
    ax1.legend(fontsize=10, loc='best')
    ax1.grid(True, alpha=0.3)
    fig1.tight_layout()
    fig1.savefig(os.path.join(args.output, 'fig1_R2_vs_cutoff.png'), dpi=150)
    print("  Figure 1 saved: fig1_R2_vs_cutoff.png")
# ── Figure 2: Best method scatter matrix ─────────────────────────
    best_method = None
    best_r2 = -1
    for method_name in selected_methods:
        mrows = all_df[all_df['method'] == method_name]
        if not mrows.empty:
            mr2 = mrows['R²'].max()
            if mr2 > best_r2:
                best_r2 = mr2
                best_method = method_name

    print(f"\n  Best method overall: {best_method} (R²={best_r2:.4f})")

    if best_method:
        N_CUTOFFS = len(CUTOFFS)
        N_COLS = 5
        N_ROWS = int(np.ceil(N_CUTOFFS / N_COLS))
        fig2, axes2 = plt.subplots(N_ROWS, N_COLS,
                                   figsize=(5 * N_COLS, 4 * N_ROWS))
        axes2 = axes2.flatten()
        for ax in axes2[N_CUTOFFS:]:
            ax.set_visible(False)

        for idx, cutoff in enumerate(CUTOFFS):
            ax = axes2[idx]
            if cutoff == 9999:
                sub, label = df.copy(), "All data"
            else:
                sub, label = df[df['ΔE'] < cutoff].copy(), f"ΔE < {cutoff}"

            if sub.empty:
                ax.text(0.5, 0.5, 'No data', ha='center', va='center',
                        fontsize=12)
                ax.set_title(f"{label}  (N=0)", fontsize=10)
                continue

            groups = sub.groupby('Compound', sort=False)
            agg_rows = []
            for _, g in groups:
                try:
                    agg_rows.append(aggregate_method(g, best_method, args))
                except Exception:
                    continue
            if not agg_rows:
                continue
            agg = pd.DataFrame(agg_rows).reset_index(drop=True)
            N = len(agg)
            X = agg['Final_Score'].values
            y = agg[pIC50_col].values

            ax.scatter(X, y, alpha=0.7, s=40, c='#1f77b4',
                       edgecolors='k', linewidths=0.5)
            if N >= 3:
                slope, icpt, r_val, p_val, _ = stats.linregress(X, y)
                xl = np.linspace(X.min(), X.max(), 100)
                ax.plot(xl, slope * xl + icpt, 'r-', linewidth=2)
                ax.text(0.05, 0.95, f"R²={r_val**2:.3f}\nN={N}",
                        transform=ax.transAxes, fontsize=10,
                        verticalalignment='top',
                        bbox=dict(boxstyle='round', facecolor='wheat',
                                  alpha=0.5))
            ax.set_xlabel('Final_Score', fontsize=9)
            ax.set_ylabel('pIC50', fontsize=9)
            ax.set_title(f"{label}  (N={N})", fontsize=11, fontweight='bold')
            ax.grid(True, alpha=0.3)

        fig2.suptitle(f'Figure 2: Final_Score vs pIC50 '
                      f'(Method: {best_method}, R²={best_r2:.3f})',
                      fontsize=14, y=1.02)
        fig2.tight_layout()
        fig2.savefig(os.path.join(args.output, 'fig2_best_scatter.png'),
                     dpi=150, bbox_inches='tight')
        print("  Figure 2 saved: fig2_best_scatter.png")

    print(f"\n{'=' * 60}")
    print(f"All outputs saved to: {args.output}")
    print(f"{'=' * 60}")


if __name__ == '__main__':
    main()