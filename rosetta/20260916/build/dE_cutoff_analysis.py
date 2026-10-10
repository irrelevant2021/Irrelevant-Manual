#!/usr/bin/env python3
"""
ΔE cutoff scanning analysis for PROTAC ternary structure modeling quality.

Analyzes the relationship between Final_Score and pIC50 across different ΔE cutoffs,
using three pose aggregation strategies.

使用方式 (与其他 build/ 脚本风格一致):
  cd /path/to/project/build
  mamba run -n openfe_env python3 dE_cutoff_analysis.py \\
      --input  ../output/enhanced_summary.csv \\
      --output ../output/cutoff_analysis

输出:
  {output}/
    ├── cutoff_summary.csv          — 汇总表 (R², Pearson r, N, p)
    ├── fig1_R2_vs_cutoff.png       — 图1: R² vs ΔE cutoff (三种策略叠加)
    ├── fig2_scatter_matrix.png     — 图2: Final_Score vs pIC50 (逐cutoff)
    └── fig3_dE_vs_FinalScore.png   — 图3: ΔE vs Final_Score (按pIC50着色)
"""

import argparse
import os
import sys
import warnings
import random

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
        description='ΔE cutoff scanning for PROTAC ternary structure modeling quality')
    parser.add_argument('--input', required=True,
                        help='输入 CSV 路径 (enhanced_summary.csv, 需含 pIC50 列)')
    parser.add_argument('--output', required=True,
                        help='输出目录 (自动创建)')
    parser.add_argument('--pIC50-col', default='pIC50_Vcap',
                        help='pIC50 列名 (默认: pIC50_Vcap)')
    parser.add_argument('--backbone-rms-max', type=float, default=2.0,
                        help='Backbone_RMS 过滤阈值 (默认: 2.0)')
    parser.add_argument('--ligand-rmsd-max', type=float, default=1.0,
                        help='Ligand_RMSD 过滤阈值 (默认: 1.0)')
    parser.add_argument('--final-score-max', type=float, default=None,
                        help='Final_Score 过滤上限 (默认: None, 即不进行过滤)')
    parser.add_argument('--seed', type=int, default=1,
                        help='随机种子 (默认: 1)')
    parser.add_argument('--cutoffs', type=str, default=None,
                        help='ΔE cutoff 列表, 逗号分隔 (默认: 0,20,40,...,200,9999)')
    return parser.parse_args()


def main(args):
    os.makedirs(args.output, exist_ok=True)

    # ── Random seed ──
    random.seed(args.seed)
    np.random.seed(args.seed)

    # ── ΔE cutoffs ──
    if args.cutoffs:
        CUTOFFS = sorted(set(int(c) for c in args.cutoffs.split(',') if c.strip()))
    else:
        CUTOFFS = sorted(set([0, 20, 40, 60, 80, 100] + list(range(110, 201, 10)) + [9999]))

    # ─── 1. Load & QC ────────────────────────────────────────────────────
    print("=" * 60)
    print("Step 1: Loading data and applying QC filter")
    print("=" * 60)

    pIC50_col = args.pIC50_col
    df = pd.read_csv(args.input)
    print(f"Total rows loaded: {len(df)}")

    # Use 'Code' column as compound identifier (MAE is empty for some rows)
    df.rename(columns={'Code': 'Compound'}, inplace=True)

    # Drop rows with empty pIC50
    before = len(df)
    df = df.dropna(subset=[pIC50_col])
    print(f"Dropped {before - len(df)} rows with missing pIC50")

    # QC filter
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

    # Ensure ΔE is numeric
    df['ΔE'] = pd.to_numeric(df['ΔE'], errors='coerce')
    df = df.dropna(subset=['ΔE'])

    # Sort for deterministic output
    df = df.sort_values(['Compound', 'PDB']).reset_index(drop=True)

    unique_compounds = df['Compound'].unique()
    print(f"Unique compounds after QC: {len(unique_compounds)}")
    print(f"pIC50 range: {df[pIC50_col].min():.3f} - {df[pIC50_col].max():.3f}")
    print(f"DE  range:   {df['ΔE'].min():.3f} - {df['ΔE'].max():.3f}")
# ─── 2. Aggregation helpers (closure over args) ─────────────────────────

    def aggregate_best_final_score(group):
        """Scenario A: pick the row with minimal Final_Score per compound."""
        return group.loc[group['Final_Score'].idxmin()]

    def aggregate_best_dE(group):
        """Scenario B: pick the row with minimal DE per compound."""
        return group.loc[group['ΔE'].idxmin()]

    def aggregate_random(group):
        """Scenario C: pick a random row per compound."""
        return group.sample(n=1, random_state=args.seed)

    # ─── 3. DE cutoff scanning ──────────────────────────────────────────
    print("=" * 60)
    print("Step 2 & 3: DE cutoff scanning with three aggregation strategies")
    print("=" * 60)

    strategies = {
        'Final_Score_best': aggregate_best_final_score,
        'ΔE_best':          aggregate_best_dE,
        'Random':           aggregate_random,
    }

    results = []  # list of dicts

    for strategy_name, agg_func in strategies.items():
        print(f"  Strategy: {strategy_name}")

        for cutoff in CUTOFFS:
            # Filter by DE < cutoff
            if cutoff == 9999:
                sub = df.copy()
            else:
                sub = df[df['ΔE'] < cutoff].copy()

            if sub.empty:
                print(f"    cutoff={cutoff:5d}: no data -> skip")
                continue

            # Aggregate: one row per compound
            aggregated = sub.groupby('Compound', sort=False).apply(agg_func)
            aggregated = aggregated.reset_index(drop=True)

            N = len(aggregated)
            if N < 3:
                print(f"    cutoff={cutoff:5d}: N={N} (<3) -> skip")
                continue

            # Linear regression: pIC50_Vcap ~ Final_Score
            X = aggregated['Final_Score'].values.reshape(-1, 1)
            y = aggregated[pIC50_col].values

            reg = LinearRegression()
            reg.fit(X, y)
            y_pred = reg.predict(X)

            r2 = r2_score(y, y_pred)
            pearson_r, p_value = stats.pearsonr(y.flatten(), y_pred.flatten())

            results.append({
                'cutoff':       cutoff,
                'N_compounds':  N,
                'R²':           round(r2, 4),
                'Pearson_r':    round(pearson_r, 4),
                'p_value':      f"{p_value:.4e}",
                'aggregation':  strategy_name,
            })

            print(f"    cutoff={cutoff:5d}: N={N:2d},  R2={r2:.4f},  "
                  f"Pearson_r={pearson_r:.4f},  p={p_value:.4e}")

    results_df = pd.DataFrame(results)
    print(f"\n{'=' * 60}")
    print("Summary table:")
    print(results_df[['cutoff', 'N_compounds', 'R²', 'Pearson_r', 'p_value', 'aggregation']].to_string(index=False))
    # Save summary table
    table_path = os.path.join(args.output, 'cutoff_summary.csv')
    results_df.to_csv(table_path, index=False)
    print(f"\nSummary table saved to: {table_path}")

    # ─── 4. Visualization ────────────────────────────────────────────────────
    print(f"\n{'=' * 60}")
    print("Step 4: Generating figures")
    print("=" * 60)

    # Color palette for strategies
    STRATEGY_COLORS = {
        'Final_Score_best': '#1f77b4',
        'ΔE_best':          '#d62728',
        'Random':           '#2ca02c',
    }
    STRATEGY_LABELS = {
        'Final_Score_best': 'Min Final_Score',
        'ΔE_best':          'Min ΔE',
        'Random':           'Random',
    }
    print()
    # ── Figure 1: R² vs ΔE cutoff (three strategies overlaid) ────────────────
    fig1, ax1 = plt.subplots(figsize=(10, 6))

    # Baseline: R2 on all data (cutoff=9999) for each strategy
    baselines = {}
    for s in strategies:
        row = results_df[(results_df['aggregation'] == s) &
                         (results_df['cutoff'] == 9999)]
        if not row.empty:
            baselines[s] = row['R²'].values[0]

    for strategy_name in strategies:
        sdata = results_df[results_df['aggregation'] == strategy_name].copy()
        if sdata.empty:
            continue
        # Exclude cutoff=9999 for the main line (show as separate baseline)
        plot_data = sdata[sdata['cutoff'] != 9999].sort_values('cutoff')
        full_data = sdata[sdata['cutoff'] == 9999]

        ax1.plot(plot_data['cutoff'], plot_data['R²'],
                 marker='o', linestyle='-', linewidth=2, markersize=8,
                 color=STRATEGY_COLORS[strategy_name],
                 label=STRATEGY_LABELS[strategy_name])

        # Annotate N
        for _, row in plot_data.iterrows():
            ax1.annotate(f"N={int(row['N_compounds'])}",
                         (row['cutoff'], row['R²']),
                         textcoords="offset points", xytext=(0, 12),
                         fontsize=7, ha='center',
                         color=STRATEGY_COLORS[strategy_name])

        # Baseline horizontal line
        if strategy_name in baselines:
            ax1.axhline(y=baselines[strategy_name], linestyle='--',
                        linewidth=1, alpha=0.5,
                        color=STRATEGY_COLORS[strategy_name])
            # Label baseline
            ax1.text(plot_data['cutoff'].max() + 5, baselines[strategy_name],
                     f"All data\nR²={baselines[strategy_name]:.3f}",
                     fontsize=8, color=STRATEGY_COLORS[strategy_name],
                     va='center')

    ax1.set_xlabel('ΔE Cutoff (kJ/mol)', fontsize=12)
    ax1.set_ylabel('R² (pIC50 ~ Final_Score)', fontsize=12)
    ax1.set_title('Figure 1: R² vs ΔE Cutoff by Aggregation Strategy', fontsize=14)
    ax1.legend(fontsize=10, loc='best')
    ax1.grid(True, alpha=0.3)
    fig1.tight_layout()
    fig1.savefig(os.path.join(args.output, 'fig1_R2_vs_cutoff.png'), dpi=150)
    print("  Figure 1 saved: fig1_R2_vs_cutoff.png")
    # ── Figure 2: Scatter matrix per cutoff interval (Scenario A only) ──────
    N_CUTOFFS = len(CUTOFFS)
    N_COLS = 5
    N_ROWS = int(np.ceil(N_CUTOFFS / N_COLS))
    fig2, axes2 = plt.subplots(N_ROWS, N_COLS, figsize=(5 * N_COLS, 4 * N_ROWS))
    axes2 = axes2.flatten()

    # Hide unused subplots
    for ax in axes2[N_CUTOFFS:]:
        ax.set_visible(False)

    # Only Scenario A
    for idx, cutoff in enumerate(CUTOFFS):
        ax = axes2[idx]

        if cutoff == 9999:
            sub = df.copy()
            cutoff_label = "All data"
        else:
            sub = df[df['ΔE'] < cutoff].copy()
            cutoff_label = f"ΔE < {cutoff}"

        if sub.empty:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', fontsize=12)
            ax.set_title(f"{cutoff_label}  (N=0)", fontsize=10)
            continue

        # Aggregate Scenario A
        agg = sub.groupby('Compound', sort=False).apply(aggregate_best_final_score)
        agg = agg.reset_index(drop=True)

        N = len(agg)
        X = agg['Final_Score'].values
        y = agg[pIC50_col].values

        ax.scatter(X, y, alpha=0.7, s=40, c='#1f77b4', edgecolors='k', linewidths=0.5)

        if N >= 3:
            # Regression line
            slope, intercept, r_val, p_val, std_err = stats.linregress(X, y)
            x_line = np.linspace(X.min(), X.max(), 100)
            y_line = slope * x_line + intercept
            ax.plot(x_line, y_line, 'r-', linewidth=2)
            ax.text(0.05, 0.95, f"R²={r_val**2:.3f}\nN={N}",
                    transform=ax.transAxes, fontsize=10,
                    verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        ax.set_xlabel('Final_Score', fontsize=9)
        ax.set_ylabel('pIC50', fontsize=9)
        ax.set_title(f"{cutoff_label}  (N={N})", fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3)

    fig2.suptitle('Figure 2: Final_Score vs pIC50 at Different ΔE Cutoffs (Scenario A: Min Final_Score)',
                  fontsize=14, y=1.02)
    fig2.tight_layout()
    fig2.savefig(os.path.join(args.output, 'fig2_scatter_matrix.png'), dpi=150, bbox_inches='tight')
    print("  Figure 2 saved: fig2_scatter_matrix.png")
    # ── Figure 3: Joint distribution of ΔE vs Final_Score ──────────────────
    fig3, ax3 = plt.subplots(figsize=(10, 7))

    # Use all QC-passed data
    sc = ax3.scatter(df['ΔE'], df['Final_Score'],
                     c=df[pIC50_col], cmap='viridis_r',
                     alpha=0.4, s=15, edgecolors='none')

    cbar = plt.colorbar(sc, ax=ax3, label='pIC50')
    ax3.set_xlabel('ΔE (kJ/mol)', fontsize=12)
    ax3.set_ylabel('Final_Score', fontsize=12)
    ax3.set_title('Figure 3: ΔE vs Final_Score (colored by pIC50)', fontsize=14)
    ax3.grid(True, alpha=0.3)

    # Pearson correlation between DE and Final_Score
    r_val, p_val = stats.pearsonr(df['ΔE'], df['Final_Score'])
    ax3.text(0.95, 0.95, f"Pearson r = {r_val:.4f}\np = {p_val:.4e}",
             transform=ax3.transAxes, fontsize=11,
             verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

    fig3.tight_layout()
    fig3.savefig(os.path.join(args.output, 'fig3_dE_vs_FinalScore.png'), dpi=150)
    print("  Figure 3 saved: fig3_dE_vs_FinalScore.png")

    print(f"\n{'=' * 60}")
    print(f"All outputs saved to: {args.output}")
if __name__ == '__main__':
    main(parse_args())
print(f"{'=' * 60}")