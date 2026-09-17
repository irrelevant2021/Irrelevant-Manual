#!/usr/bin/env python
"""
计算PROTAC三元复合物对接结果中 CBN (链X) 与 POI (链Y) 配体之间的距离，
并从 score.sc 中提取 I_sc (界面打分) 合并到输出。

输出: cbn_poi_distances.csv
"""

import os
import glob
import shutil
import warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# ============ 配置 ============
PPD_DIR = "/data/AR/PROTAC_ternary-master/20260825/ppd"
TRASH_DIR = "/data/AR/PROTAC_ternary-master/20260825/ppd/trash"
SCORE_FILE = "/data/AR/PROTAC_ternary-master/20260825/score.sc"
OUTPUT_CSV = "/data/AR/PROTAC_ternary-master/20260825/cbn_poi_distances.csv"

# 筛选阈值
MIN_DIST_THRESHOLD = 25.0  # min_dist > 此值则移动
ISC_THRESHOLD = 0.0        # I_sc >= 此值则移动


def parse_score_sc(score_file):
    """
    解析 Rosetta score.sc 文件，提取每个模型的 I_sc。
    返回 dict: {description: I_sc}
    """
    iscores = {}

    with open(score_file, "r") as f:
        # 找到列名行
        header_line = None
        for line in f:
            if line.startswith("SCORE:") and "description" in line:
                header_line = line
                break

        if header_line is None:
            raise ValueError(f"无法在 {score_file} 中找到 SCORE 表头")

        # 解析列名
        cols = header_line.split()
        # 去掉第一个 "SCORE:"
        cols = cols[1:]

        # 找到 I_sc 和 description 的列索引
        try:
            isc_idx = cols.index("I_sc")
        except ValueError:
            raise ValueError("score.sc 中未找到 I_sc 列")

        try:
            desc_idx = cols.index("description")
        except ValueError:
            raise ValueError("score.sc 中未找到 description 列")

        # 读取数据行
        for line in f:
            if not line.startswith("SCORE:"):
                continue
            if "description" in line:
                continue

            parts = line.split()
            parts = parts[1:]  # 去掉 "SCORE:"

            if len(parts) <= max(isc_idx, desc_idx):
                continue

            description = parts[desc_idx]
            try:
                isc = float(parts[isc_idx])
            except (ValueError, IndexError):
                isc = np.nan

            iscores[description] = isc

    return iscores


def parse_pdb_fast(pdb_path):
    """
    快速解析PDB文件，直接提取CBN和POI的重原子坐标。
    """
    cbn_heavy = []
    poi_heavy = []

    with open(pdb_path, "r") as f:
        for line in f:
            if not line.startswith("HETATM"):
                continue
            resname = line[17:20].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            element = line[76:78].strip()

            if element == "H":
                continue

            if resname == "CBN":
                cbn_heavy.append([x, y, z])
            elif resname == "POI":
                poi_heavy.append([x, y, z])

    if len(cbn_heavy) == 0 or len(poi_heavy) == 0:
        return None

    return np.array(cbn_heavy), np.array(poi_heavy)


def calc_min_distance(coords1, coords2):
    """计算两组坐标之间的最小距离"""
    diff = coords1[:, np.newaxis, :] - coords2[np.newaxis, :, :]
    dists = np.sqrt(np.sum(diff ** 2, axis=-1))
    return dists.min()


def process_file(pdb_path):
    """处理单个PDB文件，返回最近重原子距离"""
    result = parse_pdb_fast(pdb_path)

    if result is None:
        return np.nan

    cbn_heavy, poi_heavy = result
    return round(calc_min_distance(cbn_heavy, poi_heavy), 4)


def main():
    # 1. 解析 score.sc 获取 I_sc
    print("正在解析 score.sc ...")
    iscores = parse_score_sc(SCORE_FILE)
    print(f"  从 score.sc 中读取了 {len(iscores)} 个模型的 I_sc")

    # 2. 计算距离
    pdb_files = sorted(glob.glob(os.path.join(PPD_DIR, "*.pdb")))
    print(f"找到 {len(pdb_files)} 个PDB文件，开始计算距离 ...")

    results = []
    total = len(pdb_files)

    for i, pdb_path in enumerate(pdb_files):
        filename = os.path.basename(pdb_path)
        model_name = filename.replace(".pdb", "")

        min_dist = process_file(pdb_path)
        isc = iscores.get(model_name, np.nan)

        results.append({
            "filename": filename,
            "min_dist": min_dist,
            "I_sc": isc,
        })

        if (i + 1) % 500 == 0 or (i + 1) == total:
            print(f"  进度: {i+1}/{total} ({(i+1)/total*100:.1f}%)")

    df = pd.DataFrame(results)

    # 3. 统计摘要
    print(f"\n{'='*60}")
    print(f"处理完成! 共 {len(df)} 个模型")
    print(f"{'='*60}")

    valid_dist = df["min_dist"].dropna()
    print(f"\n--- min_dist (CBN-POI 最近重原子距离, Å) ---")
    print(f"  mean={valid_dist.mean():.2f}, median={valid_dist.median():.2f}, "
          f"min={valid_dist.min():.2f}, max={valid_dist.max():.2f}")

    valid_isc = df["I_sc"].dropna()
    print(f"\n--- I_sc (界面打分, REU) ---")
    print(f"  mean={valid_isc.mean():.2f}, median={valid_isc.median():.2f}, "
          f"min={valid_isc.min():.2f}, max={valid_isc.max():.2f}")

    n_missing_dist = df["min_dist"].isna().sum()
    n_missing_isc = df["I_sc"].isna().sum()
    if n_missing_dist > 0:
        print(f"\n  ⚠ {n_missing_dist} 个文件缺少 CBN 或 POI 配体")
    if n_missing_isc > 0:
        print(f"\n  ⚠ {n_missing_isc} 个模型在 score.sc 中未找到 I_sc")

    # 4. 保存CSV
    df.to_csv(OUTPUT_CSV, index=False, float_format="%.4f")
    print(f"\n结果已保存至: {OUTPUT_CSV}")

    # 5. 显示前10行
    print(f"\n--- 前10个结果 ---")
    print(df.head(10).to_string(index=False))

    # 6. 移动不合格文件到 trash
    print(f"\n{'='*60}")
    print(f"开始筛选不合格模型 (min_dist > {MIN_DIST_THRESHOLD} 或 I_sc >= {ISC_THRESHOLD}) ...")
    
    # 创建 trash 目录
    os.makedirs(TRASH_DIR, exist_ok=True)
    
    # 找出需要移动的文件
    mask_move = (df["min_dist"] > MIN_DIST_THRESHOLD) | (df["I_sc"] >= ISC_THRESHOLD)
    to_move = df[mask_move]
    
    moved_count = 0
    for _, row in to_move.iterrows():
        src = os.path.join(PPD_DIR, row["filename"])
        dst = os.path.join(TRASH_DIR, row["filename"])
        if os.path.exists(src):
            shutil.move(src, dst)
            moved_count += 1
    
    print(f"  移动了 {moved_count} 个文件到 {TRASH_DIR}")
    print(f"  剩余 {len(df) - len(to_move)} 个合格模型")
    
    # 统计原因
    n_dist = (df["min_dist"] > MIN_DIST_THRESHOLD).sum()
    n_isc = (df["I_sc"] >= ISC_THRESHOLD).sum()
    n_both = ((df["min_dist"] > MIN_DIST_THRESHOLD) & (df["I_sc"] >= ISC_THRESHOLD)).sum()
    print(f"    - min_dist > {MIN_DIST_THRESHOLD}: {n_dist} 个")
    print(f"    - I_sc >= {ISC_THRESHOLD}: {n_isc} 个")
    print(f"    - 两者同时满足: {n_both} 个")


if __name__ == "__main__":
    main()
