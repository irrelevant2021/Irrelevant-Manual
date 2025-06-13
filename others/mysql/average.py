import argparse
parser = argparse.ArgumentParser(description="average the numbers in cell")
parser.add_argument("-i", "--input", type=str, help="input")
parser.add_argument("-o", "--output", type=str, help="output")
args = parser.parse_args()

i = args.input
o = args.output

import pandas as pd
import re
import numpy as np

def average_numbers_in_cell(cell):
    if pd.isna(cell):
        return cell
    # 提取所有数字（支持负号、小数点）
    numbers = re.findall(r'-?\d+\.?\d*', str(cell)) #从开头匹配正则表达式，匹配成功后，匹配区域不再参与下一次匹配
    if len(numbers) <= 0:
        return cell  # 单个数字或无数字，保持原样
    nums = list(map(float, numbers))
    return np.mean(nums)

def process_csv(input_path, output_path):
    df = pd.read_csv(input_path, dtype=str)  # 读为字符串避免类型问题
    for col in df.columns:
        if col in ("External_Code", "COMPOUNDS", "Batch", "Structure"):
            continue
        df[col] = df[col].apply(average_numbers_in_cell)
    df.to_csv(output_path, index=False)

process_csv(i, o)
