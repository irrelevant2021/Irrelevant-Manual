import glob
import pandas as pd
import argparse
parser = argparse.ArgumentParser(description="Multiple molecules, having the same prefix, like 'pre-*', and references compare the differences in residue interaction.")
parser.add_argument("-ref", "--reference", type=str, help="ref")
parser.add_argument("-pre", "--prefix", type=str, help="pre")
args = parser.parse_args()


ref = args.reference
pre = args.prefix
energy = 'TOTAL'

df = pd.read_csv(ref + '/mmpbsa/FINAL_DECOMP_MMGBSA_deltas.csv', skiprows=2)
table = pd.DataFrame()
table = pd.concat([table, df['Residue']], axis=1)

for exp in glob.glob(pre):
    df1 = pd.read_csv(exp + '/mmpbsa/FINAL_DECOMP_MMGBSA_deltas.csv', skiprows=2)
    ddf = df[energy].subtract(df1[energy])
    ddf = pd.DataFrame(ddf)
    ddf = pd.concat([df['Residue'], ddf], axis=1)
    ddf['TOTAL_abs'] = ddf[energy].abs()
    ddf.columns = ['Residue', exp, exp + '_abs']
    ddf = ddf.drop(columns='Residue')
    table = pd.concat([table, ddf], axis=1)

f0 = table.filter(regex='_abs$')  #输出一张abs的表
f1 = f0 > 1  #判断abs列的cell是否大于1,是为true,输出一张f1表
f2 = f1.any(axis=1)  #判断f1中的每一行里是否有ture(即用或筛选),输出一列f2数据
table = table[f2]  #输出table中f2为true的表
table = table.drop(columns=f0.columns)  #删除abs列
table.to_csv('ddcomp_summary.csv', index=False)
print(table)
