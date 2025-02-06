import pandas as pd
import scipy.stats as stats

#B1-1是实验组第一次实验,B1-2是实验组第二次实验,B2-1是对照组第一次实验,B2-2是对照组第二次实验
df1_1 = pd.read_csv('B1-1_whole.csv')
df1_2 = pd.read_csv('B1-2_whole.csv')
df2_1 = pd.read_csv('B2-1_whole.csv')
df2_2 = pd.read_csv('B2-2_whole.csv')


#比较平行组情况
#对HitIndex计数,比较平行实验之间hit化合物的种类数量
count1_1 = df1_1[['HitIndex']].value_counts().reset_index()
count1_2 = df1_2[['HitIndex']].value_counts().reset_index()
count2_1 = df2_1[['HitIndex']].value_counts().reset_index()
count2_2 = df2_2[['HitIndex']].value_counts().reset_index()

B1 = pd.merge(left=count1_1, right=count1_2, left_on='HitIndex', right_on='HitIndex', how='inner')
print('B1, how many kinds of molecules hit after merged, exp1, and exp2')
print(B1.shape, count1_1.shape, count1_2.shape)
B2 = pd.merge(left=count2_1, right=count2_2, left_on='HitIndex', right_on='HitIndex', how='inner')
print('B2, how many kinds of molecules hit after merged, exp1, and exp2')
print(B2.shape, count2_1.shape, count2_2.shape)

#p值
t, p = stats.ttest_ind(B1['count_x'].head(50), B1['count_y'].head(50))
print('between experiments, top 50 HitIndex of B1 has p-value:')
print(p)
t = None
p = None
t, p = stats.ttest_ind(B2['count_x'].head(50), B2['count_y'].head(50))
print('between experiments, top 50 HitIndex of B2 has p-value:')
print(p)
t = None
p = None
print('以上, 平行实验之间是否有显著差异(top50)?')
print('   ')



#开始数据分析
#合并平行实验
df1 = pd.concat([df1_1, df1_2], axis=0)
df1 = df1.drop(['HitCount', 'BarCount', 'mark1', 'mark2', 'mark3'], axis=1)
print('Total number of B1 hits:')
print(df1.shape)
df2 = pd.concat([df2_1, df2_2], axis=0)
df2 = df2.drop(['HitCount', 'BarCount', 'mark1', 'mark2', 'mark3'], axis=1)
print('Total number of B2 hits:')
print(df2.shape)

#实验组-对照组
count2 = df2['HitIndex'].value_counts()
mask = df1.groupby('HitIndex').cumcount() >= df1['HitIndex'].map(count2).fillna(0)
df = df1[mask]
print('Total number of B1 - B2 hits:')
print(df.shape)

#添加R基团关联信息
df['r1r2_smiles'] = df[['r1_smiles', 'r2_smiles']].apply(lambda x: '-'.join(x.astype(str)), axis=1)
df['r1r3_smiles'] = df[['r1_smiles', 'r3_smiles']].apply(lambda x: '-'.join(x.astype(str)), axis=1)
df['r2r3_smiles'] = df[['r2_smiles', 'r3_smiles']].apply(lambda x: '-'.join(x.astype(str)), axis=1)
df['r1r2r3_smiles'] = df[['r1_smiles', 'r2_smiles', 'r3_smiles']].apply(lambda x: '-'.join(x.astype(str)), axis=1)

#分别计数
r1hit = df['r1_smiles'].value_counts()
r2hit = df['r2_smiles'].value_counts()
r3hit = df['r3_smiles'].value_counts()
r1r2hit = df['r1r2_smiles'].value_counts()
r1r3hit = df['r1r3_smiles'].value_counts()
r2r3hit = df['r2r3_smiles'].value_counts()
r1r2r3hit = df['r1r2r3_smiles'].value_counts()

#208256/x是平均水平,208256来自wc -l AllData.csv(df),x来自wc -l r*.csv
df['r1HitCount'] = df['r1_smiles'].map(r1hit) / 208256 * 30
df['r2HitCount'] = df['r2_smiles'].map(r2hit) / 208256 * 1262
df['r3HitCount'] = df['r3_smiles'].map(r3hit) / 208256 * 3534
df['r1r2HitCount'] = df['r1r2_smiles'].map(r1r2hit) / 208256 * 24294
df['r1r3HitCount'] = df['r1r3_smiles'].map(r1r3hit) / 208256 * 45100
df['r2r3HitCount'] = df['r2r3_smiles'].map(r2r3hit) / 208256 * 73218
df['r1r2r3HitCount'] = df['r1r2r3_smiles'].map(r1r2r3hit) / 208256 * 83138

#去重排序保存,打印top10
df.drop_duplicates(subset='r1_smiles').sort_values(by='r1HitCount', ascending=False).to_csv('r1.csv')
df.drop_duplicates(subset='r2_smiles').sort_values(by='r2HitCount', ascending=False).to_csv('r2.csv')
df.drop_duplicates(subset='r3_smiles').sort_values(by='r3HitCount', ascending=False).to_csv('r3.csv')
df.drop_duplicates(subset='r1r2_smiles').sort_values(by='r1r2HitCount', ascending=False).to_csv('r1r2.csv')
df.drop_duplicates(subset='r1r3_smiles').sort_values(by='r1r3HitCount', ascending=False).to_csv('r1r3.csv')
df.drop_duplicates(subset='r2r3_smiles').sort_values(by='r2r3HitCount', ascending=False).to_csv('r2r3.csv')
df.drop_duplicates(subset='r1r2r3_smiles').sort_values(by='r1r2r3HitCount', ascending=False).to_csv('r1r2r3.csv')
df.to_csv('AllData.csv')
pd.concat([df1[['r1_smiles', 'r1HitCount']].head(100), df2[['r2_smiles', 'r2HitCount']].head(100), df3[['r3_smiles', 'r3HitCount']].head(100), df12[['r1r2_smiles', 'r1r2HitCount']].head(100), df13[['r1r3_smiles', 'r1r3HitCount']].head(100), df23[['r2r3_smiles', 'r2r3HitCount']].head(100), df123[['r1r2r3_smiles', 'r1r2r3HitCount']].head(100)], axis=1).to_csv('r_top100.csv')

print('the R1 smiles hits number / the average of R1 smlies hits, or named Feature Intensity:')
print(df.drop_duplicates(subset='r1_smiles').sort_values(by='r1HitCount', ascending=False).head(10)[['r1_smiles', 'r1HitCount']])
print('the R2 smiles hits number / the average of R2 smlies hits, or named Feature Intensity:')
print(df.drop_duplicates(subset='r2_smiles').sort_values(by='r2HitCount', ascending=False).head(10)[['r2_smiles', 'r2HitCount']])
print('the R3 smiles hits number / the average of R3 smlies hits, or named Feature Intensity:')
print(df.drop_duplicates(subset='r3_smiles').sort_values(by='r3HitCount', ascending=False).head(10)[['r3_smiles', 'r3HitCount']])
print('the R1R2 smiles hits number / the average of R1R2 smlies hits, or named Feature Intensity:')
print(df.drop_duplicates(subset='r1r2_smiles').sort_values(by='r1r2HitCount', ascending=False).head(10)[['r1r2_smiles', 'r1r2HitCount']])
print('the R1R3 smiles hits number / the average of R1R3 smlies hits, or named Feature Intensity:')
print(df.drop_duplicates(subset='r1r3_smiles').sort_values(by='r1r3HitCount', ascending=False).head(10)[['r1r3_smiles', 'r1r3HitCount']])
print('the R2R3 smiles hits number / the average of R2R3 smlies hits, or named Feature Intensity:')
print(df.drop_duplicates(subset='r2r3_smiles').sort_values(by='r2r3HitCount', ascending=False).head(10)[['r2r3_smiles', 'r2r3HitCount']])
print('the R1R2R3 smiles hits number / the average of R1R2R3 smlies hits, or named Feature Intensity:')
print(df.drop_duplicates(subset='r1r2r3_smiles').sort_values(by='r1r2r3HitCount', ascending=False).head(10)[['r1r2r3_smiles', 'r1r2r3HitCount']])
print('以上,合并两次实验组,两次对照组,并作差后对smiles进行count')
