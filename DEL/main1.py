import pandas as pd

#读取原始RCODE表，及大小
raw = pd.read_csv('raw.csv', header=None)
raw.columns = ['R1', 'R2', 'R3', 'scafford']
print('raw shape =')
rate1 = raw.shape[0]  
print(rate1)

#读取scafford参照表
s_lib = pd.read_csv('../../UPDEL/format_codon/meta.txt', delimiter=' ', header=None)
s_lib.columns = ['s_name', 's_smiles', 'scafford']

#合并原始表与scafford参照，记s_raw
raw = pd.merge(left=raw, right=s_lib, on='scafford', how='inner')
print('scaf match =')
rate2 = raw.shape[0]
print(rate2)
RATE = rate2 / rate1
print('perfect match:')
print(RATE)

result_data = []

#切割raw，分别在对应骨架下匹配R
grouped = raw.groupby('s_name')
raw = None

for name, group in grouped:
	print('in' + name + '...')

    #在raw的00XX中，读00XX的R参照表
	lib_add = '../../UPDEL/format_codon/' + name + '.txt'
	r_lib = pd.read_csv(lib_add,  delimiter='\t', header=None)
	r_lib.columns = ['r_name', 'R', 'r_smiles', 'BARCODE']

	# 使用 duplicated 函数判断 'BARCODE' 中是否有重复的值，没有判断scaf
	r_lib['mrak'] = r_lib['BARCODE'].duplicated()
	#按R123，切割R参照表，赋存在r_lib_grouped_data中
	r_lib = r_lib.groupby('R')
	r_lib_grouped_data = {}
	for r_lib_name, r_lib_group in r_lib:
		r_num = str(r_lib_name)
		r_lib_grouped_data[r_num] = r_lib_group
	r_lib = None
		
		
	#开始匹配
	r1s_raw = pd.merge(left=group, right=r_lib_grouped_data['1'], left_on='R1', right_on='BARCODE', how='inner')
	r2s_raw = pd.merge(left=r1s_raw, right=r_lib_grouped_data['2'], left_on='R2', right_on='BARCODE', how='inner')
	r3s_raw = pd.merge(left=r2s_raw, right=r_lib_grouped_data['3'], left_on='R3', right_on='BARCODE', how='inner')
	print('R1,R2,R3 match =')
	print(r3s_raw.shape)

	#删除多余列  
	r3s_raw = r3s_raw.drop('BARCODE_x', axis=1)
	r3s_raw = r3s_raw.drop('BARCODE_y', axis=1)
	r3s_raw = r3s_raw.drop('BARCODE', axis=1)
	r3s_raw = r3s_raw.drop('R_x', axis=1)
	r3s_raw = r3s_raw.drop('R_y', axis=1)
	r3s_raw = r3s_raw.drop('R', axis=1)
	#更新表头
	r3s_raw.columns = ['R1', 'R2', 'R3', 'scafford', 's_name', 's_smiles', 'r1_name', 'r1_smiles', 'mark1', 'r2_name', 'r2_smiles', 'mark2', 'r3_name', 'r3_smiles', 'mark3']
	r3s_raw['HitIndex'] = r3s_raw[['s_name', 'r1_name', 'r2_name', 'r3_name']].apply(lambda x: '-'.join(x.astype(str)), axis=1)
	r3s_raw['BarCode'] = r3s_raw[['R1', 'R2', 'R3', 'scafford']].apply(lambda y: '-'.join(y.astype(str)), axis=1)
	r3s_raw['HitIndex'] = r3s_raw['HitIndex'] #去掉了Open前缀,因为现在的SMILES align ref没有Open前缀,20240522


	#lib_smiles_add = '../Raw_Lib/SMEILESlab/Open' + name + '_Fully_Enumerated_Structures.txt'
	#smiles_lib = pd.read_csv(lib_smiles_add, delimiter='\t', header=None)
	#smiles_lib.columns = ['HitIndex', 'Smiles']
	#smiles_raw = pd.merge(left=r3s_raw, right=smiles_lib, on='HitIndex', how='left')      #merge base on r3s_raw, using df[df.isnull().any(axis=1)] to find data which can not find smiles

	#result_data.append(smiles_raw)
	result_data.append(r3s_raw)
#for test
	#test = str(name) 
	#if test == "DEL0001":
	#	break
group = None
r_lib_grouped_data = None
r3s_raw = None

result = pd.concat(result_data, ignore_index=True)
result_data = None

#处理总结果
print('Counting Hit...')
hitcount = result['HitIndex'].value_counts()
result['HitCount'] = result['HitIndex'].map(hitcount)
barcount = result['BarCode'].value_counts()
result['BarCount'] = result['BarCode'].map(barcount)
print('total match =')
print(result.shape)

#result.sort_values(by='HitCount', ascending=False)
print('problem data...')
print(result['mark1'].value_counts())
print(result['mark2'].value_counts())
print(result['mark3'].value_counts())

result.to_csv('whole.csv', index=False)
result.drop_duplicates(subset='HitIndex').to_csv('summary.csv')


