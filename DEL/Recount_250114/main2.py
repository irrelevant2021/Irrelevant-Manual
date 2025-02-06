import pandas as pd

#读取scafford参照表
s_lib = pd.read_csv('../UPDEL/format_codon/meta.txt', delimiter=' ', header=None)
s_lib.columns = ['s_name', 's_smiles', 'scafford']

#开始匹配Smiles
result = pd.read_csv('r1r2r3.csv')

output_data = []

s_lib.columns = ['name_s', 'smiles_s', 'scafford']
result = pd.merge(left=result, right=s_lib, on='scafford', how='inner')
result = result.drop('name_s', axis=1)
result = result.drop('smiles_s', axis=1)

result = result.groupby('s_name')

for s_name, s_group in result:
    print('in' + s_name + '...')
    #print('catch Smiles ...')
    lib_smiles_add = '../UPDEL/format_codon/' + s_name + '.tag.txt'
    smiles_lib = pd.read_csv(lib_smiles_add, delimiter='\t')
    smiles_lib.columns = ['HitIndex', 'Smiles']
    smiles_raw = pd.merge(left=s_group, right=smiles_lib, on='HitIndex', how='left')      #merge base on r3s_raw, using df[df.isnull().any(axis=1)] to find data which can not find smiles

    output_data.append(smiles_raw)
#for test
#    test = str(s_name)
#    if test == "DEL0001":
#        break

output = pd.concat(output_data, ignore_index=True)
output_data = None
output.to_csv('output.csv', index=False)



