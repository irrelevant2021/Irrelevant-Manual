import pandas as pd

df = pd.read_csv('output.csv')
df = pd.DataFrame(df, columns=['HitIndex', 'HitCount', 'r1_name', 'r2_name', 'r3_name', 'Smiles', 's_smiles', 'r1_smiles', 'r2_smiles', 'r3_smiles', 'BarCount', 'BarCode'])

df.to_csv('final.csv', index=False)
