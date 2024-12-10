import pandas as pd
df = pd.read_csv('Tanimoto.csv')
df = df.sort_values(by='Tanimoto similarity ECF6 (raw)', ascending=False)
df[['SMILES', 'Comment', 'Tanimoto similarity ECF6 (raw)']].to_csv('Tanimoto.csv', index=None)
