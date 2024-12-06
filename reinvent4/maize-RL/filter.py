import pandas as pd
df = pd.read_csv('staged_learning_1.csv')

df = df.query('Score >= 0.6')
df = df.query('`docking with MAIZE (raw)` >= 0.9')
df = df.query('`QED (raw)` >= 0.6')
df = df.query('`Number of rotatable bonds (raw)` <= 5')
df = df.query('`Number of stereo centers (raw)` <= 1')
df = df.query('`SlogP (RDKit) (raw)` <= 3')
df = df.query('`SA score (raw)` <= 4')

df['Smiles'].to_csv('filter.csv', header=None, index=False, sep='\t')
