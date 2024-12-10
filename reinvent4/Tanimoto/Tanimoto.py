# Compare molecular similarity with Tanimoto to predict the likelihood of off-target.
# This script establishes a pathway for each molecule you pick out to compare with the active molecule from ChemBL.
# ChemBL683052.smi is too big to upload (copy at autoDL).
import os
import shutil
import pandas as pd

# final.csv is the final document after my training and MD.
df = pd.read_csv('final.csv')
df = df[['Title', 'smiles', 'ligand efficiency']]
df = df.sort_values(by='ligand efficiency', ascending=True)
df = df.query('`ligand efficiency` <= -1.4')

toml = 'scoring_components_example.toml'

for index, row in df.iterrows():
    dir_name = row['Title']
    replace_text = row['smiles']

    os.makedirs(dir_name, exist_ok=True)

    shutil.copy(toml, os.path.join(dir_name, 'scoring_components_example.toml'))

    file_path = os.path.join(dir_name, 'scoring_components_example.toml')
    with open(file_path, 'r') as file:
        content = file.read()

    content = content.replace('SMILES_EDIT', replace_text)

    with open(file_path, 'w') as file:
        file.write(content)

# bash command to run:
# for i in $(ls -d ./*/);do echo $i;cd $i;reinvent scoring_components_example.toml > /dev/null 2>&1 & cd ..;done
