import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv('mmpbsa_result.csv')
data['ligand'] = data['ligand'].astype(str)

plt.bar(data['ligand'], data['dG'])
plt.tick_params(axis='x', which='both', top=True, labelbottom=False, labeltop=True)
plt.xticks(rotation=90)
plt.title('MMGBSA', fontsize=20)
plt.ylabel('dG', fontsize=16)
plt.grid(True) 
plt.ylim([-50, -30])
plt.savefig('FINAL_MMPBSA.png', bbox_inches='tight')
