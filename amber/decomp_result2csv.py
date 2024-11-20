import matplotlib.pyplot as plt
import pandas as pd

with open('FINAL_DECOMP_MMPBSA_bak.dat', 'r') as f:
	dat = f.read()

lines = dat.split('\n\n')


for i, line in enumerate(lines):
	if 'DELTAS' in line:
		start = i
		break


deltas = ''.join(lines[start:start+1])


with open('FINAL_DECOMP_MMPBSA_deltas.csv', 'w') as f:
	f.write(deltas)

data = pd.read_csv('FINAL_DECOMP_MMPBSA_deltas.csv', skiprows=2)
data = data.dropna(subset=['Residue'])

data['TOTAL_abs']=data['TOTAL'].abs()
data1 = data.query('TOTAL_abs > 0.5')
data1 = data1.query('TOTAL_abs < 20')

plt.bar(data1['Residue'], data1['TOTAL'])
plt.tick_params(axis='x', which='both', top=True, labelbottom=False, labeltop=True)
plt.xticks(rotation=90)
plt.title('DECOMP_MMGBSA', fontsize=20)
plt.ylabel('dG', fontsize=16)
#plt.ylim([0, -10])
plt.savefig('FINAL_DECOMP_MMPBSA.png', bbox_inches='tight')

#plt.bar(data['ligand'], data['dG'])
#plt.tick_params(axis='x', which='both', top=True, labelbottom=False, labeltop=True)
#plt.xticks(rotation=80)
#plt.title('MM/GBSA', fontsize=20)
#plt.ylabel('dG', fontsize=16)
#plt.ylim([-38, -20])
#plt.savefig('CRBN_binder_MMGBSA.png', dpi=600, bbox_inches='tight')
