import argparse 
parser = argparse.ArgumentParser(description="dcom_compare")
parser.add_argument("-ref", "--reference", type=str, help="ref")
parser.add_argument("-exp", "--expriment", type=str, help="exp")
parser.add_argument("-ene", "--energy", type=str, help="van der Waals, Electrostatic, Polar Solvation, Non-Polar Solv., TOTAL")
args = parser.parse_args()

import matplotlib.pyplot as plt
import pandas as pd

ref = args.reference
exp = args.expriment
#van der Waals, Electrostatic, Polar Solvation, Non-Polar Solv., TOTAL
energy = args.energy

df = pd.read_csv(ref + '/mmpbsa/FINAL_DECOMP_MMGBSA_deltas.csv', skiprows=2)
df1 = pd.read_csv(exp + '/mmpbsa/FINAL_DECOMP_MMGBSA_deltas.csv', skiprows=2)
ddf = df[energy].subtract(df1[energy])
ddf = pd.DataFrame(ddf)
ddf = pd.concat([df['Residue'], ddf], axis=1)

ddf['TOTAL_abs'] = ddf[energy].abs()
data1 = ddf.query('TOTAL_abs > 0.2')
#data1 = data1.query('TOTAL_abs < 3')
print(data1)

plt.bar(data1['Residue'], data1[energy])
plt.tick_params(axis='x', which='both', top=True, labelbottom=False, labeltop=True)
plt.xticks(rotation=90)
plt.title('dDECOMP_MMGBSA_ref' + ref + 'vs' + exp, fontsize=20)
plt.ylabel('ddG', fontsize=16)
#plt.ylim([0, -10])
plt.text(0.5, 0.9, 'ddG < 0 means the residue has better interaction with ref', fontsize=10, ha='center', va='center')
plt.text(0.5, 0.85, '(details in FINAL_DECOMP_MMGBSA_deltas.csv)', fontsize=10, ha='center', va='center')
#plt.savefig('dDECOMP_MMGBSA_ref' + ref + 'vs' + exp, bbox_inches='tight')


