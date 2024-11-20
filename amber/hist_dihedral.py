#parm SYS_gaff2.prmtop
#trajin prot_lig_prod_2.dcd
#dihedral @4976 @4980 @4981 @4991 type psi range360 out 1.dat
#run

import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("1.dat", delimiter="    ")
plt.hist(df['  Dih_00001'], bins=30, density=True)
plt.xlim([0, 360])
plt.savefig('dihedral.png', bbox_inches='tight')
