import pytraj as pt
import pandas as pd
import matplotlib.pyplot as plt

traj = pt.load('prot_lig_prod_1.dcd', 'SYS_gaff2.prmtop')
#sele0 = traj.top.select('(:LIG)')
#sele1 = traj.top.select('(:LIG)')
#sele0 = pt.center_of_mass(traj, mask=sele0)
dis = pt.distance(traj, ':9-381 :412-706')
df = pd.DataFrame(dis)

plt.plot(df.index, df[0]) 
plt.xlabel('frame')
plt.ylabel('distance')
plt.show()
