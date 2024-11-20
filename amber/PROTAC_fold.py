import argparse
import pandas as pd
import pytraj as pt
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="does PROTAC fold?")
parser.add_argument("-cpd", "--compound", type=str, help="compoundID")
args = parser.parse_args()
cpd = args.compound

dcd = cpd + '.dcd'
top = cpd + '.prmtop'
traj = pt.load(dcd, top)

rmsd = pt.rmsd(traj, ref = 0, mask = ":LIG")
plt.plot(rmsd)
plt.savefig(cpd + 'rmsd_lig.png')
plt.clf()

surf = pt.surf(traj, '(:LIG@O*,:LIG@N*)')
plt.hist(surf, density=True)
plt.savefig(cpd + 'surf.png')
plt.clf()

R = pt.radgyr(traj, ':LIG')
plt.hist(R, density=True)
plt.savefig(cpd + 'R.png')
plt.clf()

IMHB = pt.search_hbonds(traj, ':LIG')
df = pd.DataFrame(IMHB.values).T
plt.plot(df)
plt.savefig(cpd + 'IMHB.png')
plt.clf()

cluster_data = pt.cluster.kmeans(traj, mask=':LIG', n_clusters=5, maxit=1000, metric='dme')
df_cluster = pd.DataFrame(cluster_data.cluster_index)
df_cluster.columns = ['col']
plt.plot(df_cluster)
plt.savefig(cpd + 'cluster.png')
plt.clf()
print(cluster_data.summary())
print(cluster_data.population)
traj0 = pt.Trajectory(top=traj.top)
traj1 = pt.Trajectory(top=traj.top)
traj2 = pt.Trajectory(top=traj.top)
traj3 = pt.Trajectory(top=traj.top)
traj4 = pt.Trajectory(top=traj.top)
df0 = df_cluster.query('col == 0')
df1 = df_cluster.query('col == 1')
df2 = df_cluster.query('col == 2')
df3 = df_cluster.query('col == 3')
df4 = df_cluster.query('col == 4')
traj0.append(traj[df0.index])
traj1.append(traj[df1.index])
traj2.append(traj[df2.index])
traj3.append(traj[df3.index])
traj4.append(traj[df4.index])
traj0.save(cpd + 'cluster0.dcd')
traj1.save(cpd + 'cluster1.dcd')
traj2.save(cpd + 'cluster2.dcd')
traj3.save(cpd + 'cluster3.dcd')
traj4.save(cpd + 'cluster4.dcd')
