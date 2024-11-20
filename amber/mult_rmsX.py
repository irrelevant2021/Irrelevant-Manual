import pandas as pd  
import matplotlib.pyplot as plt  
import os

#folders = [f for f in os.listdir() if os.path.isdir(f)] 
folders = ['20241115-1', '20241115-3']


for folder in folders:
    if os.path.isdir(folder): 
        folder_path = os.path.join(os.getcwd(), folder)
        file_path = os.path.join(folder_path, "rmsd_lig.csv")
        data = pd.read_csv(file_path)
        x = data["Unnamed: 0"].values
        y = data["0"].values
        plt.plot(x, y, label=folder)
plt.legend()
plt.xlabel("Frame", fontsize = 14, fontweight = 'bold')
plt.ylabel("RMSD_LIG [$\AA$]", fontsize = 14, fontweight = 'bold')
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
#plt.savefig('RMSD_LIG_compare.png')
plt.show()

for folder in folders:
    if os.path.isdir(folder):
        folder_path = os.path.join(os.getcwd(), folder)
        file_path = os.path.join(folder_path, "rmsd_ca.csv")
        data = pd.read_csv(file_path)
        x = data["Unnamed: 0"].values
        y = data["0"].values
        plt.plot(x, y, label=folder)
plt.legend()
plt.xlabel("Frame", fontsize = 14, fontweight = 'bold')
plt.ylabel("RMSD_CA [$\AA$]", fontsize = 14, fontweight = 'bold')
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
#plt.savefig('RMSD_CA_compare.png')
plt.show()

for folder in folders:
    if os.path.isdir(folder):
        folder_path = os.path.join(os.getcwd(), folder)
        file_path = os.path.join(folder_path, "rmsf_ca.csv")
        data = pd.read_csv(file_path)
        x = data["Unnamed: 0"].values
        y = data["1"].values
        plt.plot(x, y, label=folder)
plt.legend()
plt.xlabel("Residue", fontsize = 14, fontweight = 'bold')
plt.ylabel("RMSF ($\AA$)", fontsize = 14, fontweight = 'bold')
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
#plt.savefig('RMSF_CA_compare.png')
plt.show()
