Install reinvent4. https://github.com/MolecularAI/REINVENT4 #some big package in requirements-linux-64-2.lock using `pip install <package_name> -i https://pypi.tuna.tsinghua.edu.cn/simple` or download by vpn before pip install.  
Install maize. https://github.com/MolecularAI/REINVENT4 #

TL:  
At 'smi', dwnload csv file from ChemBL, and get compounds.smi by `df[['Smiles', 'pChEMBL Value', 'Molecule ChEMBL ID']].to_csv('compounds.smi', header=None, index=False, sep='\t')`   
At 'TL',  run `reinvent transfer_learning.toml` , and check running state by `tensorboard --bind_all --logdir tb_TL`.  
  
RL:  
(change py)   
At 'maize-RL'
