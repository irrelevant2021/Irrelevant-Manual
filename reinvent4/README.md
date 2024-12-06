Install reinvent4. https://github.com/MolecularAI/REINVENT4 #some big package in requirements-linux-64-2.lock using `pip install <package_name> -i https://pypi.tuna.tsinghua.edu.cn/simple` or download by vpn before pip install.  
Install maize. https://github.com/MolecularAI/REINVENT4 #

TL:  
At 'smi', dwnload csv file from ChemBL, then get compounds.smi and validation_compounds.smi by   
`df[['Smiles', 'pChEMBL Value', 'Molecule ChEMBL ID']].to_csv('temp.csv', header=None, index=False, sep='\t')`    
`sed -i '/\./d' temp.csv` #reinvnet do not recognize salt   
`shuf -n 10 temp.csv > validation_compounds.smi`  
`grep -v -F -f validation_compounds.smi temp.csv > compounds.smi`    

At 'TL',  run `reinvent transfer_learning.toml` , and check running state by `tensorboard --bind_all --logdir tb_TL`.   
#lots ref molecules in mol2mol model will explode memory, <= 25.   
#batch_size are roughly between 50 and 200
  
RL:  
(change py)   
At 'maize-RL'
