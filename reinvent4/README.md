# Install
Install reinvent4. https://github.com/MolecularAI/REINVENT4 # some big package in requirements-linux-64-2.lock using `pip install <package_name> -i https://pypi.tuna.tsinghua.edu.cn/simple` or download by vpn before pip install.  
Install maize. https://github.com/MolecularAI/REINVENT4   
install gnina. https://github.com/gnina/gnina  
install QSARtuna. https://github.com/MolecularAI/QSARtuna

# TL:  
At 'smi', dwnload csv file from ChemBL, then get compounds.smi and validation_compounds.smi by   
`df['Smiles'].to_csv('temp.csv', header=None, index=False, sep='\t')`    
`sed -i '/\./d' temp.csv` #reinvnet do not recognize salt   
`shuf -n 10 temp.csv > validation_compounds.smi`  
`grep -v -F -f validation_compounds.smi temp.csv > compounds.smi`    

At 'TL',  run `reinvent transfer_learning.toml` , and check running state by `tensorboard --bind_all --logdir tb_TL`.   
#lots ref molecules in mol2mol model will explode memory, <= 25.   
#batch_size are roughly between 50 and 200   
#I would choose the lowest point of loss_Validation Loss   
  
# RL:  

At 'maize-RL', run `reinvent staged_learning.toml` , and check running state by `tensorboard --bind_all --logdir tb_logs_0`  
#staged_learning.toml-[$PATH_to_executable/maize, gnina.yml, maize.toml];  
#gnina.yml-[config of maize for gnina];  
#maize.toml-[$PATH_to_executable/gnina, $PATH_to_executable/temp]  
#receptor.pdb should be changed to your file, and search_center should be changed to your number in gnina.yml  
#PRIMARY_SCORE_TAG = "minimizedAffinity" should be changed to "CNNscore" in $YOUR_CONDA_PATH/envs/maize/lib/python3.10/site-packages/maize/steps/mai/docking/gnina.py    
#[CNN_VS is CNNaffinity*CNN_score; CNNaffinity is a binding affinity prediction in pK units (e.g. 9 is 1 nM); CNNscore is a predicted probability that the pose is good (<2A from correct); minimizedAffinity is the Vina score]  

filter.* for deal with the results of RL.  
`python filter.py` to get top molecules' smiles  
`bash filter.sh` to find their docking pose  
#filter_sdf.txt in filter.sh is from Sublime after grep  


#### By this way, generating a batch of pose-realistically predictable ligands, which go through a series of screenings(MW, logP, SAcore et al.), in the target binding pocket.   
#### Base on this, larger molecules are developed(Tanimoto), or molecules with stronger affinity(MM/GBSA-QSARtuna), and so on.  

## note  
cuda in conda(auto) is dfferent with /usr/local/cuda-x(install by .run file, which should set environment variables by yourself), while gnina seems have to use cuda in /usr/local/cuda-x. I still haven't figured out the relationship between their adaptations, and here are some commands I used:  
1. install cuda with .run file at /usr/local/cuda-x from https://developer.nvidia.com/cuda-toolkit (don't chose driver at cuda installtion, because it has been have installed)   
2. set environment variables temporarily, ONLY at reinvnet4! It affects the variables of CUDA in other environments!!!   
`export LD_LIBRARY_PATH=/usr/local/cuda-x/lib64:$LD_LIBRARY_PATH`   

3. others:  
  remove them  
  `sudo apt-get --purge remove "*cublas*" "*cudnn*" "*cuda*"` or run cuda-uninstaller in /usr/local/cuda-12.0/bin (installed by .run file)  
  `sudo apt-get remove --purge '^nvidia-.*'`  
  `sudo apt-get autoremove`  
  install NVIDIA driver from ubuntu software updater is good. 
