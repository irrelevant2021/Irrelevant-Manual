# Install
Install reinvent4. https://github.com/MolecularAI/REINVENT4 # some big package in requirements-linux-64-2.lock using `pip install <package_name> -i https://pypi.tuna.tsinghua.edu.cn/simple` or download by vpn before pip install.  
Install maize. https://github.com/MolecularAI/REINVENT4   
install gnina. https://github.com/gnina/gnina

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

filter.* for deal with the results of RL.  
`python filter.py` to get top molecules' smiles  
`bash filter.sh` to find their docking pose  
#filter_sdf.txt in filter.sh is from Sublime after grep  


### By this way, generating a batch of pose-realistically predictable ligands, which go through a series of screenings(MW, logP, SAcore et al.), in the target binding pocket, now evaluate them with MM/GBSA!  

## note  
cuda in conda is dfferent with /usr/local/cuda-x(install by .run file, which should set environment variables by yourself), while gnina seems have to use cuda in /usr/local/cuda-x. I still haven't figured out the relationship between their adaptations, and here are some commands I used:  
remove them  
`sudo apt-get --purge remove "*cublas*" "*cudnn*" "*cuda*"` or run cuda-uninstaller in /usr/local/cuda-12.0/bin (installed by .run file)  
`sudo apt-get remove --purge '^nvidia-.*'`  
`sudo apt-get autoremove`  
install them  
from ubuntu software updater and https://developer.nvidia.com/cuda-toolkit (don't chose driver at cuda installtion, because it has been have installed)   
set environment variables   
`export LD_LIBRARY_PATH=/usr/local/cuda-x/lib64:$LD_LIBRARY_PATH`   
#`export PATH=/usr/local/cuda-x/bin:$PATH`  
#`export CUDA_HOME=/usr/local/cuda-x`

