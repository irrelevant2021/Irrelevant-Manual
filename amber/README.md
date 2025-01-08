THANKS TO https://github.com/pablo-arantes/making-it-rain `Arantes P.R., Depólo Polêto M., Pedebos C., Ligabue-Braun R. Making it rain: cloud-based molecular simulations for everyone. Journal of Chemical Information and Modeling 2021. DOI: 10.1021/acs.jcim.1c00998.`

## amber  
### amber_scripts  
all scripts put here, copy them to your working_dir for MD
### PROTEIN_ID  
protein.pdb #your protein file, after `pdbfixer your_protein.pdb --output=protein.pdb --keep-heterogens=none --add-residues --replace-nonstandard --add-atoms=none`;  

LIGANDS_ID #including mult-your_ligands.pdb, here is working_dir, start MD with `cp ../../amber_scripts/123.sh .; bash 123.sh`  


## note  
1. `conda install mamba -c conda-forge`, mamba is great!  
2. you can change your MD parameters at 123.sh(mainly), mmgbsa.in and MD_amber_openmm.py.  
3. check your version of `nvidia-smi` and `nvcc -V`, before`conda install -c conda-forge openmm cudatoolkit=1x.x` would be useful!   
4. `python -m openmm.testInstallation`  
