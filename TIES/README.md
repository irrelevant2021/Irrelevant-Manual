https://github.com/UCL-CCS/TIES  
https://github.com/UCL-CCS/TIES_MD  
  
## prepare  
after align ligands and save to .pdb files (including the protein file) at here  
`bash prepare.sh`  
`ties -action create -l ligand1.pdb ligand2.pdb -lff gaff2 -nc 0 -p protein.pdb -md openmm`  
change your filename, check which atoms turn up (1),and which turn off (-1) at 'ties/ties-ligand1-ligand2/lig/build/complex.pdb'  
  
## modification  
`bash mod.sh`  
modify running time et al. see mod.sh  

## run  
`ties_md --config_file=./ties/ties-ligand1-ligand2/com/TIES.cfg`   
`ties_md --config_file=./ties/ties-ligand1-ligand2/lig/TIES.cfg`   
change your filename, watch the progress at 'ties/ties-GSK843-IN05019/XXX/LAMBDA_X.XX/repX/X/X.log'  

## analysis  
`ties_ana --run_type=setup`   
write 'legs = com, lig' in analysis.cfg, and chose the pairs of ligands in exp.dat,then  
`ties_ana`  
