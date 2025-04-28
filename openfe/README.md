[https://github.com/OpenFreeEnergy/openfe]  
  
`pdbfixer tyk2_mae.pdb --output=tyk2_protein.pdb --keep-heterogens=none --add-residues --replace-nonstandard --add-atoms=H`  
`openfe plan-rbfe-network -M tyk2_ligands.sdf -p tyk2_protein.pdb -o network_setup` OR `python transformations.py`  
`openfe view-ligand-network network_setup/ligand_network.graphml`  
`openfe quickrun path/to/transformation.json -o results.json -d working-directory`  OR `bash run.py`  
`openfe gather results/ --report dg -o final_results.tsv`  
