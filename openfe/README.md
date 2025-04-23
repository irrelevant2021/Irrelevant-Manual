[https://github.com/OpenFreeEnergy/openfe]  

`openfe plan-rbfe-network -M tyk2_ligands.sdf -p tyk2_protein.pdb -o network_setup`  
`openfe view-ligand-network network_setup/ligand_network.graphml`  
`openfe quickrun path/to/transformation.json -o results.json -d working-directory`  
`openfe gather results/ --report dg -o final_results.tsv`  
