#!/bin/bash
cp ./sdf/*.mol2 .
$ROSETTA/main/source/scripts/python/public/molfile_to_params.py CBN.mol2 -n CBN -p CBN --conformers-in-one-file --chain=X
$ROSETTA/main/source/scripts/python/public/molfile_to_params.py POI.mol2 -n POI -p POI --conformers-in-one-file --chain=Y
sed -i '$ d' CBN.params
sed -i '$ d' POI.params
cat PP.pdb CBN.pdb POI.pdb > PPP.pdb
$ROSETTA_BIN/docking_prepack_protocol.mpi.linuxgccrelease -s PPP.pdb -use_input_sc -extra_res_fa CBN.params POI.params
mv PPP_0001.pdb TER.pdb
mpirun -np 12 $ROSETTA_BIN/docking_protocol.mpi.linuxgccrelease -s TER.pdb -nstruct 2 -use_input_sc -spin -dock_pert 3 8 -partners BX_CY -ex1 -ex2aro -load_PDB_components false -extra_res_fa CBN.params POI.params -out:file:scorefile score.sc -score:docking_interface_score
