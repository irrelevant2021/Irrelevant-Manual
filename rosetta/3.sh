#!/bin/bash
sed '/LG1/!s/HETATM/ATOM  /g' ternary0.pdb > ternary_0.pdb
#pdb_list.txt:
#ternary_0.pdb
python ../../Ternary_modifier/ternary_modify.py -s pdb_list.txt
$ROSETTA/main/source/scripts/python/public/molfile_to_params.py ternary_0_lig.mol2 -n TRN
