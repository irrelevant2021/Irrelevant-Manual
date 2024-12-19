#!/bin/bash
python ../../ternary_model_prediction.py -d TER_0001.pdb -l LNK.pdb -da decoy_atom_list.txt -la linker_atom_list.txt -c 0.4 -wd decoy_atom_list_delete.txt -ld linker_atom_list_delete.txt -t default -r rmsd.txt
