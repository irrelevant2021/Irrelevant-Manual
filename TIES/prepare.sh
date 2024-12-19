#!/bin/bash
mkdir prepare
for i in $(ls *.pdb)
do
        mv $i prepare
        cd prepare
        pdb4amber -i $i -o ligand_h_$i.pdb
        mv ligand_h_$i.pdb ../$i
        cd ..
done

pdbfixer myprotein.pdb --output=protein.pdb --keep-heterogens=none --add-residues --replace-nonstandard --add-atoms=none
