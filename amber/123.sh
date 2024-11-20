#!/bin/bash

ScriptPath=/media/star/UProtech1/zhennan/amber_scripts/
WorkPath=\/media\/star\/UProtech1\/zhennan\/RIPK2-0815\/amber\/RIPK2\/CHEMBL5136571

for file in $(ls *.pdb | grep -E 'CHEMBL5172982|CHEMBL5175396|CHEMBL5185309|CHEMBL5187179|CHEMBL5189173|CHEMBL5191796|CHEMBL5204246|CHEMBL5205764|CHEMBL5206121')
do
        echo ${file%.pdb}
        mkdir ${file%.pdb}
        cd ${file%.pdb}
        mv ../$file .
        cp ../../protein.pdb .
        cp $ScriptPath/leap* .
        pdb4amber -i $file -o ligand_h.pdb
        antechamber -i ligand_h.pdb -fi pdb -o ligand.mol2 -fo mol2 -c bcc -nc 0 -rn LIG -at gaff2
        parmchk2 -i ligand.mol2 -f mol2 -o ligand.frcmod -s gaff2
        tleap -f leap_lig.in
        #mv ./*_protein.pdb protein_mae.pdb
        #pdbfixer protein_mae.pdb --output=protein.pdb --keep-heterogens=all --replace-nonstandard --add-atoms=none 
        sed -i '/END/d' protein.pdb
        cat protein.pdb ligand_gaff.pdb > protein_ligand.pdb
        tleap -f leap_noi.in
        cd ..
done

for file in *.pdb
do
	echo ${file%.pdb}
	mkdir ${file%.pdb}
	cd ${file%.pdb}
	mv ../$file .
	cp ../../protein.pdb .
	cp $ScriptPath/leap* .
	pdb4amber -i $file -o ligand_h.pdb
	antechamber -i ligand_h.pdb -fi pdb -o ligand.mol2 -fo mol2 -c bcc -nc 1 -rn LIG -at gaff2
	parmchk2 -i ligand.mol2 -f mol2 -o ligand.frcmod -s gaff2
	tleap -f leap_lig.in
	#mv ./*_protein.pdb protein_mae.pdb
	#pdbfixer protein_mae.pdb --output=protein.pdb --keep-heterogens=all --replace-nonstandard --add-atoms=none 
	sed -i '/END/d' protein.pdb
	cat protein.pdb ligand_gaff.pdb > protein_ligand.pdb
	tleap -f leap_noi.in
	cd ..
done

echo "ligand,dG,rms" > mmpbsa_result.csv
for dir in $(ls -d */ | tr -d '/' )
#for dir in $(ls -d */ | grep -E -v '1024-02|1026-02|1030-11|1030-14' | tr -d '/')
do
	echo $dir
	cp $ScriptPath/MD_amber_openmm.py ./$dir/MD_amber_openmm.py
	sed -i 's/EditWorkDir/'$WorkPath'\/'$dir'/g'  ./$dir/MD_amber_openmm.py
	sed -i 's/EditEqTime/1/g' ./$dir/MD_amber_openmm.py
	sed -i 's/EditProdTime/10/g' ./$dir/MD_amber_openmm.py
	sed -i 's/EditProdNum/2/g' ./$dir/MD_amber_openmm.py
	sed -i 's/EditProdSkip/1/g' ./$dir/MD_amber_openmm.py
	cd ./$dir
	python MD_amber_openmm.py
	mkdir mmpbsa
        cd mmpbsa
        cp $ScriptPath/mmgbsa.in .
	cp $ScriptPath/decom.py .
        ante-MMPBSA.py  -p ../SYS_gaff2.prmtop -c com.prmtop -r rec.prmtop -l ligand.prmtop -s :WAT:Na+:Cl-:Mg+:K+ -n :LIG --radii mbondi2
        mpirun -np 50 MMPBSA.py.MPI -O -i mmgbsa.in -o FINAL_RESULTS_MMGBSA.dat -sp ../SYS_gaff2.prmtop -cp com.prmtop -rp rec.prmtop -lp ligand.prmtop -y ../prot_lig_prod_2.dcd
        result=$(grep 'DELTA TOTAL' FINAL_RESULTS_MMGBSA.dat | awk '{print $3}')
        rms=$(awk -F ',' '{sum += $2} END {print sum/NR}' ../rmsd_lig.csv)
        echo -e "$dir,\c" >> ../../mmpbsa_result.csv
        echo -e "$result,\c" >> ../../mmpbsa_result.csv
        echo -e "$rms" >> ../../mmpbsa_result.csv
	python decom.py
        cd ../../
done
