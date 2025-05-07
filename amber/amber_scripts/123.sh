#!/bin/bash

ScriptPath='/data/Irrelevant-Manual/amber/amber_scripts'
WorkPath='\/data\/AKT\/amber\/8uw9\/20250506'
num_gpu=8
current_list_index=0
echo "ligand,dG,rms" > mmpbsa_result.csv

#Leap
for file in *.pdb
#for file in $(ls *.pdb | grep -E 'CHEMBL5172982|CHEMBL5175396|CHEMBL5185309|CHEMBL5187179|CHEMBL5189173|CHEMBL5191796|CHEMBL5204246|CHEMBL5205764|CHEMBL5206121')
do
        echo ${file%.pdb}
        mkdir ${file%.pdb}
        cd ${file%.pdb}
        mv ../$file .
        cp ../../protein.pdb .
        cp $ScriptPath/leap* .
        $SCHRODINGER/utilities/prepwizard -WAIT -noepik -noprotassign -nometaltreat $file $file
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

#protocol & multi-GPU & run
total=$(ls -d */ | tr -d '/' | wc -l)
per=$((total / num_gpu))
echo "------"
echo "total tasks: $total"
echo "gpu num setting: $num_gpu"
echo "running $per tasks per gpu"
echo "results write in mmpbsa_result.csv"
echo "------"

for subdir in $(ls -d */ | tr -d '/' )
do
        list_index=$((current_list_index % num_gpu))
	#if $var match x), then execute
        case $list_index in
		0) list0+=("$subdir") ;;
	        1) list1+=("$subdir") ;;
        	2) list2+=("$subdir") ;;
	        3) list3+=("$subdir") ;;
        	4) list4+=("$subdir") ;;
	        5) list5+=("$subdir") ;;
        	6) list6+=("$subdir") ;;
	        7) list7+=("$subdir") ;;
        esac
        current_list_index=$((current_list_index + 1))
done

#for i in ${list0[@]};do echo $i;done
device=$(seq 0 $((num_gpu - 1)))
echo "using device: $device"
for i in $device
do
	(
        echo "working on gpu device $i ..."
        declare -n current_list="list$i"
        for dir in "${current_list[@]}"
        do
                echo $dir
		cd $dir
                cp $ScriptPath/MD_amber_openmm.py ./MD_amber_openmm.py
                sed -i 's/EditWorkDir/'$WorkPath'\/'$dir'/g'  ./MD_amber_openmm.py
                sed -i 's/EditEqTime/1/g' ./MD_amber_openmm.py
                sed -i 's/EditProdTime/10/g' ./MD_amber_openmm.py
                sed -i 's/EditProdNum/2/g' ./MD_amber_openmm.py
                sed -i 's/EditProdSkip/1/g' ./MD_amber_openmm.py
                sed -i 's/WhichCUDA/'$i'/g' ./MD_amber_openmm.py
		python MD_amber_openmm.py
		mkdir mmpbsa
		cd mmpbsa
	        cp $ScriptPath/mmgbsa.in .
        	cp $ScriptPath/decom.py .
	        ante-MMPBSA.py  -p ../SYS_gaff2.prmtop -c com.prmtop -r rec.prmtop -l ligand.prmtop -s :WAT:Na+:Cl-:Mg+:K+ -n :LIG --radii mbondi2
        	mpirun -np 10 MMPBSA.py.MPI -O -i mmgbsa.in -o FINAL_RESULTS_MMGBSA.dat -sp ../SYS_gaff2.prmtop -cp com.prmtop -rp rec.prmtop -lp ligand.prmtop -y ../prot_lig_prod1-2_whole.dcd
	        result=$(grep 'DELTA TOTAL' FINAL_RESULTS_MMGBSA.dat | awk '{print $3}')
        	rms=$(awk -F ',' '{sum += $2} END {print sum/NR}' ../rmsd_lig.csv)
	        echo -e "$dir,\c" >> ../../mmpbsa_result.csv
        	echo -e "$result,\c" >> ../../mmpbsa_result.csv
	        echo -e "$rms" >> ../../mmpbsa_result.csv
        	python decom.py
	        cd ..
		cd ..
        done
	) > $i.log 2>&1 < /dev/null &
done
