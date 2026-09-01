#!/bin/bash

ScriptPath='/data/Irrelevant-Manual/amber/amber_scripts'
WorkPath='\/data\/AR\/amber\/ROSETTA\/20260827'
num_gpu=4
current_list_index=0
echo "ligand,dG,rms" > mmpbsa_result.csv

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
        list_index=$((current_list_index % num_gpu + 4))
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
device=$(seq 4 $((num_gpu + 3)))
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
#                cp $ScriptPath/MD_amber_openmm.py ./MD_amber_openmm.py
		cp $ScriptPath/analysis.py ./analysis.py
		sed -i 's/EditWorkDir/'$WorkPath'\/'$dir'/g'  ./analysis.py
                sed -i 's/EditEqTime/1/g' ./analysis.py
                sed -i 's/EditProdTime/10/g' ./analysis.py
                sed -i 's/EditProdNum/10/g' ./analysis.py
                sed -i 's/EditProdSkip/1/g' ./analysis.py
#                sed -i 's/WhichCUDA/'$i'/g' ./analysis.py
		python analysis.py
		mkdir mmpbsa
		cd mmpbsa
	        cp ../../../mmgbsa.in .
        	cp $ScriptPath/decom.py .
	        ante-MMPBSA.py  -p ../SYS_gaff2.prmtop -c com.prmtop -r rec.prmtop -l ligand.prmtop -s :WAT:Na+:Cl-:Mg+:K+ -n :LIG --radii mbondi2
        	mpirun -np 10 MMPBSA.py.MPI -O -i mmgbsa.in -o FINAL_RESULTS_MMGBSA.dat -sp ../SYS_gaff2.prmtop -cp com.prmtop -rp rec.prmtop -lp ligand.prmtop -y ../prot_lig_prod1-10_whole.dcd
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
