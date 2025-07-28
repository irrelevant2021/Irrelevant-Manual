#!/bin/bash

for i in $(seq 2 39)
do
	nam=$(awk -F "," "NR == $i { print \$1 }" 20250612.csv)
	smi=$(awk -F "," "NR == $i { print \$2 }" 20250612.csv)
	echo $nam
	echo $smi
	cp ./example.yaml ./$nam.yaml
	sed -i 's|EDIT_SMILES|'$smi'|g' ./$nam.yaml
done

for file in $(ls *.yaml | grep -E 'DY')
do
	boltz predict $file --use_msa_server --cache /data/.boltz/
done

echo "ligand,boltz_value,boltz_probability" > boltz_result.csv

for ii in $(ls DY*.yaml)
do
        echo ${ii%.yaml}
        result1=$(awk -F ": " "NR == 2 { print \$2 }" boltz_results_${ii%.yaml}/predictions/*/affinity_*.json)
        result2=$(awk -F ": " "NR == 3 { print \$2 }" boltz_results_${ii%.yaml}/predictions/*/affinity_*.json)
        echo -e "${ii%.yaml},\c" >> ./boltz_result.csv
        echo -e "$result1\c" >> ./boltz_result.csv
        echo -e "${result2%,}" >> ./boltz_result.csv
done
