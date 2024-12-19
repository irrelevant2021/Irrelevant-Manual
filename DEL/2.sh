#!/bin/bash
mkdir output
for target in $(ls *.csv)
do
        echo ${target%.csv}
	mkdir DEL_${target%.csv}
	mv $target DEL_${target%.csv}
	cp ./*.py DEL_${target%.csv}
	cd DEL_${target%.csv}
	mv $target raw.csv
	python -u main1.py > log1.txt 2>&1 
	cp summary.csv ../output/$target 
	cd ..
done
