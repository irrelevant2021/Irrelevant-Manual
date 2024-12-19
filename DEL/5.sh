#!/bin/bash
for target in $(ls -d DEL_*/ | tr -d '/')
do
        echo $target
	cd $target
	cp ../main3.py .
	python main3.py
	mv final.csv final_$target.csv
	cd ..
done


