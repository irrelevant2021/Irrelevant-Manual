#!/bin/bash
for target in $(ls -d DEL_*/ | tr -d '/')
do
        echo $target
	cd $target
	#sed -i 's/Open//g' whole.csv
	python  main2.py  
	cd ..
done


