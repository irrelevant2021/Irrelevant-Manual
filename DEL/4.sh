#!/bin/bash
for target in $(ls -d DEL_*/ | tr -d '/')
do
        echo $target
	cd $target
	mv whole.csv whole_$target.csv
	mv output.csv top10000_$target.csv
	cd ..
done


