#!/bin/bash
cd preprocess
for target in $(ls *.trimed.fq)
do
	echo ${target%.trimed.fq}
	awk 'NR%4==2' $target | awk '{print substr($1,1,13)","substr($1,14,13)","substr($1,27,13)","substr($1,40,9)}' > ../${target%.trimed.fq}.csv
	grep "^[ATCG]" ../${target%.trimed.fq}.csv | wc -l
	wc -l ../${target%.trimed.fq}.csv
	#grep "^[ATCG]" $target.trimed.fq | awk '{print substr($1,1,13)","substr($1,14,13)","substr($1,27,13)","substr($1,40,9)}' > ../$target.csv
done
cd ..
