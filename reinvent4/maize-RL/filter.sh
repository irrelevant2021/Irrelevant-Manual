
#!/bin/bash
mkdir filter
for i in $(ls *sdf)
do
	echo $i
	awk '/\$\$\$\$/ {n++; next} {print > "filter/" substr(FILENAME, 1, length(FILENAME)-3) n ".sdf"}' "$i"
done

#filter_smiles.csv is from df['SMILES'] of filter.csv 
#for i in $(cat filter_smiles.csv) ;do grep -x -F $i filter/*sdf; done

#filter_sdf.txt is from Sublime after grep
#for i in $(cat filter_sdf.txt) ;do cat $i >> filter/result.sdf; echo '$$$$' >> filter/result.sdf; done
