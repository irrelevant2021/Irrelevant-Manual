#!/bin/bash

num_gpu=6
current_list_index=0

#multi-GPU & run
total=$(ls ls ./transformations/*json | wc -l)
per=$((total / num_gpu))
echo "------"
echo "total tasks: $total"
echo "gpu num setting: $num_gpu"
echo "running $per tasks per gpu"
echo "------"

for subdir in $(ls ./transformations/*json)
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
#	(
        echo "working on gpu device $i ..."
	#declare that "list$i" is a $var
        declare -n current_list="list$i"
        for json in "${current_list[@]}"
        do
                echo $json
		file=${json:18}  # strip off "./transformations/"
		dirpath=${file%.*}  # strip off final ".json"
		echo $file
		echo $dirpath
        	sed -i 's/"gpu_device_index": null/"gpu_device_index": ['$i']/g' $json
#	        openfe quickrun $json -o results/$file -d results/$dirpath
        done
#	) > $i.log 2>&1 < /dev/null &
done
