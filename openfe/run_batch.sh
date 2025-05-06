#!/bin/bash
list1=(
rbfe_DY-0137_complex_DY-0139_complex.json
rbfe_DY-0137_complex_DY-0141_complex.json
rbfe_DY-0137_solvent_DY-0139_solvent.json
rbfe_DY-0137_solvent_DY-0141_solvent.json
)
list2=(
rbfe_DY-0138_complex_DY-0142_complex.json
rbfe_DY-0138_solvent_DY-0142_solvent.json
rbfe_DY-0142_complex_DY-0140_complex.json
rbfe_DY-0142_solvent_DY-0140_solvent.json
)
list3=(
rbfe_DY-0143_complex_DY-0137_complex.json
rbfe_DY-0143_solvent_DY-0137_solvent.json
rbfe_DY-0144_complex_DY-0142_complex.json
rbfe_DY-0144_solvent_DY-0142_solvent.json
)
list4=(
rbfe_DY-0145_complex_DY-0151_complex.json
rbfe_DY-0145_solvent_DY-0151_solvent.json
rbfe_DY-0146_complex_DY-0148_complex.json
rbfe_DY-0146_solvent_DY-0148_solvent.json
)
list5=(
rbfe_DY-0147_complex_DY-0151_complex.json
rbfe_DY-0147_solvent_DY-0151_solvent.json
rbfe_DY-0149_complex_DY-0151_complex.json
rbfe_DY-0149_solvent_DY-0151_solvent.json
)
list6=(
rbfe_DY-0150_complex_DY-0148_complex.json
rbfe_DY-0150_solvent_DY-0148_solvent.json
rbfe_DY-0152_complex_DY-0148_complex.json
rbfe_DY-0152_solvent_DY-0148_solvent.json
)

(
echo "working on gpu device 1 ..."
for file in ${list1[@]}
do
        dirpath=${file%.*}  # strip off final ".json"
	echo $dirpath
	sed -i 's/"gpu_device_index": null/"gpu_device_index": [1]/g' ./transformations/$file
	openfe quickrun transformations/$file -o results/$file -d results/$dirpath
done
) > 1.log 2>&1 < /dev/null &

(
echo "working on gpu device 2 ..."
for file in ${list2[@]}
do
        dirpath=${file%.*}  # strip off final ".json"
        echo $dirpath
        sed -i 's/"gpu_device_index": null/"gpu_device_index": [2]/g' ./transformations/$file
        openfe quickrun transformations/$file -o results/$file -d results/$dirpath
done
) > 2.log 2>&1 < /dev/null &

(
echo "working on gpu device 3 ..."
for file in ${list3[@]}
do
        dirpath=${file%.*}  # strip off final ".json"
        echo $dirpath
        sed -i 's/"gpu_device_index": null/"gpu_device_index": [3]/g' ./transformations/$file
        openfe quickrun transformations/$file -o results/$file -d results/$dirpath
done
) > 3.log 2>&1 < /dev/null &

(
echo "working on gpu device 4 ..."
for file in ${list4[@]}
do
        dirpath=${file%.*}  # strip off final ".json"
        echo $dirpath
        sed -i 's/"gpu_device_index": null/"gpu_device_index": [4]/g' ./transformations/$file
        openfe quickrun transformations/$file -o results/$file -d results/$dirpath
done
) > 4.log 2>&1 < /dev/null &

(
echo "working on gpu device 5 ..."
for file in ${list5[@]}
do
        dirpath=${file%.*}  # strip off final ".json"
        echo $dirpath
        sed -i 's/"gpu_device_index": null/"gpu_device_index": [5]/g' ./transformations/$file
        openfe quickrun transformations/$file -o results/$file -d results/$dirpath
done
) > 5.log 2>&1 < /dev/null &

(
echo "working on gpu device 6 ..."
for file in ${list6[@]}
do
        dirpath=${file%.*}  # strip off final ".json"
        echo $dirpath
        sed -i 's/"gpu_device_index": null/"gpu_device_index": [6]/g' ./transformations/$file
        openfe quickrun transformations/$file -o results/$file -d results/$dirpath
done
) > 6.log 2>&1 < /dev/null &
