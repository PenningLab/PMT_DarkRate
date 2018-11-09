#!/bin/bash
if [ $# -lt 3 ]
then
 echo "Missing arguments!!! To run"
 echo "bash PMT_batch_no_Trigger.sh [number_of_sample] [working dir] [last file #] <-i [startingindex]> <-b [baseline_samples]> <-p [pulseThreshold]>"
 return
fi
number_samples=$1
dirname=$2
num=$3
init=0
bs=-1
pth=-1
while getopts ":i:b:p:" do
	case ${opt} in
		i )
			init=$OPTARG
			;;
		b )
			bs=$OPTARG
			;;
		p )
			pth=$OPTARG
done

for j in `seq ${init} ${num}`
do
 #echo "starting file ${j}"
 cmdarr=(./DDC10_data_readout -wd ${dirname} -i ${j}_Trigger.txt -o ${j}_PMT_Only.root -n ${number_samples} -e)

 if [ ${bs} -ne -1 ]; then
  cmdarr+=(-b ${bs})
 fi
 if [ ${pth} -ne -1 ]; then
  cmdarr+=(-pt ${pth})
 fi
 echo "${cmdarr[@]}"
 "${cmdarr[@]}"
 #echo "file ${j} complete"
done
