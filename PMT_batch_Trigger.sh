#!/bin/bash
if [ $# -lt 3 ]
then
 echo "Missing arguments!!! To run"
 echo "bash PMT_batch_Trigger.sh [number_of_sample] [dir] [last file #] <-i [startingindex]> <-b [baseline_samples]>"
 exit 1
fi
number_samples=$1
dirname=$2
num=$3
init=0
bs=-1
while getopts ":i:b:e"
do
	case ${opt} in
		i )
			init=$OPTARG
			;;
		b )
			bs=$OPTARG
			;;
	esac
done

for j in `seq ${init} ${num}`
do
 cmdarr=(./DDC10_data_readout -wd ${dirname} -i ${j}_PMT.txt -t ${j}_Trigger.txt -o ${j}_PMT_Trigger.root -n ${number_samples})

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
