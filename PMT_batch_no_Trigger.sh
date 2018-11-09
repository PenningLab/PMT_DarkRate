#!/bin/bash
if [ $# -lt 3 ]
then
 echo "Missing arguments!!! To run"
 echo "bash PMT_batch_no_Trigger.sh [number_of_sample] [working dir] [last file #] <-i [startingindex]> <-b [baseline_samples]> <-p [pulseThreshold]>"
 exit 1
fi
number_samples=$1
dirname=$2
num=$3
init=0
bs=-1
pth=-1
shift 3
while getopts 'i:b:p:' opt; do
	case "${opt}" in
		i )
			init=${OPTARG}
			echo "initial file is ${init}"
			;;
		b )
			bsam=${OPTARG}
			echo "using ${bsam} baseline samples"
			;;
		p )
			pth=${OPTARG}
			echo "pulsethreshold is ${pth}"
			;;

	esac
done

for j in `seq ${init} ${num}`
do
 #echo "starting file ${j}"
 cmdarr=(./DDC10_data_readout -wd ${dirname} -i ${j}_Trigger.txt -o ${j}_PMT_Only.root -n ${number_samples} -e)

 if [ "${bsam}" -ne -1 ]; then
  cmdarr+=(-bs ${bsam})
 fi
 if [ "${pth}" -ne -1 ]; then
  cmdarr+=(-pt ${pth})
 fi
 echo "${cmdarr[@]}"
 "${cmdarr[@]}"
 #echo "file ${j} complete"
done
