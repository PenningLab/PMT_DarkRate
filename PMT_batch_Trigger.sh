#!/bin/bash
if [ $# -lt 3 ]
then
 echo "Missing arguments!!! To run"
 echo "bash PMT_batch_Trigger.sh [number_of_sample] [dir] [last file #] <-i [startingindex]> <-b [baseline_samples]> <-T>"
 exit 1
fi
number_samples=$1
dirname=$2
num=$3
init=0
bsam=-1
pth=-1
trig=0
shift 3
while getopts 'i:b:p:T' opt; do
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
	        T )
		        trig=1
			echo "No trigger in use"
			;;

	esac
done

for j in `seq ${init} ${num}`
do
 cmdarr=(./DDC10_data_readout -wd ${dirname} -i ${j}_PMT.txt -o ${j}_PMT_Trigger.root -n ${number_samples} -invert -e)

 if [ "${bsam}" -ne -1 ]; then
  cmdarr+=(-bs ${bsam})
 fi
 if [ "${pth}" -ne -1 ]; then
  cmdarr+=(-pt ${pth})
 fi
 if [ "${trig}" -ne 1 ]; then
  cmdarr+=(-t ${j}_Trigger.txt)
 fi
 echo "${cmdarr[@]}"
 "${cmdarr[@]}"
 #echo "file ${j} complete"
done
