#!/bin/bash
if [ $# -lt 3 ]
then
 echo "Missing arguments!!! To run"
 echo "bash PMT_batch_Trigger.sh [waveform channel] [dir] [last file #] <-i [startingindex]> <-b [baseline_samples]> <-T>"
 exit 1
fi
wformchan=$1
dirname=$2
num=$3
init=0
bsam=-1
pth=-1
trig=-1
shift 3
while getopts 'i:b:p:T:' opt; do
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
			trig=${OPTARG}
			echo "Using channel ${trig} as trigger"
			;;

	esac
done

for j in `seq ${init} ${num}`
do
 cmdarr=(./DDC10_bin_data_readout -wd ${dirname} -i ${j}.bin -o ${j}_PMT_Trigger.root -wform ${wformchan} -invert -e)

 if [ "${bsam}" -ne -1 ]; then
  cmdarr+=(-bs ${bsam})
 fi
 if [ "${pth}" -ne -1 ]; then
  cmdarr+=(-pt ${pth})
 fi
 if [ "${trig}" -ne -1 ]; then
  cmdarr+=(-t ${trig})
 fi
 echo "${cmdarr[@]}"
 "${cmdarr[@]}"
 #echo "file ${j} complete"
done
