#!/bin/bash
if [ $# -lt 3 ]
then
 echo "Missing arguments!!! To run"
 echo "bash PMT_bin_batch.sh [waveform channel] [dir] [last file #] <-i [startingindex]> <-b [baseline_samples]> <-T [trigger channel]> <--pmt> <--win> <--bfile>"
 exit 1
fi
wformchan=$1
dirname=$2
num=$3
init=0
bsam=-1
pth=-1
trig=-1
invert=-1
win=-1
usebase=-1
it=-1
outname="PMT_Trigger.root"
shift 3
TEMP=`getopt -o i:b:b:p:T:o: --long pmt,tri,sit:,win:,bfile: -n 'PMT_bin_batch.sh' -- "$@"`
eval set -- "$TEMP"
while true; do
	case "$1" in
		-i )
			init=$2
			echo "initial file is ${init}"
			shift 2
			;;
		-b )
			bsam=$2
			echo "using ${bsam} baseline samples"
			shift 2
			;;
		-p )
			pth=$2
			echo "pulsethreshold is ${pth}"
			shift 2
			;;
		-T )
			trig=$2
			echo "Using channel ${trig} as trigger"
			shift 2
			;;
	    -o )
		    outname="${2}.root"
			echo "Using output name ${outname}"
			shift 2
			;;
		--pmt )
			invert=0
			echo "using pmt as trigger, will not invert waveform"
			shift
			;;
		--win )
			win=$2
			echo "using window size ${win}"
			shift
			;;
		--bfile )
			usebase=0
			base_filename=$2
			echo "using baseline file ${base_filename}"
			shift
			;;
		--sit )
			it=$2
			echo "using baseline file ${base_filename}"
			shift 2
			;;
		-- )
			shift ;
			break
			;;
		* )
			shift
			;;
	esac
done

for j in `seq ${init} ${num}`
do
 cmdarr=(./DDC10_bin_data_readout -wd ${dirname} -i ${j}.bin -o ${j}_${outname} -wform ${wformchan} -e -tri)

 if [ "${bsam}" -ne -1 ]; then
  cmdarr+=(-bs ${bsam})
 fi
 if [ "${pth}" -ne -1 ]; then
  cmdarr+=(-pt ${pth})
 fi
 if [ "${trig}" -ne -1 ]; then
  cmdarr+=(-t ${trig})
 fi
 if [ "${win}" -ne -1 ]; then
  cmdarr+=(-win ${win})
 fi
 if [ "${usebase}" -ne -1 ]; then
  cmdarr+=(-bf ${base_filename})
 fi
 if [ "${invert}" -eq -1 ]; then
  cmdarr+=(-invert)
 fi
 if [ "${it}" -ne -1 ]; then
  cmdarr+=(-sit ${it})
 fi
 echo "${cmdarr[@]}"
 "${cmdarr[@]}"
 #echo "file ${j} complete"
done
