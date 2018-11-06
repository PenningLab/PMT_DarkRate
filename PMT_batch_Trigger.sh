#!/bin/sh
if [ $# -lt 2 ]
then
 echo "Missing arguments!!! To run"
 echo "bash One_trigger_PMT_batch.sh [number_of_sample] [dir] [# files] <starting index>"
 return 
fi
number_samples=$1
dirname=$2
num=`expr $3 - 1`
init=0
if [ $# -eq 4 ]
then
 init=$4
 num=`expr ${init} + ${num}`
fi

for j in `seq ${init} ${num}`
do
 #echo "starting file ${j}"
 echo "./DDC10_data_readout -wd ${dirname} -i ${j}.txt -o ${j}_${outfilename}.root"
 ./DDC10_data_readout_trigger -wd ${dirname} -i ${j}_PMT.txt -t ${j}_Trigger.txt -o ${j}_PMT_Trigger.root -n ${number_samples}
 echo "file ${j} complete"
done

