# PMT_DarkRate


This repository contains codes that produce ready-to-analysis ntuple root file for PMT testing.

For typical use (random trigger on channel one and PMT signal on channel 3) :

  Copy txt files to PDSF /global/homes/w/wangryan/lz_project_data/wangbtc
  Then use PMT_batch_Trigger.sh to extract info from txt files into root file for each run.
  (Usage :  bash PMT_batch_Trigger.sh [number_of_sample] [dir] [# files] <starting index>)
  
  When it is done, use ana_combine_root.cc to combine all the files into one single root files
  (Usage : ./ana_combine_root -wd [working dir] -i [file name : PMT_Trigger.root] -o [output filename] -n [number of files to be combined] -ct [upper charge threshold] <-t>)


For Triggering on PMT signal, use PMT_batch_PMT_Trigger.sh instead in the first step and the rest is the same
