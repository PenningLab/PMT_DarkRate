# PMT_DarkRate


This repository contains codes that produce ready-to-analysis ntuple root file for PMT testing.

For typical use (random trigger on channel one and PMT signal on channel 3) :

  Copy txt files to PDSF /global/homes/w/wangryan/lz_project_data/wangbtc
  Then use One_trigger_PMT_batch.sh to extract info from txt files into root file for each run.
  (Usage :  bash One_trigger_PMT_batch.sh [number_of_sample] [dir] [# files] <starting index>)
  
  When it is done, use ana_combine_root.cc to combine all the files into one single root files
  (Usage : ./ana_combine_root -wd [working dir] -i [file name : PMT_Trigger.root] -o [output filename] -n [number of files to be combined])


For Triggering on PMT signal, use One_trigger_PMT_batch_PMT_Trigger.sh instead in teh first step and the rest is the same
