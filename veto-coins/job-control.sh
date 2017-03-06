#!/bin/bash
# C. Wiseman, USC.

function runThresh()
{
  # DS map (from DataSetInfo.hh): {0,76},{1,51},{2,7},{3,24},{4,22},{5,80}
  dsNum=3
  dsMax=24
  echo "running MkCookie"
  MkCookie mjdb.phy.ornl.gov 443 majorana Pump\$\&G\@ug3\$
  echo "make -s"
  make -s
  for ((i=0;i<=${dsMax};i++)); do
    echo $i
    # ./skim-job.sh ${dsNum} $i
    qsub skim-job.sh ${dsNum} $i
  done
}

function mergeFiles()
{
  # DS map (from DataSetInfo.hh): {0,76},{1,51},{2,7},{3,24},{4,22},{5,80}
  dsNum=5
  dsMax=80
  outFile="skimDS${dsNum}_HitsOver2630.root"
  args=()
  for ((i=0; i<=${dsMax}; i++)); do
      args+=("skimDS${dsNum}_$i.root")
  done
  echo "${args[@]}"
  thisDir=`pwd`
  cd data
  hadd $outFile ${args[@]}
  cd $thisDir
}

function findExposure()
{
  arr=(0 1 3 4 5)
  for i in "${arr[@]}"; do
    echo $i
    make -s && ./thresholds -e final/thresholdsDS$i.root $i
  done
}

function applyThresholds()
{
  make -s
  # DS map (from DataSetInfo.hh): {0,76},{1,51},{2,7},{3,24},{4,22},{5,80}
  dsNum=3
  # dsMax=6
  for ((i=19; i<=24; i++)); do
    ./skim_mjd_data ${dsNum} $i ./skim -cgw final/thresholdsDS${dsNum}.root -t 0.9
  done
}

# =================================
# runThresh
mergeFiles
# findExposure
# applyThresholds