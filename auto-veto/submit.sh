#! /bin/bash
#$ -cwd
#$ -j y
#$ -o /global/homes/w/wisecg/vetoProcessing/logs
#$ -P majorana
source /global/homes/w/wisecg/env/EnvBatch.sh
cd /global/homes/w/wisecg/vetoProcessing/data

export dataPath=/project/projectdirs/majorana/data/mjd/surfmjd/data/raw/
# export dataDir=${dataPath}P3LQK/Data
# export dataDir=${dataPath}P3KJR/Data
# export dataDir=${dataPath}P3K93/Data
export dataDir=${dataPath}P3JDY/Data

echo "Job Start:"
date
echo "Node:  "$HOSTNAME
echo "Job ID:  "$JOB_ID
echo " "
if [ -n "$1" ]; then
   run=$1
else
   echo "Didn't get a run number.  Exiting ..."
   exit
fi
if [ -n "$2" ]; then
   opt=$2
fi
if [ ! "$opt" = "y" ]; then
   opt="N"
fi
echo "Checking run $run.  Force reprocess? $opt"
if [ -f OR_run${run}.root ]; then
   echo "Run has already been processed."
   if [ "$opt" = "y" ]; then
      echo "Force reprocess ..."
	  majorveroot --vetoloader $dataDir/*$run
   fi
else
   echo "Now processing run $run ..."
   majorveroot --vetoloader $dataDir/*$run
fi

echo "Job Complete:"
date
