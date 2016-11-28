#! /bin/bash
#$ -cwd
#$ -j y
#$ -o /global/homes/w/wisecg/auto-veto/logs
#$ -P majorana
source /global/homes/w/wisecg/env/EnvBatch.sh
cd /global/homes/w/wisecg/auto-veto

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

echo "Now processing run $run ..."
./auto-veto $run

echo "Job Complete:"
date
