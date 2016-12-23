#! /bin/bash
#$ -cwd
#$ -j y
#$ -o /global/homes/w/wisecg/auto-veto/logs/
#$ -P majorana
source /global/homes/w/wisecg/env/EnvBatch.sh
cd /global/homes/w/wisecg/auto-veto

echo "Job Start:"
date
echo "Node:  "$HOSTNAME
echo "Job ID:  "$JOB_ID
echo " "
echo "Auto-multijob got this many runs: "$#

for var in "$@"
do
    echo "Now processing run $var ..."
    # ./auto-veto $var -o avout/DS5/
    ./auto-veto $var
done

echo "Job Complete:"
date
