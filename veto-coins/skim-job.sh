#!/bin/bash
#$ -cwd
#$ -j y
#$ -o /global/homes/w/wisecg/veto-coins/logs/
#$ -P majorana
source /global/homes/w/wisecg/env/EnvBatch.sh
cd /global/homes/w/wisecg/veto-coins

echo "Job Start:"
date
echo "Node:  "$HOSTNAME
echo "Job ID:  "$JOB_ID

echo "./skim-coins $1 $2"
./skim-coins $1 $2 ./data

echo "Job Complete:"
date
