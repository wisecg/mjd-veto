#! /bin/bash
#$ -cwd
#$ -j y
#$ -o /global/homes/w/wisecg/auto-veto/skimout
#$ -P majorana
source /global/homes/w/wisecg/env/EnvBatch.sh
cd /global/homes/w/wisecg/auto-veto

# suggestion:  This scans ge data for a long time,
# so try submitting with:
# qsub -l h_vmem=2G skim-job.sh

echo "Job Start:"
date
echo "Node:  "$HOSTNAME
echo "Job ID:  "$JOB_ID
echo " "

make -s && ./skim-coins -l 5 runs/ds5-complete.txt skimout/

echo "Job Complete:"
date
