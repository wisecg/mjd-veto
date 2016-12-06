#!/bin/bash
# run auto-veto over a run list, and optionally make me press a key
# to advance each run, so I can examine the errors.
# C. Wiseman, USC.

# NOTE: Make sure your run lists have a newline at the end!
# NOTE: To change the output directory, must edit auto-job.sh or auto-multijob.sh

make -s

# 1. debug - login node
# cat runs/test.txt | while read -r line; do echo $line; ./auto-veto $line; done

# 2. interactive version - login node
# cat runs/test.txt | while read -r line;
# do echo $line;
# ./auto-veto $line;
# echo ""
# read -p "Press Enter to continue..." </dev/tty
# done

# 3. grid version - qsub the jobs
# neat bash trick: see how many are still running:
# qstat -u wisecg | wc -l
# cat runs/ds1-complete.txt | while read -r line; do echo $line; qsub auto-job.sh $line; done

# 4. login node mode - save the output to logfiles
# cat runs/ds1-complete.txt | while read -r line; do echo $line; ./auto-veto $line > ./logs/auto-job.sh.o$line; done

# 5. grid multi-job - Scan 300 runs per job
mapfile RunList < ./runs/ds5-incomplete.txt
i=0
runs=${#RunList[@]}
tempArray=
while ((i < ${#RunList[@]})); do
  tempArray+=" "${RunList[$i]};
  if [[ $i%300 -eq 0 ]]; then
    qsub auto-multijob.sh $tempArray
    tempArray=
  fi
  ((i++))
done
qsub auto-multijob.sh $tempArray

# 6. clean up
# cat runs/ds5-incomplete.txt | while read -r line; do mv ./avout/veto_run$line.root ./avout/DS5/; done