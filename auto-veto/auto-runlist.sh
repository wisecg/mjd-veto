#!/bin/bash
# run auto-veto over a run list, and optionally make me press a key
# to advance each run, so I can examine the errors.
# C. Wiseman, USC.

# 1. debug - login node
cat runs/test.txt | while read -r line; do echo $line; ./auto-veto $line; done

# 2. interactive version - login node
# cat runs/test2.txt | while read -r line;
# do echo $line;
# make --quiet && ./auto-veto $line;
# echo ""
# read -p "Press Enter to continue..." </dev/tty
# done

# 3. grid version - qsub the jobs
# neat bash trick: see how many are still running:
# qstat -u wisecg | wc -l
# cat runs/ds3-complete.txt | while read -r line; do echo $line; qsub auto-job.sh $line; done
