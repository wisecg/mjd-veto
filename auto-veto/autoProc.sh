#!/bin/bash

# run auto-veto over a run list, and optionally make me press a key
# to advance each run, so I can examine the errors.

# debug on login node
# cat runs/p3lqk-complete.txt | while read -r line; do echo $line; ./auto-veto $line; done

# interactive version
cat runs/test.txt | while read -r line;
do echo $line;
./auto-veto $line;
echo ""
read -p "Press Enter to continue..." </dev/tty
done

# grid version
# cat runs/p3lqk-complete.txt | while read -r line; do echo $line; qsub autoSubmit.sh $line; done