#!/bin/bash
# This is a handy little tool for searching multiple patterns in the production
# log files, and figuring out the run they come from.
# C. Wiseman, 11/3/2016

# export logdir="/project/projectdirs/majorana/data/production/mjdProcessingLogs"
export logdir="./logs"

# grep -H: give the filename
# grep -h: don't give the filename

# grep -rIl "Error\[18\]:" $logdir | xargs grep -hn "processing run\|Error\[18\]:"
grep -rIl "Error\[18\]:" $logdir | xargs grep -hn "processing run"

# grep -rIl "segmentation" $logdir | xargs grep -Hn "processing run"

# grep -rIl "Warning: Bad Duration" $logdir | xargs grep -Hn "processing run"

# grep -rIl "LED may be" $logdir | xargs grep -hn "processing run\|LED may be"

# how many hits?
# hitType="2+ panels"
# hitType="vertical"
# hitType="side+bottom"
# hitType="top+sides"
# hitType="compound"
# grep -rIl "$hitType" $logdir | xargs grep -hn "$hitType" | wc -l
