#!/bin/csh
#$ -cwd
#$ -j y
#$ -o /global/u1/w/wisecg/vetoScan/output
#$ -P majorana

cd /global/u1/w/wisecg/vetoScan/
echo job start
date
./vetoScan ./runs/JDY_01.txt
echo job end
date
