#!/bin/csh -f
#$ -cwd
#$ -j y
#$ -o /global/u1/w/wisecg/vetoScan/output
#$ -P majorana

source /etc/profile.d/modules.csh
source /common/majorana/scripts/setupMajorana.csh
source /global/u1/w/wisecg/ManualSetup.csh

# Example usage options
# 1. qsub RunScan.csh Your_File.txt
# 2. qsub -l debug=1 RunScan.csh Your_File.txt
# 3. qsub -l h_vmem=2G -q mndl_prod.q RunScan.csh

echo Running vetoScan with input file: 
echo $1

cd /global/u1/w/wisecg/vetoScan/
echo job start
date
make clean
make --quiet
./vetoScan ./runs/$1
echo job end
date
echo Why doesnt the script run the last line
