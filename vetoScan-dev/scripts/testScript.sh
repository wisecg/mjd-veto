#!/bin/csh -f
#$ -cwd
#$ -j y
#$ -o /global/u1/w/wisecg/vetoScan/output
#$ -P majorana

echo The $0 command is called with $#argv parameters
echo   parameter 1 is $1
echo   parameter 2 is $2
echo   parameter 3 is $3
echo   parameter 4 is $4
echo All Parameters are \"$argv\"
echo 2nd and on parameters are \"$argv[2-]\" 

#setenv ROOTSYS /global/project/projectdirs/majorana/software/sl64/root/root_v5.34.34
#cd $ROOTSYS
#source bin/thisroot.csh

source /etc/profile.d/modules.csh
source /common/majorana/scripts/setupMajorana.csh
source /global/u1/w/wisecg/ManualSetup.csh

echo ./vetoScan ./runs/$1.txt

#setenv

echo $CHOS
echo $ROOTSYS
echo $MGDODIR
echo $MJDDATADIR
echo $MJSWDIR

which root
ldd $ROOTSYS/bin/root.exe
file $ROOTSYS/bin/root.exe
file $ROOTSYS/lib/libCore.so

cd $HOME/vetoScan
echo HEY CLINT
echo
root-config --ldflags

mjdiskspace

./vetoScan ./runs/Debug.txt

echo HAY CLETUS
