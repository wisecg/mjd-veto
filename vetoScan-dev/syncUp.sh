#!/bin/bash

# This uses atom's "script" package.
# Execute with cmd-i.

echo "Syncing vetoScan to PDSF"
rsync -av --progress --exclude='*.DS_Store' --exclude='*sublime*' --exclude='*.o' --exclude='*.root' ./ wisecg@dtn01.nersc.gov:/global/u1/w/wisecg/vetoScan