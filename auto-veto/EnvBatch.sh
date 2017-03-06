#!/bin/bash
# Create an environment just for batch jobs.
# Copied mostly from SetupEnv.sh
# Clint Wiseman, USC

if [ "${CHOS}" != "sl64" ]; then
	env CHOS=sl64 chos
fi

export HOMEDIR=/global/homes/w/wisecg
export SWDIR=/project/projectdirs/majorana/software/sl64
export CXXFLAGS='-std=c++11'
export PATH=${HOMEDIR}:${PATH}	# rmate
export PATH=${HOMEDIR}/channelSel/cookies:${PATH} # MkCookie

alias rootmj="root -b -l ~/env/rootlogon.C"
alias qloginmj="qlogin -l h_vmem=2G -q mndl_prod.q"
alias qsubmj="qsub -l h_vmem=2G -q mndl_prod.q"

# MJ software version
# Check for new versions in $SWDIR/mjsw
# If you want to run a manual version:
#   1. Uncomment MANUALSW   2. Comment out SETUP_MJSW****
export MJSWDIR=${SWDIR}/mjsw/mjsw201612Prod
export SETUP_MJSW201612PROD="clint is great"
# export MANUALSW="clint is great"

source /etc/profile.d/modules.sh
source /common/majorana/scripts/setupMajorana.sh
source /project/projectdirs/majorana/software/sl64/setupMajoranaSL64.sh

# MANUAL SETUP --------------------------------------------------

# Choose custom ($HOMEDIR) or default ($MJSWDIR) installs of MJSW.
if [ -n "$MANUALSW" ]; then

	# ROOT --------------------------------------------------
	# Look for MJD's current ROOT version in:
	# /global/project/projectdirs/majorana/software/sl64/root
	# export ROOTSYS=${SWDIR}/root/root-6.06.04
	export ROOTSYS=${SWDIR}/root/root-6.06.08_py2.7.11
	source ${ROOTSYS}/bin/thisroot.sh
	# -------------------------------------------------------

	# CLHEP -------------------------------------------------
	# Look for MJD's current CLHEP version in:
	# /global/project/projectdirs/majorana/software/sl64/CLHEP
	# export CLHEP_BASE_DIR ${SWDIR}/CLHEP/2.2.0.5/CLHEP
	# export CLHEP_BASE_DIR=${SWDIR}/CLHEP/2.3.2.2/CLHEP
	export CLHEP_BASE_DIR=${HOMEDIR}/mgsw/CLHEP/2.3.2.2/install
	export CLHEP_INCLUDE_DIR=${CLHEP_BASE_DIR}/include
	export CLHEP_LIB_DIR=${CLHEP_BASE_DIR}/lib
	export CLHEP_LIB=CLHEP
	export PATH=$CLHEP_BASE_DIR/bin:${PATH}
	export LD_LIBRARY_PATH=${CLHEP_LIB_DIR}:${LD_LIBRARY_PATH}
	# -------------------------------------------------------

	# MJSW --------------------------------------------------
  export MGDODIR=${HOMEDIR}/mgsw/MGDO
  export MJORDIR=${HOMEDIR}/mgsw/MJOR
	export GATDIR=${HOMEDIR}/mgsw/GAT
  export ORDIR=${HOMEDIR}/mgsw/OrcaRoot

  # export MGDODIR=${MJSWDIR}/MGDO
	# export GATDIR=${MJSWDIR}/GAT
  # export MJORDIR=${MJSWDIR}/MJOR
	# export ORDIR=${MJSWDIR}/OrcaRoot
  # export MAGEDIR=${MJSWDIR}/MaGe

  export TAMDIR=${MGDODIR}/tam
	export PATH=$MGDODIR/install/bin:$MJSWDIR/bin:$ORDIR/Applications:${MJORDIR}:${GATDIR}/Apps:${GATDIR}/Scripts:${PATH}
	export LD_LIBRARY_PATH=${MJSWDIR}/lib:${MGDODIR}/install/lib:${TAMDIR}/lib:${ORDIR}/lib:${GATDIR}/lib:${LD_LIBRARY_PATH}
	export ROOT_INCLUDE_PATH=$MGDODIR/Base:$MGDODIR/Gerda:$MGDODIR/GerdaTransforms:$MGDODIR/Majorana:$MGDODIR/MJDB:$MGDODIR/Root:$MGDODIR/Tabree:$MGDODIR/Tools:$MGDODIR/Transforms:$TAMDIR/inc:$MGDODIR/tam:$CLHEP_BASE_DIR/CLHEP/lib:${GATDIR}/BaseClasses:${GATDIR}/MGOutputMCRunProcessing
	# -------------------------------------------------------
fi