#!/bin/bash



. /cvmfs/hadessoft.gsi.de/install/debian10/6.24.02/hydra2-6.3/defall.sh
. /cvmfs/hadessoft.gsi.de/install/debian10/root-6.24.02/bin/thisroot.sh


export HADDIR=/cvmfs/hadessoft.gsi.de/install/debian10/6.24.02/hydra2-6.3


#export CERN_ROOT=/cvmfs/hades.gsi.de/install/cernlib_gfortran/2005/
#export HGEANT_DIR=/lustre/nyx/hades/user/rlalik/fwdet/install/hgeant2-fwdet

#export PLUTODIR=/cvmfs/hades.gsi.de/install/5.34.34/pluto/v5_44/

export LC_ALL=C
export ROOTLOGON=$HADDIR/macros/rootlogon.C

export MYHADDIR=

#INSTALL_DIR=/lustre/nyx/hades/user/rlalik/hades/install/5.34.34
#export PATH=${INSTALL_DIR}/bin:${MYHADDIR}/bin:${PATH}
export LD_LIBRARY_PATH=${INSTALL_DIR}/lib:${MYHADDIR}/lib:${PLUTODIR}:${HADDIR}/lib:${LD_LIBRARY_PATH}
