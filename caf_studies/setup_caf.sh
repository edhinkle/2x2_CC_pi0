#!/usr/bin/env bash

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup edepsim v3_2_0c -q e20:prof

# Below is copied from NDCAFMaker setup script
setup cmake v3_22_2
setup gcc v9_3_0
setup pycurl
setup ifdhc
setup geant4 v4_11_0_p01c -q e20:debug
setup dk2nugenie   v01_10_01k -q debug:e20
setup genie_xsec   v3_04_00 -q AR2320i00000:e1000:k250
setup genie_phyopt v3_04_00 -q dkcharmtau
setup jobsub_client
setup eigen v3_3_5
setup srproxy v00.44 -q py3913 #Added 
setup duneanaobj v03_06_01b -q e20:prof
setup hdf5 v1_10_5a -q e20
setup fhiclcpp v4_15_03 -q debug:e20
#setup edepsim v3_2_0c -q debug:e20

# edep-sim needs to know where a certain GEANT .cmake file is...
G4_cmake_file=`find ${GEANT4_FQ_DIR}/lib64 -name 'Geant4Config.cmake'`
export Geant4_DIR=`dirname $G4_cmake_file`

# edep-sim needs to have the GEANT bin directory in the path
export PATH=$PATH:$GEANT4_FQ_DIR/bin

# shut up ROOT include errors
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$GENIE_INC/GENIE

# add SRProxy path to ROOT to allow for use in ROOT macros
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$SRPROXY_INC


# nusystematics paths
#export NUSYST=${PWD}/nusystematics
#export LD_LIBRARY_PATH=${NUSYST}/build/Linux/lib:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=${NUSYST}/build/nusystematics/artless:$LD_LIBRARY_PATH
#export FHICL_FILE_PATH=${NUSYST}/nusystematics/fcl:$FHICL_FILE_PATH

# Add pyGeoEff to pythonpath
export PYTHONPATH=${PYTHONPATH}:${PWD}/DUNE_ND_GeoEff/lib/

# duneananobj needs to be in the libs too
export LD_LIBRARY_PATH=${DUNEANAOBJ_LIB}:$LD_LIBRARY_PATH
# finally, add our lib & bin to the paths
mydir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export LD_LIBRARY_PATH=$mydir/lib:$LD_LIBRARY_PATH
export PATH=$mydir/bin:$PATH

# our FCL needs to be findable too
export FHICL_FILE_PATH="$FHICL_FILE_PATH:$mydir/cfg"

source ../2x2_sim/util/init.inc.sh

