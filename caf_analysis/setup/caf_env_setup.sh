#!/usr/bin/env bash

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup edepsim v3_2_0c -q e20:prof
setup duneanaobj v03_06_01b -q e20:prof
setup dk2nugenie   v01_10_01k -q debug:e20
setup cmake v3_27_4
setup geant4 v4_11_2_p02 -q "e26:prof"

# Below is copied from NDCAFMaker setup script
setup gcc v9_3_0
setup pycurl
setup ifdhc
setup genie_xsec   v3_04_00 -q AR2320i00000:e1000:k250
setup genie_phyopt v3_04_00 -q dkcharmtau
setup jobsub_client
setup eigen v3_3_5
setup srproxy v00.44 -q py3913 #Added 
setup hdf5 v1_10_5a -q e20
setup fhiclcpp v4_15_03 -q debug:e20

# edep-sim needs to know where a certain GEANT .cmake file is...
G4_cmake_file=$(find "${GEANT4_FQ_DIR}" -type f -name 'Geant4Config.cmake' -print -quit)
if [ -z "$G4_cmake_file" ]; then
  echo "ERROR: Geant4Config.cmake not found under ${GEANT4_FQ_DIR}" >&2
else
  export Geant4_DIR=$(dirname "$G4_cmake_file")
fi

# edep-sim needs to have the GEANT bin directory in the path
export PATH=$PATH:$GEANT4_FQ_DIR/bin

# shut up ROOT include errors
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$GENIE_INC/GENIE

# add SRProxy path to ROOT to allow for use in ROOT macros
# can use duneanaobj SRProxy
#export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$SRPROXY_INC

########### TODO: build own version of nusystematics and point to it here
# nusystematics paths
#export NUSYST=/global/u2/r/rdiurba/nusystematics/nusystematics
#export LD_LIBRARY_PATH=${NUSYST}/build/Linux/lib:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=${NUSYST}/build/nusystematics/artless:$LD_LIBRARY_PATH
#export FHICL_FILE_PATH=${NUSYST}/nusystematics/fcl:$FHICL_FILE_PATH

# Add pyGeoEff to pythonpath
export PYTHONPATH=${PYTHONPATH}:${PWD}/DUNE_ND_GeoEff/lib/

# duneananobj needs to be in the libs too
export LD_LIBRARY_PATH=${DUNEANAOBJ_LIB}:$LD_LIBRARY_PATH

# Add yaml-cpp to the paths
export LD_LIBRARY_PATH=${PWD}/../spine_ana/yaml-cpp/install/lib64:$LD_LIBRARY_PATH

# finally, add our lib & bin to the paths
mydir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export LD_LIBRARY_PATH=$mydir/lib:$LD_LIBRARY_PATH
export PATH=$mydir/bin:$PATH

# our FCL needs to be findable too -- we don't have a FHICL
#export FHICL_FILE_PATH="$FHICL_FILE_PATH:$mydir/cfg"

########### TODO: build own version of ROOUNFOLD and point to it here
#export ROOUNFOLD=/global/homes/r/rdiurba/cafAna/RooUnfold
#export LD_LIBRARY_PATH=$ROOUNFOLD:$LD_LIBRARY_PATH
#export C_INCLUDE_PATH="$ROOUNFOLD:$C_INCLUDEPATH"
#export CPLUS_INCLUDE_PATH=$ROOUNFOLD:$CPLUS_INCLUDE_PATH
#export PATH="$ROOUNFOLD:$PATH"

#export LD_LIBRARY_PATH=$ROOUNFOLD/src:$LD_LIBRARY_PATH
#export C_INCLUDE_PATH="$ROOUNFOLD/src:$C_INCLUDEPATH"
#export CPLUS_INCLUDE_PATH=$ROOUNFOLD/src:$CPLUS_INCLUDE_PATH
#export PATH="$ROOUNFOLD/src:$PATH"

########### TODO: build own version of nusystematics and point to it here
#export mywd=/global/u2/r/rdiurba/nusystematics

setup boost v1_80_0 -q e26:prof
setup tbb v2021_7_0 -q e26
setup sqlite v3_39_02_00

#export fhiclcpp_ROOT=${mywd}/fhicl-cpp-standalone/build/fhicl-cpp-install/
#export PATH=${fhiclcpp_ROOT}/bin/:${PATH}
#export LD_LIBRARY_PATH=${fhiclcpp_ROOT}/lib/:${LD_LIBRARY_PATH}
#
#export cetlib_ROOT=${mywd}/fhicl-cpp-standalone/build/cetlib-install/
#export PATH=${cetlib_ROOT}/bin/:${PATH}
#export LD_LIBRARY_PATH=${cetlib_ROOT}/lib/:${LD_LIBRARY_PATH}
#
#export cetlib_except_ROOT=${mywd}/fhicl-cpp-standalone/build/cetlib-except-install/
#export PATH=${cetlib_except_ROOT}/bin/:${PATH}
#export LD_LIBRARY_PATH=${cetlib_except_ROOT}/lib/:${LD_LIBRARY_PATH}
#
#export hep_concurrency_ROOT=${mywd}/fhicl-cpp-standalone/build/hep-concurrency-install
#export PATH=${hep_concurrency_ROOT}/bin/:${PATH}
#export LD_LIBRARY_PATH=${hep_concurrency_ROOT}/lib/:${LD_LIBRARY_PATH}

########### TODO: build own version of nusystematics and point to it here
#source ${NUSYST}/build/Linux/bin/setup.nusystematics.sh
#export PATH=/global/u2/r/rdiurba/nusystematics/ndnusyst/ndnusyst-src/build/Linux/bin/:${PATH}

#source ../../2x2_sim/util/init.inc.sh -- not necessary for analysis, just for running 2x2_sim

