source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup duneanaobj v03_06_01b -q e20:prof
setup dk2nugenie   v01_10_01k -q debug:e20
setup cmake v3_27_4
setup geant4 v4_11_2_p02 -q "e26:prof"

# edep-sim needs to have the GEANT bin directory in the path
export PATH=$PATH:$GEANT4_FQ_DIR/bin

# shut up ROOT include errors
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$GENIE_INC/GENIE

# nusystematics paths
export NUSYST=/global/u2/r/rdiurba/nusystematics/nusystematics
export LD_LIBRARY_PATH=${NUSYST}/build/Linux/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${NUSYST}/build/nusystematics/artless:$LD_LIBRARY_PATH
export FHICL_FILE_PATH=${NUSYST}/nusystematics/fcl:$FHICL_FILE_PATH

# Add pyGeoEff to pythonpath
export PYTHONPATH=${PYTHONPATH}:${PWD}/DUNE_ND_GeoEff/lib/

# duneananobj needs to be in the libs too
export LD_LIBRARY_PATH=${DUNEANAOBJ_LIB}:$LD_LIBRARY_PATH

# finally, add our lib & bin to the paths
mydir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export LD_LIBRARY_PATH=$mydir/lib:$LD_LIBRARY_PATH
export PATH=$mydir/bin:$PATH


export ROOUNFOLD=/global/homes/r/rdiurba/cafAna/RooUnfold
export LD_LIBRARY_PATH=$ROOUNFOLD:$LD_LIBRARY_PATH
export C_INCLUDE_PATH="$ROOUNFOLD:$C_INCLUDEPATH"
export CPLUS_INCLUDE_PATH=$ROOUNFOLD:$CPLUS_INCLUDE_PATH
export PATH="$ROOUNFOLD:$PATH"

export LD_LIBRARY_PATH=$ROOUNFOLD/src:$LD_LIBRARY_PATH
export C_INCLUDE_PATH="$ROOUNFOLD/src:$C_INCLUDEPATH"
export CPLUS_INCLUDE_PATH=$ROOUNFOLD/src:$CPLUS_INCLUDE_PATH
export PATH="$ROOUNFOLD/src:$PATH"

export mywd=/global/u2/r/rdiurba/nusystematics

setup boost v1_80_0 -q e26:prof
setup tbb v2021_7_0 -q e26
setup sqlite v3_39_02_00

export fhiclcpp_ROOT=${mywd}/fhicl-cpp-standalone/build/fhicl-cpp-install/
export PATH=${fhiclcpp_ROOT}/bin/:${PATH}
export LD_LIBRARY_PATH=${fhiclcpp_ROOT}/lib/:${LD_LIBRARY_PATH}

export cetlib_ROOT=${mywd}/fhicl-cpp-standalone/build/cetlib-install/
export PATH=${cetlib_ROOT}/bin/:${PATH}
export LD_LIBRARY_PATH=${cetlib_ROOT}/lib/:${LD_LIBRARY_PATH}

export cetlib_except_ROOT=${mywd}/fhicl-cpp-standalone/build/cetlib-except-install/
export PATH=${cetlib_except_ROOT}/bin/:${PATH}
export LD_LIBRARY_PATH=${cetlib_except_ROOT}/lib/:${LD_LIBRARY_PATH}

export hep_concurrency_ROOT=${mywd}/fhicl-cpp-standalone/build/hep-concurrency-install
export PATH=${hep_concurrency_ROOT}/bin/:${PATH}
export LD_LIBRARY_PATH=${hep_concurrency_ROOT}/lib/:${LD_LIBRARY_PATH}


source ${NUSYST}/build/Linux/bin/setup.nusystematics.sh
export PATH=/global/u2/r/rdiurba/nusystematics/ndnusyst/ndnusyst-src/build/Linux/bin/:${PATH}