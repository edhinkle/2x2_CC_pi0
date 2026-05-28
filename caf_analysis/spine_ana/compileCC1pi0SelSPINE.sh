#!/usr/bin/env bash
# Compile the C++ code using g++. Partially based on https://github.com/rdiurba/2x2_trackMultStudies/blob/master/spineAna/compile.sh

# Compile code
g++ -std=c++17 -O2 -g -pedantic -Wall -I. -Iinclude `root-config --cflags --glibs` \
  src/analysis/spine_CC1pi0_sel.cc \
  src/config/ConfigLoader.cc \
  src/io/CAFUtils.cc \
  src/plotting/HistogramManager.cc \
  src/plotting/CutFlowManager.cc \
  src/plotting/RecoHists.cc \
  src/plotting/TruthHists.cc \
  src/plotting/TruthMatchedHists.cc \
  src/selection/RecoSelection.cc \
  src/selection/TruthSelection.cc \
  -lyaml-cpp \
  -o spine_CC1pi0_sel

#g++ -std=c++17 -O2 -g -pedantic -Wall `root-config --cflags --glibs` dlp_sel.cc ../fluxSyst/FluxCovarianceReweight.cc -o dlp_sel -I/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_06_01b/include -I$SRPROXY_INC -L$DUNEANAOBJ_LIB -lduneanaobj_StandardRecordProxy -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord_dict

#g++ -std=c++17 -O2 -g -pedantic -Wall `root-config --cflags --glibs` dlp_genieRw_sel.cc ../fluxSyst/FluxCovarianceReweight.cc  -o dlp_genieRw_sel -I/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_06_01b/include -I$SRPROXY_INC -L$DUNEANAOBJ_LIB -lduneanaobj_StandardRecordProxy -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord_dict

#g++ -std=c++17 -O2 -g -pedantic -Wall `root-config --cflags --glibs` dlp_fluxRw_sel.cc  ../fluxSyst/FluxCovarianceReweight.cc -o dlp_fluxRw_sel -I/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_06_01b/include -I$SRPROXY_INC -L$DUNEANAOBJ_LIB -lduneanaobj_StandardRecordProxy -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord_dict


if [ $? -eq 0 ]; then
    echo "Compilation successful"
else
    echo "Compilation failed"
fi