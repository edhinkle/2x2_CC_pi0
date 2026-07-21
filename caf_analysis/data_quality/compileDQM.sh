#!/usr/bin/env bash
# Compile the C++ code using g++. Partially based on https://github.com/rdiurba/2x2_trackMultStudies/blob/master/spineAna/compile.sh

# Compile code
g++ -std=c++17 -O2 -g -gdwarf-4 -pedantic -Wall -I. -Iinclude \
  -I/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_06_01b/include \
  -I./../spine_ana/yaml-cpp/install/include \
  `root-config --cflags --glibs` \
  src/analysis/physicsDQM.cc \
  src/io/CAFUtils.cc \
  src/plotting/HistogramManager.cc \
  src/plotting/BeamHists.cc \
  -L$DUNEANAOBJ_LIB \
  -lduneanaobj_StandardRecord \
  -L./../spine_ana/yaml-cpp/install/lib64 \
  -lyaml-cpp \
  -o physicsDQM

#g++ -std=c++17 -O2 -g -pedantic -Wall `root-config --cflags --glibs` dlp_sel.cc ../fluxSyst/FluxCovarianceReweight.cc -o dlp_sel -I/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_06_01b/include -I$SRPROXY_INC -L$DUNEANAOBJ_LIB -lduneanaobj_StandardRecordProxy -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord_dict

#g++ -std=c++17 -O2 -g -pedantic -Wall `root-config --cflags --glibs` dlp_genieRw_sel.cc ../fluxSyst/FluxCovarianceReweight.cc  -o dlp_genieRw_sel -I/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_06_01b/include -I$SRPROXY_INC -L$DUNEANAOBJ_LIB -lduneanaobj_StandardRecordProxy -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord_dict

#g++ -std=c++17 -O2 -g -pedantic -Wall `root-config --cflags --glibs` dlp_fluxRw_sel.cc  ../fluxSyst/FluxCovarianceReweight.cc -o dlp_fluxRw_sel -I/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_06_01b/include -I$SRPROXY_INC -L$DUNEANAOBJ_LIB -lduneanaobj_StandardRecordProxy -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord_dict


if [ $? -eq 0 ]; then
    echo "Compilation successful"
else
    echo "Compilation failed"
fi