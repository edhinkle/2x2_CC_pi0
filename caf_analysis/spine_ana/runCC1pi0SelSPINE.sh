#!/usr/bin/env bash
# Run CC1pi0 analysis. 
# Partially based on https://github.com/rdiurba/2x2_trackMultStudies/blob/master/spineAna/runAnalysis.sh

# Set parameters for running analysis
INPUTFILELIST="$(pwd)/include/config_data/MiniRun65FileList.txt"
OUTPUTFILE="testCC1pi0SelSPINE.root"
MCONLYSTRING="1"
CONFIGFILEPATH="$(pwd)/include/config_data/CC1pi0.yaml"
OUTPUTDIR="$(pwd)/outputs/"

# Compile analysis code
echo "Compiling CC1pi0 selection code..."
source compileCC1pi0SelSPINE.sh

# Run the analysis
./spine_CC1pi0_sel ${INPUTFILELIST} ${OUTPUTFILE} ${MCONLYSTRING} ${CONFIGFILEPATH}
#./dlp_fluxRw_sel ${INPUTFILELIST} testFluxRWSPINE.root 1
#./dlp_genieRw_sel ${INPUTFILELIST} testGENIERWSPINE.root 1

#./spine_CC1pi0_sel ${INPUTFILELIST} ${OUTPUTFILE} ${MCONLYSTRING} ${CONFIGFILEPATH}
#./fakedata_dlp_sel ${INPUTFILELIST} testFakeDataHighEReweight.root 1
#./fakedata_dlp_selLowE ${INPUTFILELIST} testFakeDataLowEReweight.root 1

mv ${OUTPUTFILE} ${OUTPUTDIR}
#root -l -b -q unfoldGENIERW.C
#root -l -b -q unfoldFluxRW.C 
#root -l -b -q testEff.C
#root -l -b -q plotSPINE.C