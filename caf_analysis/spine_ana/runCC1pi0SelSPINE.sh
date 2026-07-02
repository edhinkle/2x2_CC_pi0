#!/usr/bin/env bash
# Run CC1pi0 analysis and plotting scripts.
# Partially based on https://github.com/rdiurba/2x2_trackMultStudies/blob/master/spineAna/runAnalysis.sh

####################################################################################
################################ PARAMETERS + SETUP ################################
####################################################################################

# Run analysis, plotting, or both? 
RUN_ANALYSIS=1
RUN_PLOTTING=1

# Set parameters for running analysis
RUN_NAME="testCC1pi0SelSPINE_MATCHOLD"
INPUTFILELIST="$(pwd)/include/config_data/MiniRun65FileList.txt"
OUTPUTFILE="${RUN_NAME}.root"
MCONLYSTRING="1"
CONFIGFILEPATH="$(pwd)/include/config_data/CC1pi0.yaml"
OUTPUTDIR="$(pwd)/outputs/${RUN_NAME}/"

# Name for saving configuration file used for this run (for reproducibility)
CONFIGFILE_SAVE_NAME="${RUN_NAME}_config.yaml"

# Additional plotting-only parameters
OUTPUTFILEPATH="${OUTPUTDIR}${OUTPUTFILE}"
THEME="default"

# Pre-container invocation setup
source ../setup/pre_container_setup.sh
setup_cuda

# Make output directory if it doesn't exist
if [[ ! -d "${OUTPUTDIR}" ]]; then
    mkdir -p "${OUTPUTDIR}"
fi

####################################################################################
##################################### ANALYSIS #####################################
####################################################################################
if [[ $RUN_ANALYSIS -eq 1 ]]; then
    # Invoke Fermilab SL7 container to compile and run analysis code on CAFs (old container name: fermilab/fnal-wn-sl7:latest)
    shifter --image=fermilab/fnal-dev-sl7:latest --module=cvmfs,gpu /bin/bash << EOF1

        # Set up environment for analysis code
        source /environment
        source ../setup/caf_env_setup.sh

        # Compile analysis code
        echo "Compiling CC1pi0 selection code..."
        source compileCC1pi0SelSPINE.sh 

        # Run the analysis
        echo "Running CC1pi0 analysis code..."
        ./spine_CC1pi0_sel ${INPUTFILELIST} ${OUTPUTFILE} ${MCONLYSTRING} ${CONFIGFILEPATH} 
# EOF must be at the beginning of the line with no spaces before it or after it
EOF1

    #./dlp_fluxRw_sel ${INPUTFILELIST} testFluxRWSPINE.root 1
    #./dlp_genieRw_sel ${INPUTFILELIST} testGENIERWSPINE.root 1

    #./spine_CC1pi0_sel ${INPUTFILELIST} ${OUTPUTFILE} ${MCONLYSTRING} ${CONFIGFILEPATH}
    #./fakedata_dlp_sel ${INPUTFILELIST} testFakeDataHighEReweight.root 1
    #./fakedata_dlp_selLowE ${INPUTFILELIST} testFakeDataLowEReweight.root 1

    # Delete output file in output directory if it already exists (e.g. from a previous run)
    if [[ -f "${OUTPUTDIR}/${OUTPUTFILE}" ]]; then
        rm -f "${OUTPUTDIR}/${OUTPUTFILE}"
    fi

    # Move new output file to output directory
    mv ${OUTPUTFILE} ${OUTPUTDIR} 

    # Delete config file in output directory if it already exists (e.g. from a previous run)
    if [[ -f "${OUTPUTDIR}/${CONFIGFILE_SAVE_NAME}" ]]; then
        rm -f "${OUTPUTDIR}/${CONFIGFILE_SAVE_NAME}"
    fi

    # Copy new config file to output directory
    cp ${CONFIGFILEPATH} ${OUTPUTDIR}/${CONFIGFILE_SAVE_NAME}

fi

####################################################################################
##################################### PLOTTING #####################################
####################################################################################
if [[ $RUN_PLOTTING -eq 1 ]]; then

    # Invoke DLP container to run plotting scripts on output analysis file
    # NOTE: DLP container is used because it has the necessary python dependencies for plotting,
    #       and I've had issues with python inside the SL7 container conflicting with NERSC
    #       python. Ideally, we only use one container and/or don't need a container for plotting.
    shifter --image=deeplearnphysics/larcv2:ub2204-cu121-torch251-larndsim bash << EOF2

        # Install python plotting code in CC1pi0 repo
        echo "Installing Python plotting dependencies..."
        pip install -e .

        # Run plotting script
        echo "Plotting histograms from output file ${OUTPUTFILEPATH} with theme ${THEME}..."
        python plotanapy/scripts/plot_histograms.py --root_file ${OUTPUTFILEPATH} --config_file ${CONFIGFILEPATH} --theme ${THEME} --all_plots_dir ${OUTPUTDIR} --mcOnly ${MCONLYSTRING}

# EOF must be at the beginning of the line with no spaces before it or after it
EOF2

fi


#root -l -b -q unfoldGENIERW.C
#root -l -b -q unfoldFluxRW.C 
#root -l -b -q testEff.C
#root -l -b -q plotSPINE.C