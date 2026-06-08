# 2x2_CC_pi0

This repository is used to run analyses related to CC $\nu_\mu$ events in the 2x2 Demonstrator which include at least one $\pi^0$ in the final state. The `truth_studies` directory includes code used to study such interactions using truth-level information from 2x2 simulation files. The `reco_studies` directory includes code used to study shower reconstruction in SPINE using 2x2 simulation files. The `caf_studies` directory includes code used to perform initial analyses running over 2x2 CAF files. The `caf_analysis` directory includes the main code used to run the final analysis. The latter two directories include code which use some scripts from https://github.com/rdiurba/2x2_trackMultStudies/tree/master/spineAna as reference points.

## CAF Analysis

The analysis directory is split into `setup` and `spine_ana` (analysis) subdirectories. This analysis is designed to be run on the NERSC Perlmutter system. Prior to running the analysis, `yaml-cpp` (https://github.com/jbeder/yaml-cpp.git) must be cloned inside the `caf_analysis/spine_ana` directory and built: 

    cd yaml-cpp
    mkdir build && cd build
    cmake -DCMAKE_INSTALL_PREFIX=../install ..
    make
    make install

### `setup`

The two scripts in this directory are used to set up the environment necessary to run the CAF Analysis. `pre_container_setup.sh` defines how Python and CUDA are set up on NERSC. `caf_env_setup.sh` loads environment variables for running the analysis in a Fermilab SL7 container. Both scripts are run as part of the main analysis script (`runCC1pi0SelSPINE.sh`) and don't need to be run separately by the analyzer. 

### `spine_ana`

This is the main analysis directory. The analysis is run in the direcotory using `runCC1pi0SelSPINE.sh`. This script first compiles all necessary code using `compileCC1pi0SelSPINE.sh` and then runs the main analysis script, `src/analysis/spine_CC1pi0_sel.cc` with inputs configurable in `runCC1pi0SelSPINE.sh`. Analysis code is written in C++ and organized in the `include` (header files) and `src` (source code) subdirectories. The scripts are further organized as follows:

#### `config_data` and `config`

All configurable detector boundaries, beam information, and cuts are recorded in `include/config_data/CC1pi0.yaml`. The `include/config_data/` directory also includes the text file list of files over which to run the analysis, e.g. `MiniRun65FileList.txt` includes all MiniRun6.5 files, with each line of the file listing the location of a single file on the NERSC Perlmutter system. There is no `src/config_data`.

The configurable quantities are split into selection, beam, and detector categories. Each category has an associated struct object, defined in `include/config/SelectionConfig.h`, `include/config/BeamConfig.h`, and `include/config/DetectorConfig.h`, respectively. The information from `include/config_data/CC1pi0.yaml` is loaded into the structs as directed by the loading methods in the `ConfigLoader` namespace, defined by `include/config/ConfigLoader.h` and `src/config/ConfigLoader.cc`. Once loaded into the structs, selection, detector, and beam information are accessible by the analysis code without hard-coding of quantities such that changing the value of a cut only requires making changes in the `include/config_data/CC1pi0.yaml` file (or creating a new configuration file and changing the config file path in `runCC1pi0SelSPINE.sh`). 

#### `io`

#### `plotting`

#### `selection`

#### `cuts`

#### `analysis`
