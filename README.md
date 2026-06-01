# 2x2_CC_pi0

This repository is used to run analyses related to CC $\nu_\mu$ events in the 2x2 Demonstrator which include at least one $\pi^0$ in the final state. The `truth_studies` directory includes code used to study such interactions using truth-level information from 2x2 simulation files. The `reco_studies` directory includes code used to study shower reconstruction in SPINE using 2x2 simulation files. The `caf_studies` directory includes code used to perform initial analyses running over 2x2 CAF files. The `caf_analysis` directory includes the main code used to run the final analysis. The latter two directories include code which use some scripts from https://github.com/rdiurba/2x2_trackMultStudies/tree/master/spineAna as reference points.

## CAF Analysis

The analysis directory is split into `setup` and `spine_ana` (analysis) subdirectories. Scripts in `setup` must be run before running the analysis scripts. Prior to this, two other repositories must be installed inside the `2x2_CC_pi0` repository. First, the `2x2_sim` repository (https://github.com/DUNE/2x2_sim) must be cloned inside the main folder of `2x2_CC_pi0` (one level above `caf_analysis`). Second, `yaml-cpp` (https://github.com/jbeder/yaml-cpp.git) must be cloned inside the `caf_analysis/spine_ana` directory and built: 

    cd yaml-cpp
    mkdir build && cd build
    cmake -DCMAKE_INSTALL_PREFIX=../install ..
    make
    make install

### `setup`

The three scripts in this directory are used to set up the environment necessary to run the CAF Analysis. `setup_container.sh` sets up an SL7 container and is run first, followed by `setup_caf.sh`. Currently, `setup_general.sh` is extraneous and included for reference as the repository grows. 

### `spine_ana`

This is the main analysis directory. The analysis is run in the direcotory using `runCC1pi0SelSPINE.sh`. This script first compiles all necessary code using `compileCC1pi0SelSPINE.sh` and then runs the main analysis script, `src/analysis/spine_CC1pi0_sel.cc` with inputs configurable in `runCC1pi0SelSPINE.sh`. Analysis code is written in C++ and organized in the `include` (header files) and `src` (source code) subdirectories. The scripts are further organized as follows:

#### `config` and `config_data`

#### `io`

#### `plotting`

#### `selection`

#### `cuts`

#### `analysis`
