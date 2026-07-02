# 2x2_CC_pi0

This repository is used to run analyses related to CC $\nu_\mu$ events in the 2x2 Demonstrator which include at least one $\pi^0$ in the final state. The `truth_studies` directory includes code used to study such interactions using truth-level information from 2x2 simulation files. The `reco_studies` directory includes code used to study shower reconstruction in SPINE using 2x2 simulation files. The `caf_studies` directory includes code used to perform initial analyses running over 2x2 CAF files. The `caf_analysis` directory includes the main code used to run the final analysis. The latter two directories include code which use some scripts from https://github.com/rdiurba/2x2_trackMultStudies/tree/master/spineAna as reference points. This repository also includes code generated partially by LLMs including models managed by ChatGPT, Claude, and GitHub Copilot.

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

This is the main analysis directory. The analysis is run in the direcotory using `runCC1pi0SelSPINE.sh`. This script first compiles all necessary code using `compileCC1pi0SelSPINE.sh` and then runs the main analysis script, `src/analysis/spine_CC1pi0_sel.cc` with inputs configurable in `runCC1pi0SelSPINE.sh`. At the top of this script, the `RUN_ANALYSIS` and `RUN_PLOTTING` variables are set (1 = run this portion of the code, 0 = don't run). Other varibles pertaining to configuration scripts and outfile names are also set in `runCC1pi0SelSPINE.sh`. Plotting requires running the analysis at least once. There is also a `RUN_NAME` variable which should be changed when testing different configurations. The output files are saved in a subdirectory of the `outputs` directory named after the `RUN_NAME`. The version of the configuration YAML (see below) used for that run is also saved there. Analysis code is written in C++ and organized in the `include` (header files) and `src` (source code) subdirectories. Plotting code is written in Python and organized in `plotanapy`. `setup.py` allows for all plotting scripts to exist in a pip-installable Python package. The scripts are further organized as follows.

#### ANALYSIS: `config_data` and `config`

All configurable detector boundaries, beam information, and cuts are recorded in `include/config_data/CC1pi0.yaml`. The `include/config_data/` directory also includes the text file list of files over which to run the analysis, e.g. `MiniRun65FileList.txt` includes all MiniRun6.5 files, with each line of the file listing the location of a single file on the NERSC Perlmutter system. The flux covariance binning file and the NuMI flux systematics file, `flux_covariance_binning_NuMI_GeV.txt` and `numiFluxSystFull.root`, are located in `include/config_data/` as well. There is no `src/config_data`.

The configurable quantities are split into selection, beam, detector, and flux systematics categories. Each category has an associated struct object, defined in `include/config/SelectionConfig.h`, `include/config/BeamConfig.h`, `include/config/DetectorConfig.h`, and `include/config/FluxSystConfig.h` respectively. The information from `include/config_data/CC1pi0.yaml` is loaded into the structs as directed by the loading methods in the `ConfigLoader` namespace, defined by `include/config/ConfigLoader.h` and `src/config/ConfigLoader.cc`. Once loaded into the structs, selection, detector, and beam information are accessible by the analysis code without hard-coding of quantities such that changing the value of a cut only requires making changes in the `include/config_data/CC1pi0.yaml` file (or creating a new configuration file and changing the config file path in `runCC1pi0SelSPINE.sh`). 

#### ANALYSIS: `io`

CAF Chain building is handled in `src/io/CAFUtils`. This module allows for easy reading of a full sample of events as denoted by the `include/config_data` file list text file, e.g. `MiniRun65FileList.txt`. 

#### ANALYSIS: `plotting`

The `plotting` directory in the C++ analysis modules is using to handle histogram declaration, management, filling, and writing to the output file. The `HistogramManager` class deals with all histogram categories and owns instances of each, as well as a global histogram writing method which separates histgorams in different categories into different directories in the final output file (`truth`, `reco`, `truthMatched`, and `cutflows`). `RecoHists` handles histograms filled with purely reconstructed quantities (independent of a truth-level selection or any truth backtracking for analysis of simulation files), `TruthHists` handles histograms filled with purely truth-level quantities (independent of any selection using reconstructed quantities), and `TruthMatchedHists` handles histograms containing or informed by a mix of true and reconstructed information (e.g. reco quantities for reco events matched to a true signal event, true vs. reco response matrices, etc.). When each class is instantiated, all histograms corresponding to that class are declared. Filling of each histogram occurs using corresponding class methods which are called in the `selection` scripts (see next section). Finally, each class includes a `Write()` method which is called as part of the global `HistogramManager::Write()` method. The final class in the `plotting` subcategory, `CutFlowManager`, handles methods related to creating histograms which store how many events are kept after each selection cut. 

#### ANALYSIS: `selection`

The `selection` analysis code contains the two main classes controlling event selection, `TruthSelection` and `RecoSelection`. Both selections use information from the beam, selection, and detector configurations. The main method in `TruthSelection`, `TruthSelection::SelectTruthInteractions`, runs through truth-level interactions in a single spill, selects for true $\nu_{\mu}$ and $\bar{\nu}\_{\mu}$ -argon CC $1\pi^{0}$ events, and fills the relevant `TruthHists` and cutflow histograms. `TruthSelection::BuildTruthSummary`, the other main method in this class, saves pertinent information about a truth-level event in a `TruthInteractionSummary` struct defined in `include/src/analysis/TruthInteractionSummary.h`. The `RecoSelection` class is similar to `TruthSelection`, but for reconstructed interactions, therefore also requiring information as to whether the input file is simulation or data. The main method, `RecoSelection::SelectRecoInteractions` runs through all reco interactions in a single spill and implements cuts designed to select for $\nu_{\mu}$ and $\bar{\nu}_{\mu}$ -argon CC $1\pi^{0}$ events from reco variables. `RecoHists` and cutflow histograms are filled as part of this method, as are `TruthMatchedHists` corresponding to the best-matched truth interactions or particles for each reco interaction (if the input file is simulation and truth-level information exists). The other main class methods include `RecoSelection::BuildRecoSummary` and `RecoSelection::BuildMatchedIxnSummary`, which save information about a single reco interaction and a single truth-matched interaction in `RecoInteractionSummary` and `MatchedInteractionSummary` structs defined in `include/src/analysis/RecoInteractionSummary.h` and `include/src/analysis/MatchedInteractionSummary.h`, respectively.

#### ANALYSIS: `cuts`

The classes in the `cuts` subdivision define more complicated cuts for the truth- and reco-level selections. `DetectorCuts` uses settings from the detector configuration struct to set up cuts on vertex or point location (`DetectorCuts::InModuleVolumes` and `DetectorCuts::InFiducialVolume`) and whether or not true interactions/particles are expected to be over a kinetic energy threshold defined in `include/config_data/CC1pi0.yaml` (`DetectorCuts::TrueIxnAboveKEThreshold` and `DetectorCuts::TrueParticleAboveKEThreshold`). `Mx2Matcher` contains methods necessary to match tracks in 2x2 to tracks in downstream Mx2 in order to select for tracks starting in 2x2 and punching through downstream Mx2 which are hypothesized to be muons from $\nu_\mu$/$\bar{\nu}_\mu$-argon CC interactions. The main class method, `Mx2Matcher::MatchInteraction`, fills a `Mx2MatchResult` struct defined in `include/cuts/Mx2MatchResult.h` which holds information about the successfully matched Mx2 and 2x2 tracks, as well as relevant truth-backtracking information if the analysis is being run on simulation. 

#### ANALYSIS: `analysis`

The only code in the `src/analysis` directory is the main analysis script, `spine_CC1pi0_sel.cc`, which is explained above. `include/analysis` holds the header files defining the structs `TruthInteractionSummary`, `RecoInteractionSummary`, and `MatchedInteractionSummary` which contain relevant information about individual truth interactions, reco interactions, and best-matched truth backtracked interactions while the selections are run.

#### ANALYSIS: `systematics/flux`

Currently, only flux systematics are incorporated into the selection workflow. The `FluxCovarianceReweight` class contains methods for generating throws for a collection of universes with slightly varied flux settings. TODO: Finish this section

#### PLOTTING: `plotanapy/plotting`

The main plotting directory includes `histogram_loader.py` as well as three subdirectories, `plotters`, `styles`, and `utils`. `histogram_loader.py` defines the `HistogramLoader` class, which includes methods helpful for loading histograms from the output analysis ROOT file. 

#### PLOTTING: `plotanapy/plotting/plotters`

This subdirectory includes scripts providing class definitions for plotting histograms of each category defined during the analysis, `truth`, `reco`, `truthMatched`, `

#### PLOTTING: `plotanapy/plotting/styles`

#### PLOTTING: `plotanapy/plotting/utils`

#### PLOTTING: `plotanapy/scripts`
