// Macro to read the CAF files and select muon neutrino interactions with 1 pi0 in final state
// and save histograms for unfolding and plots.

// Author: Elise Hinkle (ehinkle@uchicago.edu)
// Adapted from: https://github.com/rdiurba/2x2_trackMultStudies/blob/master/spineAna/dlp_sel.cc

#include "TChain.h"
#include "TEfficiency.h"
#include "TFile.h"
// #include "../fluxSyst/FluxCovarianceReweight.h" // TODO: Add flux systematics

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "duneanaobj/StandardRecord/StandardRecord.h" //Ideally, this should be SRProxy.h, but there is an include error for that now. Alternatively, you can use SetBranchStatus function in TreeLoader, but it does not work for the common branch (to do)
#include <fstream>
#include <iostream>
#include <string>

// Local includes:
#include "config/ConfigLoader.h"
#include "io/CAFUtils.h"
#include "plotting/HistogramManager.h"
#include "selection/TruthSelection.h"
#include "selection/RecoSelection.h"
#include "systematics/flux/FluxSystManager.h"


// Main macro method
int main(int argc, char **argv) {

    // -------------------------------
    // 1. Check for correct number of arguments
    // -------------------------------
    if (argc != 5) {
        std::cout << "\n USAGE: " << argv[0]
            << "input_caf_file_list output_root_file mcOnlyString configFilepath\n"
            << std::endl;
        return 1;
    }

    // -------------------------------
    // 2. Get arguments and set mcOnly flag
    // -------------------------------
    std::string input_file_list = argv[1];
    std::string output_rootfile = argv[2];
    std::string mcOnlyString = argv[3];
    std::string configFilepath = argv[4];
    bool mcOnly = true;
    if (mcOnlyString == "0")
        mcOnly = false;
    std::cout << "MC Only: " << mcOnly << "," << argv[3] << std::endl;

    // -------------------------------
    // 3. Load configurations
    // -------------------------------
    SelectionConfig selectionCuts = config::LoadSelectionConfig(configFilepath);
    BeamConfig beamInfo = config::LoadBeamConfig(configFilepath);
    DetectorConfig detInfo = config::LoadDetectorConfig(configFilepath);

    FluxSystConfig fluxSystInfo = config::LoadFluxSystConfig(configFilepath);

    // -------------------------------
    // 4. Initialize systematics
    // -------------------------------   
    FluxSystManager fluxSystManager(fluxSystInfo);

    // -------------------------------
    // 5. Initialize histograms
    // -------------------------------
    HistogramManager hist(fluxSystInfo.nThrows, fluxSystManager.GetNominalIntegratedFlux()); 

    // -------------------------------
    // 6. Load CAF chain + set up branch address/total spill/POT variables
    // -------------------------------
    TChain* caf_chain = io::BuildCAFChain(input_file_list);
    
      // Check that CAF chain loaded correctly 
    if (!caf_chain) {
      std::cerr << "Failed to build CAF chain from file list: " << input_file_list << "\n" << std::endl;
      return 1;
    }

      // Define Standard Record and link it to the CAF tree branch "rec"
    auto sr = new caf::StandardRecord;
    caf_chain->SetBranchAddress("rec", &sr);

      // Set number of entries & fill histogram for total spills processed
    long Nentries = caf_chain->GetEntries();
    std::cout << "Total spills: " << Nentries << std::endl;
    hist.reco.FillTotalSpillsProcessed(Nentries);

      // Set totalPOT to 0 before looping over spills
    double totalPOT = 0;

    // -------------------------------
    // 7. Loop over spills
    // -------------------------------
    for (long i = 0; i < Nentries; ++i) {
      
      caf_chain->GetEntry(i);
      // Add the POT for this spill to the totalPOT
      totalPOT = sr->beam.pulsepot / 1e13 + totalPOT;

      // -------------------------------------------------------------
      // 7a. Truth Interaction selection (if MC)
      // -------------------------------------------------------------
      if (mcOnly) {
        TruthSelection truthSel(selectionCuts, beamInfo, detInfo, fluxSystManager);
        truthSel.SelectTruthInteractions(*sr, hist);
      }

      // -------------------------------------------------------------
      // 7b. Reco interaction selection (includes truth backtracking)
      // -------------------------------------------------------------
      RecoSelection recoSel(selectionCuts, beamInfo, detInfo, fluxSystManager, mcOnly);
      recoSel.SelectRecoInteractions(*sr, hist);
    }
    // -------------------------------
    // 8. Fill Total POT histogram after looping over spills
    // -------------------------------
    hist.reco.FillTotalPOT(totalPOT);
    // ----------------------------------------------------------------------------
    // 9. Finalize flux systematics histograms including profiles from throw hists
    // ----------------------------------------------------------------------------
    hist.fluxSyst.SetThrowFluxTitles(fluxSystManager.GetIntegratedFluxThrows());
    hist.fluxSyst.Finalize(); // builds profiles from completed throw histograms
    // -------------------------------
    // 10. Finalize + write output
    // -------------------------------
    TFile *caf_out_file = new TFile(output_rootfile.c_str(), "recreate");
    hist.Write(caf_out_file);
    caf_out_file->Close();

    return 0;
}