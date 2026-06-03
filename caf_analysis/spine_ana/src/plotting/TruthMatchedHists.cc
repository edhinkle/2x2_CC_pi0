#include "plotting/TruthMatchedHists.h"
#include "TDirectory.h"

// Define all truth match histograms 
TruthMatchedHists::TruthMatchedHists()
{

    // Edges for CosL histograms
    Double_t edges[7] = {0.90, 0.96, 0.98, 0.9887, 0.994, 0.9974, 1};
    // Truth match mx2Track histograms
    // PDG
    h_truthMatchPrim_mx2TrackPDG = new TH1D("h_truthMatchPrim_mx2TrackPDG", "h_truthMatchPrim_mx2TrackPDG", 6000, -3000, 3000);
    h_truthMatchSec_mx2TrackPDG = new TH1D("h_truthMatchSec_mx2TrackPDG", "h_truthMatchSec_mx2TrackPDG", 6000, -3000, 3000);
    
    // Energy
    h_truthMatchPrim_mx2TrackE = new TH1D("h_truthMatchPrim_mx2TrackE", "h_truthMatchPrim_mx2TrackE", 50, 0, 20);
    h_truthMatchSec_mx2TrackE = new TH1D("h_truthMatchSec_mx2TrackE", "h_truthMatchSec_mx2TrackE", 50, 0, 20);

    // Cos Angle w/ beam
    h_truthMatchPrim_mx2TrackCosL = new TH1D("h_truthMatchPrim_mx2TrackCosL", "h_truthMatchPrim_mx2TrackCosL", 100, -1, 1);
    h_truthMatchSec_mx2TrackCosL = new TH1D("h_truthMatchSec_mx2TrackCosL", "h_truthMatchSec_mx2TrackCosL", 100, -1, 1);
    h_truthMatchIxn_Reco_responseMx2MatchPrimCosL =
      new TH2D("h_truthMatchIxn_Reco_responseMx2MatchPrimCosL", "h_truthMatchIxn_Reco_responseMx2MatchPrimCosL", 6, edges, 6, edges);
    h_truthMatchIxn_Reco_responseMx2MatchSecCosL =
      new TH2D("h_truthMatchIxn_Reco_responseMx2MatchSecCosL", "h_truthMatchIxn_Reco_responseMx2MatchSecCosL", 6, edges, 6, edges);

    // Muon kinematics for best matched interaction
    h_truthMatchIxn_MuonCosL = new TH1D("h_truthMatchIxn_MuonCosL", "h_truthMatchIxn_MuonCosL", 6, edges);
    h_truthMatchIxn_MuonCosL_unbinned = new TH1D("h_truthMatchIxn_MuonCosL_unbinned", "h_truthMatchIxn_MuonCosL_unbinned", 100, -1, 1);
    h_truthMatchIxn_MuonCosL_zoomOut = new TH1D("h_truthMatchIxn_MuonCosL_zoomOut", "h_truthMatchIxn_MuonCosL_zoomOut", 50, 0.8, 1);
    h_truthMatchIxn_MuonE = new TH1D("h_truthMatchIxn_MuonE", "h_truthMatchIxn_MuonE", 50, 0, 20);
    h_truthMatchIxn_Reco_responseMuonCosL =
      new TH2D("h_truthMatchIxn_Reco_responseMuonCosL", "h_truthMatchIxn_Reco_responseMuonCosL", 6, edges, 6, edges);

    // True - Reco Diff in Start Points for Mx2 Track Match in LAr
    h_truthMatchPrim_mx2TrackRecoDiffLArStartX =
      new TH1D("h_truthMatchPrim_mx2TrackRecoDiffLArStartX", "h_truthMatchPrim_mx2TrackRecoDiffLArStartX", 60, -15, 15);
    h_truthMatchSec_mx2TrackRecoDiffLArStartX =
      new TH1D("h_truthMatchSec_mx2TrackRecoDiffLArStartX", "h_truthMatchSec_mx2TrackRecoDiffLArStartX", 60, -15, 15);

    h_truthMatchPrim_mx2TrackRecoDiffLArStartY =
      new TH1D("h_truthMatchPrim_mx2TrackRecoDiffLArStartY", "h_truthMatchPrim_mx2TrackRecoDiffLArStartY", 60, -15, 15);
    h_truthMatchSec_mx2TrackRecoDiffLArStartY =
      new TH1D("h_truthMatchSec_mx2TrackRecoDiffLArStartY", "h_truthMatchSec_mx2TrackRecoDiffLArStartY", 60, -15, 15);

    h_truthMatchPrim_mx2TrackRecoDiffLArStartZ =
      new TH1D("h_truthMatchPrim_mx2TrackRecoDiffLArStartZ", "h_truthMatchPrim_mx2TrackRecoDiffLArStartZ", 60, -15, 15);
    h_truthMatchSec_mx2TrackRecoDiffLArStartZ =
      new TH1D("h_truthMatchSec_mx2TrackRecoDiffLArStartZ", "h_truthMatchSec_mx2TrackRecoDiffLArStartZ", 60, -15, 15);

    // Truth - Reco Diff in Vertex Position for Best Truth Match
    h_truthMatch_diffTruthRecoVertex = 
      new TH1D("h_truthMatch_diffTruthRecoVertex", "h_truthMatch_diffTruthRecoVertex", 30, 0, 15);

    // Truth match interaction pi0 multiplicity and shower multiplicity histograms
    // Bins and edges for multiplicity histograms
    //int binsMult = 10;
    //Double_t edgesMult[11] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    h_truthMatchIxn_PrimPi0Multiplicity = new TH1D("h_truthMatchIxn_PrimPi0Multiplicity", "h_truthMatchIxn_PrimPi0Multiplicity", 10, 0, 10);
    h_truthMatchIxn_SecPi0Multiplicity = new TH1D("h_truthMatchIxn_SecPi0Multiplicity", "h_truthMatchIxn_SecPi0Multiplicity", 10, 0, 10);
    h_truthMatchIxn_PrimShowerMultiplicity = new TH1D("h_truthMatchIxn_PrimShowerMultiplicity", "h_truthMatchIxn_PrimShowerMultiplicity", 20, 0, 20);
    h_truthMatchIxn_SecShowerMultiplicity = new TH1D("h_truthMatchIxn_SecShowerMultiplicity", "h_truthMatchIxn_SecShowerMultiplicity", 20, 0, 20);
    h_truthMatchIxn_PrimElectronMultiplicity = new TH1D("h_truthMatchIxn_PrimElectronMultiplicity", "h_truthMatchIxn_PrimElectronMultiplicity", 20, 0, 20);
    h_truthMatchIxn_SecElectronMultiplicity = new TH1D("h_truthMatchIxn_SecElectronMultiplicity", "h_truthMatchIxn_SecElectronMultiplicity", 20, 0, 20);
    h_truthMatchIxn_PrimPhotonMultiplicity = new TH1D("h_truthMatchIxn_PrimPhotonMultiplicity", "h_truthMatchIxn_PrimPhotonMultiplicity", 20, 0, 20);
    h_truthMatchIxn_SecPhotonMultiplicity = new TH1D("h_truthMatchIxn_SecPhotonMultiplicity", "h_truthMatchIxn_SecPhotonMultiplicity", 20, 0, 20);

    h_truthMatchIxn_Reco_PrimShowerResponseMult = new TH2D("h_truthMatchIxn_Reco_PrimShowerResponseMult", "h_truthMatchIxn_Reco_PrimShowerResponseMult", 20, 0, 20, 20, 0, 20);
    h_truthMatchIxn_Reco_SecShowerResponseMult = new TH2D("h_truthMatchIxn_Reco_SecShowerResponseMult", "h_truthMatchIxn_Reco_SecShowerResponseMult", 20, 0, 20, 20, 0, 20);
    h_truthMatchIxn_Reco_PrimElectronResponseMult = new TH2D("h_truthMatchIxn_Reco_PrimElectronResponseMult", "h_truthMatchIxn_Reco_PrimElectronResponseMult", 20, 0, 20, 20, 0, 20);
    h_truthMatchIxn_Reco_SecElectronResponseMult = new TH2D("h_truthMatchIxn_Reco_SecElectronResponseMult", "h_truthMatchIxn_Reco_SecElectronResponseMult", 20, 0, 20, 20, 0, 20);
    h_truthMatchIxn_Reco_PrimPhotonResponseMult = new TH2D("h_truthMatchIxn_Reco_PrimPhotonResponseMult", "h_truthMatchIxn_Reco_PrimPhotonResponseMult", 20, 0, 20, 20, 0, 20);
    h_truthMatchIxn_Reco_SecPhotonResponseMult = new TH2D("h_truthMatchIxn_Reco_SecPhotonResponseMult", "h_truthMatchIxn_Reco_SecPhotonResponseMult", 20, 0, 20, 20, 0, 20);

    // Reco Pi0, Shower, Electron, and Photon Multiplicity for Best-Matched Interaction
    // TRUTH MATCHED SIGNAL ONLY
    h_truthMatchIxn_TruthSignal_RecoPrimShowerMult = new TH1D("h_truthMatchIxn_TruthSignal_RecoPrimShowerMult", "h_truthMatchIxn_TruthSignal_RecoPrimShowerMult", 20, 0, 20);
    h_truthMatchIxn_TruthSignal_RecoSecShowerMult = new TH1D("h_truthMatchIxn_TruthSignal_RecoSecShowerMult", "h_truthMatchIxn_TruthSignal_RecoSecShowerMult", 20, 0, 20);
    h_truthMatchIxn_TruthSignal_RecoPrimElectronMult = new TH1D("h_truthMatchIxn_TruthSignal_RecoPrimElectronMult", "h_truthMatchIxn_TruthSignal_RecoPrimElectronMult", 20, 0, 20);
    h_truthMatchIxn_TruthSignal_RecoSecElectronMult = new TH1D("h_truthMatchIxn_TruthSignal_RecoSecElectronMult", "h_truthMatchIxn_TruthSignal_RecoSecElectronMult", 20, 0, 20);
    h_truthMatchIxn_TruthSignal_RecoPrimPhotonMult = new TH1D("h_truthMatchIxn_TruthSignal_RecoPrimPhotonMult", "h_truthMatchIxn_TruthSignal_RecoPrimPhotonMult", 20, 0, 20);
    h_truthMatchIxn_TruthSignal_RecoSecPhotonMult = new TH1D("h_truthMatchIxn_TruthSignal_RecoSecPhotonMult", "h_truthMatchIxn_TruthSignal_RecoSecPhotonMult", 20, 0, 20);

    // TRUTH MATCHED SIGNAL CC QE ONLY
    h_truthMatchIxn_TruthSignalCCQE_RecoPrimShowerMult = new TH1D("h_truthMatchIxn_TruthSignalCCQE_RecoPrimShowerMult", "h_truthMatchIxn_TruthSignalCCQE_RecoPrimShowerMult", 20, 0, 20);
    h_truthMatchIxn_TruthSignalCCQE_RecoSecShowerMult = new TH1D("h_truthMatchIxn_TruthSignalCCQE_RecoSecShowerMult", "h_truthMatchIxn_TruthSignalCCQE_RecoSecShowerMult", 20, 0, 20);
    h_truthMatchIxn_TruthSignalCCQE_RecoPrimElectronMult = new TH1D("h_truthMatchIxn_TruthSignalCCQE_RecoPrimElectronMult", "h_truthMatchIxn_TruthSignalCCQE_RecoPrimElectronMult", 20, 0, 20);
    h_truthMatchIxn_TruthSignalCCQE_RecoSecElectronMult = new TH1D("h_truthMatchIxn_TruthSignalCCQE_RecoSecElectronMult", "h_truthMatchIxn_TruthSignalCCQE_RecoSecElectronMult", 20, 0, 20);
    h_truthMatchIxn_TruthSignalCCQE_RecoPrimPhotonMult = new TH1D("h_truthMatchIxn_TruthSignalCCQE_RecoPrimPhotonMult", "h_truthMatchIxn_TruthSignalCCQE_RecoPrimPhotonMult", 20, 0, 20);
    h_truthMatchIxn_TruthSignalCCQE_RecoSecPhotonMult = new TH1D("h_truthMatchIxn_TruthSignalCCQE_RecoSecPhotonMult", "h_truthMatchIxn_TruthSignalCCQE_RecoSecPhotonMult", 20, 0, 20);

    // TRUTH MATCHED SIGNAL CC MEC ONLY
    h_truthMatchIxn_TruthSignalCCMEC_RecoPrimShowerMult = new TH1D("h_truthMatchIxn_TruthSignalCCMEC_RecoPrimShowerMult", "h_truthMatchIxn_TruthSignalCCMEC_RecoPrimShowerMult", 20, 0, 20);
    h_truthMatchIxn_TruthSignalCCMEC_RecoSecShowerMult = new TH1D("h_truthMatchIxn_TruthSignalCCMEC_RecoSecShowerMult", "h_truthMatchIxn_TruthSignalCCMEC_RecoSecShowerMult", 20, 0, 20);
    h_truthMatchIxn_TruthSignalCCMEC_RecoPrimElectronMult = new TH1D("h_truthMatchIxn_TruthSignalCCMEC_RecoPrimElectronMult", "h_truthMatchIxn_TruthSignalCCMEC_RecoPrimElectronMult", 20, 0, 20);
    h_truthMatchIxn_TruthSignalCCMEC_RecoSecElectronMult = new TH1D("h_truthMatchIxn_TruthSignalCCMEC_RecoSecElectronMult", "h_truthMatchIxn_TruthSignalCCMEC_RecoSecElectronMult", 20, 0, 20);
    h_truthMatchIxn_TruthSignalCCMEC_RecoPrimPhotonMult = new TH1D("h_truthMatchIxn_TruthSignalCCMEC_RecoPrimPhotonMult", "h_truthMatchIxn_TruthSignalCCMEC_RecoPrimPhotonMult", 20, 0, 20);
    h_truthMatchIxn_TruthSignalCCMEC_RecoSecPhotonMult = new TH1D("h_truthMatchIxn_TruthSignalCCMEC_RecoSecPhotonMult", "h_truthMatchIxn_TruthSignalCCMEC_RecoSecPhotonMult", 20, 0, 20);

    // TRUTH MATCHED SIGNAL CC DIS ONLY
    h_truthMatchIxn_TruthSignalCCDIS_RecoPrimShowerMult = new TH1D("h_truthMatchIxn_TruthSignalCCDIS_RecoPrimShowerMult", "h_truthMatchIxn_TruthSignalCCDIS_RecoPrimShowerMult", 20, 0, 20);
    h_truthMatchIxn_TruthSignalCCDIS_RecoSecShowerMult = new TH1D("h_truthMatchIxn_TruthSignalCCDIS_RecoSecShowerMult", "h_truthMatchIxn_TruthSignalCCDIS_RecoSecShowerMult", 20, 0, 20);
    h_truthMatchIxn_TruthSignalCCDIS_RecoPrimElectronMult = new TH1D("h_truthMatchIxn_TruthSignalCCDIS_RecoPrimElectronMult", "h_truthMatchIxn_TruthSignalCCDIS_RecoPrimElectronMult", 20, 0, 20);
    h_truthMatchIxn_TruthSignalCCDIS_RecoSecElectronMult = new TH1D("h_truthMatchIxn_TruthSignalCCDIS_RecoSecElectronMult", "h_truthMatchIxn_TruthSignalCCDIS_RecoSecElectronMult", 20, 0, 20);
    h_truthMatchIxn_TruthSignalCCDIS_RecoPrimPhotonMult = new TH1D("h_truthMatchIxn_TruthSignalCCDIS_RecoPrimPhotonMult", "h_truthMatchIxn_TruthSignalCCDIS_RecoPrimPhotonMult", 20, 0, 20);
    h_truthMatchIxn_TruthSignalCCDIS_RecoSecPhotonMult = new TH1D("h_truthMatchIxn_TruthSignalCCDIS_RecoSecPhotonMult", "h_truthMatchIxn_TruthSignalCCDIS_RecoSecPhotonMult", 20, 0, 20);

    // TRUTH MATCHED SIGNAL CC RES ONLY
    h_truthMatchIxn_TruthSignalCCRES_RecoPrimShowerMult = new TH1D("h_truthMatchIxn_TruthSignalCCRES_RecoPrimShowerMult", "h_truthMatchIxn_TruthSignalCCRES_RecoPrimShowerMult", 20, 0, 20);
    h_truthMatchIxn_TruthSignalCCRES_RecoSecShowerMult = new TH1D("h_truthMatchIxn_TruthSignalCCRES_RecoSecShowerMult", "h_truthMatchIxn_TruthSignalCCRES_RecoSecShowerMult", 20, 0, 20);
    h_truthMatchIxn_TruthSignalCCRES_RecoPrimElectronMult = new TH1D("h_truthMatchIxn_TruthSignalCCRES_RecoPrimElectronMult", "h_truthMatchIxn_TruthSignalCCRES_RecoPrimElectronMult", 20, 0, 20);
    h_truthMatchIxn_TruthSignalCCRES_RecoSecElectronMult = new TH1D("h_truthMatchIxn_TruthSignalCCRES_RecoSecElectronMult", "h_truthMatchIxn_TruthSignalCCRES_RecoSecElectronMult", 20, 0, 20);
    h_truthMatchIxn_TruthSignalCCRES_RecoPrimPhotonMult = new TH1D("h_truthMatchIxn_TruthSignalCCRES_RecoPrimPhotonMult", "h_truthMatchIxn_TruthSignalCCRES_RecoPrimPhotonMult", 20, 0, 20);
    h_truthMatchIxn_TruthSignalCCRES_RecoSecPhotonMult = new TH1D("h_truthMatchIxn_TruthSignalCCRES_RecoSecPhotonMult", "h_truthMatchIxn_TruthSignalCCRES_RecoSecPhotonMult", 20, 0, 20);

    // TRUTH MATCHED SIGNAL CC COH ONLY
    h_truthMatchIxn_TruthSignalCCCOH_RecoPrimShowerMult = new TH1D("h_truthMatchIxn_TruthSignalCCCOH_RecoPrimShowerMult", "h_truthMatchIxn_TruthSignalCCCOH_RecoPrimShowerMult", 20, 0, 20);
    h_truthMatchIxn_TruthSignalCCCOH_RecoSecShowerMult = new TH1D("h_truthMatchIxn_TruthSignalCCCOH_RecoSecShowerMult", "h_truthMatchIxn_TruthSignalCCCOH_RecoSecShowerMult", 20, 0, 20);
    h_truthMatchIxn_TruthSignalCCCOH_RecoPrimElectronMult = new TH1D("h_truthMatchIxn_TruthSignalCCCOH_RecoPrimElectronMult", "h_truthMatchIxn_TruthSignalCCCOH_RecoPrimElectronMult", 20, 0, 20);
    h_truthMatchIxn_TruthSignalCCCOH_RecoSecElectronMult = new TH1D("h_truthMatchIxn_TruthSignalCCCOH_RecoSecElectronMult", "h_truthMatchIxn_TruthSignalCCCOH_RecoSecElectronMult", 20, 0, 20);
    h_truthMatchIxn_TruthSignalCCCOH_RecoPrimPhotonMult = new TH1D("h_truthMatchIxn_TruthSignalCCCOH_RecoPrimPhotonMult", "h_truthMatchIxn_TruthSignalCCCOH_RecoPrimPhotonMult", 20, 0, 20);
    h_truthMatchIxn_TruthSignalCCCOH_RecoSecPhotonMult = new TH1D("h_truthMatchIxn_TruthSignalCCCOH_RecoSecPhotonMult", "h_truthMatchIxn_TruthSignalCCCOH_RecoSecPhotonMult", 20, 0, 20);

    // Pi0, Shower, Electron, and Photon Multiplicity for Best-Matched Interaction
    // TRUTH MATCHED BKG ONLY
    h_truthMatchIxn_TruthBkg_RecoPrimShowerMult = new TH1D("h_truthMatchIxn_TruthBkg_RecoPrimShowerMult", "h_truthMatchIxn_TruthBkg_RecoPrimShowerMult", 20, 0, 20);
    h_truthMatchIxn_TruthBkg_RecoSecShowerMult = new TH1D("h_truthMatchIxn_TruthBkg_RecoSecShowerMult", "h_truthMatchIxn_TruthBkg_RecoSecShowerMult", 20, 0, 20);
    h_truthMatchIxn_TruthBkg_RecoPrimElectronMult = new TH1D("h_truthMatchIxn_TruthBkg_RecoPrimElectronMult", "h_truthMatchIxn_TruthBkg_RecoPrimElectronMult", 20, 0, 20);
    h_truthMatchIxn_TruthBkg_RecoSecElectronMult = new TH1D("h_truthMatchIxn_TruthBkg_RecoSecElectronMult", "h_truthMatchIxn_TruthBkg_RecoSecElectronMult", 20, 0, 20);
    h_truthMatchIxn_TruthBkg_RecoPrimPhotonMult = new TH1D("h_truthMatchIxn_TruthBkg_RecoPrimPhotonMult", "h_truthMatchIxn_TruthBkg_RecoPrimPhotonMult", 20, 0, 20);
    h_truthMatchIxn_TruthBkg_RecoSecPhotonMult = new TH1D("h_truthMatchIxn_TruthBkg_RecoSecPhotonMult", "h_truthMatchIxn_TruthBkg_RecoSecPhotonMult", 20, 0, 20);

    // TRUTH MATCHED BKG NC ONLY
    h_truthMatchIxn_TruthBkgNC_RecoPrimShowerMult = new TH1D("h_truthMatchIxn_TruthBkgNC_RecoPrimShowerMult", "h_truthMatchIxn_TruthBkgNC_RecoPrimShowerMult", 20, 0, 20);
    h_truthMatchIxn_TruthBkgNC_RecoSecShowerMult = new TH1D("h_truthMatchIxn_TruthBkgNC_RecoSecShowerMult", "h_truthMatchIxn_TruthBkgNC_RecoSecShowerMult", 20, 0, 20);
    h_truthMatchIxn_TruthBkgNC_RecoPrimElectronMult = new TH1D("h_truthMatchIxn_TruthBkgNC_RecoPrimElectronMult", "h_truthMatchIxn_TruthBkgNC_RecoPrimElectronMult", 20, 0, 20);
    h_truthMatchIxn_TruthBkgNC_RecoSecElectronMult = new TH1D("h_truthMatchIxn_TruthBkgNC_RecoSecElectronMult", "h_truthMatchIxn_TruthBkgNC_RecoSecElectronMult", 20, 0, 20);
    h_truthMatchIxn_TruthBkgNC_RecoPrimPhotonMult = new TH1D("h_truthMatchIxn_TruthBkgNC_RecoPrimPhotonMult", "h_truthMatchIxn_TruthBkgNC_RecoPrimPhotonMult", 20, 0, 20);
    h_truthMatchIxn_TruthBkgNC_RecoSecPhotonMult = new TH1D("h_truthMatchIxn_TruthBkgNC_RecoSecPhotonMult", "h_truthMatchIxn_TruthBkgNC_RecoSecPhotonMult", 20, 0, 20);

    // TRUTH MATCHED BKG ROCK ONLY
    h_truthMatchIxn_TruthBkgROCK_RecoPrimShowerMult = new TH1D("h_truthMatchIxn_TruthBkgROCK_RecoPrimShowerMult", "h_truthMatchIxn_TruthBkgROCK_RecoPrimShowerMult", 20, 0, 20);
    h_truthMatchIxn_TruthBkgROCK_RecoSecShowerMult = new TH1D("h_truthMatchIxn_TruthBkgROCK_RecoSecShowerMult", "h_truthMatchIxn_TruthBkgROCK_RecoSecShowerMult", 20, 0, 20);
    h_truthMatchIxn_TruthBkgROCK_RecoPrimElectronMult = new TH1D("h_truthMatchIxn_TruthBkgROCK_RecoPrimElectronMult", "h_truthMatchIxn_TruthBkgROCK_RecoPrimElectronMult", 20, 0, 20);
    h_truthMatchIxn_TruthBkgROCK_RecoSecElectronMult = new TH1D("h_truthMatchIxn_TruthBkgROCK_RecoSecElectronMult", "h_truthMatchIxn_TruthBkgROCK_RecoSecElectronMult", 20, 0, 20);
    h_truthMatchIxn_TruthBkgROCK_RecoPrimPhotonMult = new TH1D("h_truthMatchIxn_TruthBkgROCK_RecoPrimPhotonMult", "h_truthMatchIxn_TruthBkgROCK_RecoPrimPhotonMult", 20, 0, 20);
    h_truthMatchIxn_TruthBkgROCK_RecoSecPhotonMult = new TH1D("h_truthMatchIxn_TruthBkgROCK_RecoSecPhotonMult", "h_truthMatchIxn_TruthBkgROCK_RecoSecPhotonMult", 20, 0, 20);

    // TRUTH MATCHED BKG CC ONLY
    h_truthMatchIxn_TruthBkgCC_RecoPrimShowerMult = new TH1D("h_truthMatchIxn_TruthBkgCC_RecoPrimShowerMult", "h_truthMatchIxn_TruthBkgCC_RecoPrimShowerMult", 20, 0, 20);
    h_truthMatchIxn_TruthBkgCC_RecoSecShowerMult = new TH1D("h_truthMatchIxn_TruthBkgCC_RecoSecShowerMult", "h_truthMatchIxn_TruthBkgCC_RecoSecShowerMult", 20, 0, 20);
    h_truthMatchIxn_TruthBkgCC_RecoPrimElectronMult = new TH1D("h_truthMatchIxn_TruthBkgCC_RecoPrimElectronMult", "h_truthMatchIxn_TruthBkgCC_RecoPrimElectronMult", 20, 0, 20);
    h_truthMatchIxn_TruthBkgCC_RecoSecElectronMult = new TH1D("h_truthMatchIxn_TruthBkgCC_RecoSecElectronMult", "h_truthMatchIxn_TruthBkgCC_RecoSecElectronMult", 20, 0, 20);
    h_truthMatchIxn_TruthBkgCC_RecoPrimPhotonMult = new TH1D("h_truthMatchIxn_TruthBkgCC_RecoPrimPhotonMult", "h_truthMatchIxn_TruthBkgCC_RecoPrimPhotonMult", 20, 0, 20);
    h_truthMatchIxn_TruthBkgCC_RecoSecPhotonMult = new TH1D("h_truthMatchIxn_TruthBkgCC_RecoSecPhotonMult", "h_truthMatchIxn_TruthBkgCC_RecoSecPhotonMult", 20, 0, 20);

    // Ixn Mode for Best-Matched Interaction
    h_truthMatchIxn_IxnMode = new TH1D("h_truthMatchIxn_IxnMode", "h_truthMatchIxn_IxnMode", 1100, 0, 1100);
  
    // Is best matched truth ixn rock? (0 = no, 1 = yes)
    h_truthMatchIxn_IsRock = new TH1D("h_truthMatchIxn_IsRock", "h_truthMatchIxn_IsRock", 2, 0, 2);

    // Check if best matched truth ixn is CC/NC, nue/numu if it's rock
    h_truthMatchIxn_TruthBkgROCK_isCC = new TH1D("h_truthMatchIxn_TruthBkgROCK_isCC", "h_truthMatchIxn_TruthBkgROCK_isCC", 2, 0, 2);
    h_truthMatchIxn_TruthBkgROCK_nuPDG = new TH1D("h_truthMatchIxn_TruthBkgROCK_nuPDG", "h_truthMatchIxn_TruthBkgROCK_nuPDG", 2, 0, 2);

    // Check if best matched truth ixn is CC/NC, nue/numu if it's not rock
    h_truthMatchIxn_TruthBkgNONROCK_isCC = new TH1D("h_truthMatchIxn_TruthBkgNONROCK_isCC", "h_truthMatchIxn_TruthBkgNONROCK_isCC", 2, 0, 2);
    h_truthMatchIxn_TruthBkgNONROCK_nuPDG = new TH1D("h_truthMatchIxn_TruthBkgNONROCK_nuPDG", "h_truthMatchIxn_TruthBkgNONROCK_nuPDG", 2, 0, 2);
}


//------------------------------------------
// Histogram Filling Methods
//------------------------------------------    
void TruthMatchedHists::FillTruthMatchMx2TrackInfo(MatchedInteractionSummary& matchSummary,
                                                   RecoInteractionSummary& recoSummary)
{
    // Info from match w/ Mx2
    if (matchSummary.truthMatchMx2TrackisPrimary){
        h_truthMatchPrim_mx2TrackPDG->Fill(matchSummary.truthMatchMx2TrackPDG);
        h_truthMatchPrim_mx2TrackE->Fill(matchSummary.truthMatchMx2TrackE);
        h_truthMatchPrim_mx2TrackCosL->Fill(matchSummary.truthMatchMx2TrackCosL);

        h_truthMatchPrim_mx2TrackRecoDiffLArStartX->Fill(matchSummary.truthMatchMx2TrackLArStartPosX-
                                                         recoSummary.mx2MatchLArStartPosX);
        h_truthMatchPrim_mx2TrackRecoDiffLArStartY->Fill(matchSummary.truthMatchMx2TrackLArStartPosY-
                                                         recoSummary.mx2MatchLArStartPosY);
        h_truthMatchPrim_mx2TrackRecoDiffLArStartZ->Fill(matchSummary.truthMatchMx2TrackLArStartPosZ-
                                                         recoSummary.mx2MatchLArStartPosZ);

        h_truthMatchIxn_Reco_responseMx2MatchPrimCosL->Fill(recoSummary.muonCosL, matchSummary.truthMatchMx2TrackCosL);
    }
    else {
        h_truthMatchSec_mx2TrackPDG->Fill(matchSummary.truthMatchMx2TrackPDG);
        h_truthMatchSec_mx2TrackE->Fill(matchSummary.truthMatchMx2TrackE);
        h_truthMatchSec_mx2TrackCosL->Fill(matchSummary.truthMatchMx2TrackCosL);

        h_truthMatchSec_mx2TrackRecoDiffLArStartX->Fill(matchSummary.truthMatchMx2TrackLArStartPosX-
                                                        recoSummary.mx2MatchLArStartPosX);
        h_truthMatchSec_mx2TrackRecoDiffLArStartY->Fill(matchSummary.truthMatchMx2TrackLArStartPosY-
                                                        recoSummary.mx2MatchLArStartPosY);
        h_truthMatchSec_mx2TrackRecoDiffLArStartZ->Fill(matchSummary.truthMatchMx2TrackLArStartPosZ-
                                                        recoSummary.mx2MatchLArStartPosZ);

        h_truthMatchIxn_Reco_responseMx2MatchSecCosL->Fill(recoSummary.muonCosL, matchSummary.truthMatchMx2TrackCosL);
    }

    // Info from best-matched truth interaction muon
    if (matchSummary.passesMx2 == true) h_truthMatchIxn_MuonCosL->Fill(matchSummary.truthSummaryforBestMatch.muonCosL);
    h_truthMatchIxn_MuonCosL_unbinned->Fill(matchSummary.truthSummaryforBestMatch.muonCosL);
    h_truthMatchIxn_MuonCosL_zoomOut->Fill(matchSummary.truthSummaryforBestMatch.muonCosL);
    h_truthMatchIxn_MuonE->Fill(matchSummary.truthSummaryforBestMatch.muonEnergy);
    h_truthMatchIxn_Reco_responseMuonCosL->Fill(recoSummary.muonCosL, matchSummary.truthSummaryforBestMatch.muonCosL);

}

void TruthMatchedHists::FillTruthMatchDiffTruthRecoVertex(MatchedInteractionSummary& matchSummary)
{
    h_truthMatch_diffTruthRecoVertex->Fill(matchSummary.diffTruthRecoVertex);
}

void TruthMatchedHists::FillTruthMatchIxnShowerMultiplicity(MatchedInteractionSummary& matchSummary, RecoInteractionSummary& recoSummary)
{
    // Fill histograms for shower multiplicity for best matched interaction (all interactions passing reco cuts)
    h_truthMatchIxn_PrimPi0Multiplicity->Fill(matchSummary.truthSummaryforBestMatch.nPrimPi0);
    h_truthMatchIxn_SecPi0Multiplicity->Fill(matchSummary.truthSummaryforBestMatch.nSecPi0);

    int total_prim_showers_truthMatch = matchSummary.truthSummaryforBestMatch.nPrimElectrons + matchSummary.truthSummaryforBestMatch.nPrimPhotons;
    int total_sec_showers_truthMatch  = matchSummary.truthSummaryforBestMatch.nSecElectrons + matchSummary.truthSummaryforBestMatch.nSecPhotons;
    h_truthMatchIxn_PrimShowerMultiplicity->Fill(total_prim_showers_truthMatch);
    h_truthMatchIxn_SecShowerMultiplicity->Fill(total_sec_showers_truthMatch);

    h_truthMatchIxn_PrimElectronMultiplicity->Fill(matchSummary.truthSummaryforBestMatch.nPrimElectrons);
    h_truthMatchIxn_SecElectronMultiplicity->Fill(matchSummary.truthSummaryforBestMatch.nSecElectrons);

    h_truthMatchIxn_PrimPhotonMultiplicity->Fill(matchSummary.truthSummaryforBestMatch.nPrimPhotons);
    h_truthMatchIxn_SecPhotonMultiplicity->Fill(matchSummary.truthSummaryforBestMatch.nSecPhotons);

    // Response in Reco Shower Multiplicity for Prim and Sec Showers for Best Matched Interaction (all interactions passing reco cuts)
    int total_prim_showers_reco = recoSummary.nPrimElectrons + recoSummary.nPrimPhotons;
    int total_sec_showers_reco = recoSummary.nSecElectrons + recoSummary.nSecPhotons;

    h_truthMatchIxn_Reco_PrimShowerResponseMult->Fill(total_prim_showers_reco, total_prim_showers_truthMatch);
    h_truthMatchIxn_Reco_SecShowerResponseMult->Fill(total_sec_showers_reco, total_sec_showers_truthMatch);
    h_truthMatchIxn_Reco_PrimElectronResponseMult->Fill(recoSummary.nPrimElectrons, matchSummary.truthSummaryforBestMatch.nPrimElectrons);
    h_truthMatchIxn_Reco_SecElectronResponseMult->Fill(recoSummary.nSecElectrons, matchSummary.truthSummaryforBestMatch.nSecElectrons);
    h_truthMatchIxn_Reco_PrimPhotonResponseMult->Fill(recoSummary.nPrimPhotons, matchSummary.truthSummaryforBestMatch.nPrimPhotons);
    h_truthMatchIxn_Reco_SecPhotonResponseMult->Fill(recoSummary.nSecPhotons, matchSummary.truthSummaryforBestMatch.nSecPhotons);

    // Fill reco shower multiplicity for best matched ixn truth signal/bkg histograms
    if (matchSummary.passesLArCuts && matchSummary.passesMx2 == true)
    {
        h_truthMatchIxn_TruthSignal_RecoPrimShowerMult->Fill(total_prim_showers_reco);
        h_truthMatchIxn_TruthSignal_RecoSecShowerMult->Fill(total_sec_showers_reco);
        h_truthMatchIxn_TruthSignal_RecoPrimElectronMult->Fill(recoSummary.nPrimElectrons);
        h_truthMatchIxn_TruthSignal_RecoSecElectronMult->Fill(recoSummary.nSecElectrons);
        h_truthMatchIxn_TruthSignal_RecoPrimPhotonMult->Fill(recoSummary.nPrimPhotons);
        h_truthMatchIxn_TruthSignal_RecoSecPhotonMult->Fill(recoSummary.nSecPhotons);

        // All signal ixns are CC by definition -- check ixn mode
        if (matchSummary.truthSummaryforBestMatch.ixnMode == 1 || matchSummary.truthSummaryforBestMatch.ixnMode == 1001) {
          h_truthMatchIxn_TruthSignalCCQE_RecoPrimShowerMult->Fill(total_prim_showers_reco);
          h_truthMatchIxn_TruthSignalCCQE_RecoSecShowerMult->Fill(total_sec_showers_reco);
          h_truthMatchIxn_TruthSignalCCQE_RecoPrimElectronMult->Fill(recoSummary.nPrimElectrons);
          h_truthMatchIxn_TruthSignalCCQE_RecoSecElectronMult->Fill(recoSummary.nSecElectrons);
          h_truthMatchIxn_TruthSignalCCQE_RecoPrimPhotonMult->Fill(recoSummary.nPrimPhotons);
          h_truthMatchIxn_TruthSignalCCQE_RecoSecPhotonMult->Fill(recoSummary.nSecPhotons);
        }
        else if (matchSummary.truthSummaryforBestMatch.ixnMode == 10) {
          h_truthMatchIxn_TruthSignalCCMEC_RecoPrimShowerMult->Fill(total_prim_showers_reco);
          h_truthMatchIxn_TruthSignalCCMEC_RecoSecShowerMult->Fill(total_sec_showers_reco);
          h_truthMatchIxn_TruthSignalCCMEC_RecoPrimElectronMult->Fill(recoSummary.nPrimElectrons);
          h_truthMatchIxn_TruthSignalCCMEC_RecoSecElectronMult->Fill(recoSummary.nSecElectrons);
          h_truthMatchIxn_TruthSignalCCMEC_RecoPrimPhotonMult->Fill(recoSummary.nPrimPhotons);
          h_truthMatchIxn_TruthSignalCCMEC_RecoSecPhotonMult->Fill(recoSummary.nSecPhotons);
        }
        else if (matchSummary.truthSummaryforBestMatch.ixnMode == 3) {
          h_truthMatchIxn_TruthSignalCCDIS_RecoPrimShowerMult->Fill(total_prim_showers_reco);
          h_truthMatchIxn_TruthSignalCCDIS_RecoSecShowerMult->Fill(total_sec_showers_reco);
          h_truthMatchIxn_TruthSignalCCDIS_RecoPrimElectronMult->Fill(recoSummary.nPrimElectrons);
          h_truthMatchIxn_TruthSignalCCDIS_RecoSecElectronMult->Fill(recoSummary.nSecElectrons);
          h_truthMatchIxn_TruthSignalCCDIS_RecoPrimPhotonMult->Fill(recoSummary.nPrimPhotons);
          h_truthMatchIxn_TruthSignalCCDIS_RecoSecPhotonMult->Fill(recoSummary.nSecPhotons);
        }
        else if (matchSummary.truthSummaryforBestMatch.ixnMode == 4) {
          h_truthMatchIxn_TruthSignalCCRES_RecoPrimShowerMult->Fill(total_prim_showers_reco);
          h_truthMatchIxn_TruthSignalCCRES_RecoSecShowerMult->Fill(total_sec_showers_reco);
          h_truthMatchIxn_TruthSignalCCRES_RecoPrimElectronMult->Fill(recoSummary.nPrimElectrons);
          h_truthMatchIxn_TruthSignalCCRES_RecoSecElectronMult->Fill(recoSummary.nSecElectrons);
          h_truthMatchIxn_TruthSignalCCRES_RecoPrimPhotonMult->Fill(recoSummary.nPrimPhotons);
          h_truthMatchIxn_TruthSignalCCRES_RecoSecPhotonMult->Fill(recoSummary.nSecPhotons);
        }
        else if (matchSummary.truthSummaryforBestMatch.ixnMode == 5) {
          h_truthMatchIxn_TruthSignalCCCOH_RecoPrimShowerMult->Fill(total_prim_showers_reco);
          h_truthMatchIxn_TruthSignalCCCOH_RecoSecShowerMult->Fill(total_sec_showers_reco);
          h_truthMatchIxn_TruthSignalCCCOH_RecoPrimElectronMult->Fill(recoSummary.nPrimElectrons);
          h_truthMatchIxn_TruthSignalCCCOH_RecoSecElectronMult->Fill(recoSummary.nSecElectrons);
          h_truthMatchIxn_TruthSignalCCCOH_RecoPrimPhotonMult->Fill(recoSummary.nPrimPhotons);
          h_truthMatchIxn_TruthSignalCCCOH_RecoSecPhotonMult->Fill(recoSummary.nSecPhotons);
        }
    }
    else {
        h_truthMatchIxn_TruthBkg_RecoPrimShowerMult->Fill(total_prim_showers_reco);
        h_truthMatchIxn_TruthBkg_RecoSecShowerMult->Fill(total_sec_showers_reco);
        h_truthMatchIxn_TruthBkg_RecoPrimElectronMult->Fill(recoSummary.nPrimElectrons);
        h_truthMatchIxn_TruthBkg_RecoSecElectronMult->Fill(recoSummary.nSecElectrons);
        h_truthMatchIxn_TruthBkg_RecoPrimPhotonMult->Fill(recoSummary.nPrimPhotons);
        h_truthMatchIxn_TruthBkg_RecoSecPhotonMult->Fill(recoSummary.nSecPhotons);

        if (matchSummary.isRockIxn == true)
        {
            h_truthMatchIxn_TruthBkgROCK_RecoPrimShowerMult->Fill(total_prim_showers_reco);
            h_truthMatchIxn_TruthBkgROCK_RecoSecShowerMult->Fill(total_sec_showers_reco);
            h_truthMatchIxn_TruthBkgROCK_RecoPrimElectronMult->Fill(recoSummary.nPrimElectrons);
            h_truthMatchIxn_TruthBkgROCK_RecoSecElectronMult->Fill(recoSummary.nSecElectrons);
            h_truthMatchIxn_TruthBkgROCK_RecoPrimPhotonMult->Fill(recoSummary.nPrimPhotons);
            h_truthMatchIxn_TruthBkgROCK_RecoSecPhotonMult->Fill(recoSummary.nSecPhotons);
        }
        else if (matchSummary.truthSummaryforBestMatch.iscc == true)
        {
            h_truthMatchIxn_TruthBkgCC_RecoPrimShowerMult->Fill(total_prim_showers_reco);
            h_truthMatchIxn_TruthBkgCC_RecoSecShowerMult->Fill(total_sec_showers_reco);
            h_truthMatchIxn_TruthBkgCC_RecoPrimElectronMult->Fill(recoSummary.nPrimElectrons);
            h_truthMatchIxn_TruthBkgCC_RecoSecElectronMult->Fill(recoSummary.nSecElectrons);
            h_truthMatchIxn_TruthBkgCC_RecoPrimPhotonMult->Fill(recoSummary.nPrimPhotons);
            h_truthMatchIxn_TruthBkgCC_RecoSecPhotonMult->Fill(recoSummary.nSecPhotons);
        }
        else
        {
            h_truthMatchIxn_TruthBkgNC_RecoPrimShowerMult->Fill(total_prim_showers_reco);
            h_truthMatchIxn_TruthBkgNC_RecoSecShowerMult->Fill(total_sec_showers_reco);
            h_truthMatchIxn_TruthBkgNC_RecoPrimElectronMult->Fill(recoSummary.nPrimElectrons);
            h_truthMatchIxn_TruthBkgNC_RecoSecElectronMult->Fill(recoSummary.nSecElectrons);
            h_truthMatchIxn_TruthBkgNC_RecoPrimPhotonMult->Fill(recoSummary.nPrimPhotons);
            h_truthMatchIxn_TruthBkgNC_RecoSecPhotonMult->Fill(recoSummary.nSecPhotons);
        }
    }
}

// Fill histogram for interaction mode for best-matched interaction (all interactions passing reco cuts)
void TruthMatchedHists::FillTruthMatchIxnIxnMode(MatchedInteractionSummary& matchSummary)
{
  h_truthMatchIxn_IxnMode->Fill(matchSummary.truthSummaryforBestMatch.ixnMode);
}

// Fill histogram for counting true rock ixns selected in reco sel
void TruthMatchedHists::FillTruthMatchIxnIsRock(MatchedInteractionSummary& matchSummary)
{
  h_truthMatchIxn_IsRock->Fill(matchSummary.isRockIxn);
}

// Fill histogram for counting CC/NC, nu PDG for rock and non-rock BKG
void TruthMatchedHists::FillTruthMatchIxnIsBkgCCnuPDG(MatchedInteractionSummary& matchSummary)
{
  if (matchSummary.passesLArCuts && matchSummary.passesMx2 == true) {
    return; // only want to fill for interactions that fail reco cuts (i.e. are in bkg)
  }
  else if (matchSummary.isRockIxn == true)
  { 
    h_truthMatchIxn_TruthBkgROCK_isCC->Fill(matchSummary.truthSummaryforBestMatch.iscc);
    h_truthMatchIxn_TruthBkgROCK_nuPDG->Fill(matchSummary.truthSummaryforBestMatch.nuPDG);
  }
  else {
    h_truthMatchIxn_TruthBkgNONROCK_isCC->Fill(matchSummary.truthSummaryforBestMatch.iscc);
    h_truthMatchIxn_TruthBkgNONROCK_nuPDG->Fill(matchSummary.truthSummaryforBestMatch.nuPDG);
  }
}

//------------------------------------------
// Histogram Writing Method
//------------------------------------------   
void TruthMatchedHists::Write(TDirectory* dir)
{
    dir->cd();

    h_truthMatchPrim_mx2TrackPDG->Write();
    h_truthMatchSec_mx2TrackPDG->Write();

    h_truthMatchPrim_mx2TrackE->Write();
    h_truthMatchSec_mx2TrackE->Write();

    h_truthMatchPrim_mx2TrackCosL->Write();
    h_truthMatchSec_mx2TrackCosL->Write();
    h_truthMatchIxn_Reco_responseMx2MatchPrimCosL->Write();
    h_truthMatchIxn_Reco_responseMx2MatchSecCosL->Write();

    h_truthMatchIxn_MuonCosL->Write();
    h_truthMatchIxn_MuonCosL_unbinned->Write();
    h_truthMatchIxn_MuonCosL_zoomOut->Write();
    h_truthMatchIxn_MuonE->Write();
    h_truthMatchIxn_Reco_responseMuonCosL->Write();

    h_truthMatchPrim_mx2TrackRecoDiffLArStartX->Write();
    h_truthMatchSec_mx2TrackRecoDiffLArStartX->Write();
    h_truthMatchPrim_mx2TrackRecoDiffLArStartY->Write();
    h_truthMatchSec_mx2TrackRecoDiffLArStartY->Write();
    h_truthMatchPrim_mx2TrackRecoDiffLArStartZ->Write();
    h_truthMatchSec_mx2TrackRecoDiffLArStartZ->Write();

    h_truthMatch_diffTruthRecoVertex->Write();

    h_truthMatchIxn_PrimPi0Multiplicity->Write();
    h_truthMatchIxn_SecPi0Multiplicity->Write();
    h_truthMatchIxn_PrimShowerMultiplicity->Write();
    h_truthMatchIxn_SecShowerMultiplicity->Write();
    h_truthMatchIxn_PrimElectronMultiplicity->Write();
    h_truthMatchIxn_SecElectronMultiplicity->Write();
    h_truthMatchIxn_PrimPhotonMultiplicity->Write();
    h_truthMatchIxn_SecPhotonMultiplicity->Write();

    h_truthMatchIxn_Reco_PrimShowerResponseMult->Write();
    h_truthMatchIxn_Reco_SecShowerResponseMult->Write();
    h_truthMatchIxn_Reco_PrimElectronResponseMult->Write();
    h_truthMatchIxn_Reco_SecElectronResponseMult->Write();
    h_truthMatchIxn_Reco_PrimPhotonResponseMult->Write();
    h_truthMatchIxn_Reco_SecPhotonResponseMult->Write();

    h_truthMatchIxn_TruthSignal_RecoPrimShowerMult->Write();
    h_truthMatchIxn_TruthSignal_RecoSecShowerMult->Write();
    h_truthMatchIxn_TruthSignal_RecoPrimElectronMult->Write();
    h_truthMatchIxn_TruthSignal_RecoSecElectronMult->Write();
    h_truthMatchIxn_TruthSignal_RecoPrimPhotonMult->Write();
    h_truthMatchIxn_TruthSignal_RecoSecPhotonMult->Write();

    h_truthMatchIxn_TruthSignalCCQE_RecoPrimShowerMult->Write();
    h_truthMatchIxn_TruthSignalCCQE_RecoSecShowerMult->Write();
    h_truthMatchIxn_TruthSignalCCQE_RecoPrimElectronMult->Write();
    h_truthMatchIxn_TruthSignalCCQE_RecoSecElectronMult->Write();
    h_truthMatchIxn_TruthSignalCCQE_RecoPrimPhotonMult->Write();
    h_truthMatchIxn_TruthSignalCCQE_RecoSecPhotonMult->Write();

    h_truthMatchIxn_TruthSignalCCMEC_RecoPrimShowerMult->Write();
    h_truthMatchIxn_TruthSignalCCMEC_RecoSecShowerMult->Write();
    h_truthMatchIxn_TruthSignalCCMEC_RecoPrimElectronMult->Write();
    h_truthMatchIxn_TruthSignalCCMEC_RecoSecElectronMult->Write();
    h_truthMatchIxn_TruthSignalCCMEC_RecoPrimPhotonMult->Write();
    h_truthMatchIxn_TruthSignalCCMEC_RecoSecPhotonMult->Write();

    h_truthMatchIxn_TruthSignalCCDIS_RecoPrimShowerMult->Write();
    h_truthMatchIxn_TruthSignalCCDIS_RecoSecShowerMult->Write();
    h_truthMatchIxn_TruthSignalCCDIS_RecoPrimElectronMult->Write();
    h_truthMatchIxn_TruthSignalCCDIS_RecoSecElectronMult->Write();
    h_truthMatchIxn_TruthSignalCCDIS_RecoPrimPhotonMult->Write();
    h_truthMatchIxn_TruthSignalCCDIS_RecoSecPhotonMult->Write();

    h_truthMatchIxn_TruthSignalCCRES_RecoPrimShowerMult->Write();
    h_truthMatchIxn_TruthSignalCCRES_RecoSecShowerMult->Write();
    h_truthMatchIxn_TruthSignalCCRES_RecoPrimElectronMult->Write();
    h_truthMatchIxn_TruthSignalCCRES_RecoSecElectronMult->Write();
    h_truthMatchIxn_TruthSignalCCRES_RecoPrimPhotonMult->Write();
    h_truthMatchIxn_TruthSignalCCRES_RecoSecPhotonMult->Write();

    h_truthMatchIxn_TruthSignalCCCOH_RecoPrimShowerMult->Write();
    h_truthMatchIxn_TruthSignalCCCOH_RecoSecShowerMult->Write();
    h_truthMatchIxn_TruthSignalCCCOH_RecoPrimElectronMult->Write();
    h_truthMatchIxn_TruthSignalCCCOH_RecoSecElectronMult->Write();
    h_truthMatchIxn_TruthSignalCCCOH_RecoPrimPhotonMult->Write();
    h_truthMatchIxn_TruthSignalCCCOH_RecoSecPhotonMult->Write();

    h_truthMatchIxn_TruthBkg_RecoPrimShowerMult->Write();
    h_truthMatchIxn_TruthBkg_RecoSecShowerMult->Write();
    h_truthMatchIxn_TruthBkg_RecoPrimElectronMult->Write();
    h_truthMatchIxn_TruthBkg_RecoSecElectronMult->Write();
    h_truthMatchIxn_TruthBkg_RecoPrimPhotonMult->Write();
    h_truthMatchIxn_TruthBkg_RecoSecPhotonMult->Write();

    h_truthMatchIxn_TruthBkgNC_RecoPrimShowerMult->Write();
    h_truthMatchIxn_TruthBkgNC_RecoSecShowerMult->Write();
    h_truthMatchIxn_TruthBkgNC_RecoPrimElectronMult->Write();
    h_truthMatchIxn_TruthBkgNC_RecoSecElectronMult->Write();
    h_truthMatchIxn_TruthBkgNC_RecoPrimPhotonMult->Write();
    h_truthMatchIxn_TruthBkgNC_RecoSecPhotonMult->Write();

    h_truthMatchIxn_TruthBkgROCK_RecoPrimShowerMult->Write();
    h_truthMatchIxn_TruthBkgROCK_RecoSecShowerMult->Write();
    h_truthMatchIxn_TruthBkgROCK_RecoPrimElectronMult->Write();
    h_truthMatchIxn_TruthBkgROCK_RecoSecElectronMult->Write();
    h_truthMatchIxn_TruthBkgROCK_RecoPrimPhotonMult->Write();
    h_truthMatchIxn_TruthBkgROCK_RecoSecPhotonMult->Write();

    h_truthMatchIxn_TruthBkgCC_RecoPrimShowerMult->Write();
    h_truthMatchIxn_TruthBkgCC_RecoSecShowerMult->Write();
    h_truthMatchIxn_TruthBkgCC_RecoPrimElectronMult->Write();
    h_truthMatchIxn_TruthBkgCC_RecoSecElectronMult->Write();
    h_truthMatchIxn_TruthBkgCC_RecoPrimPhotonMult->Write();
    h_truthMatchIxn_TruthBkgCC_RecoSecPhotonMult->Write();

    h_truthMatchIxn_IxnMode->Write();

    h_truthMatchIxn_IsRock->Write();

    h_truthMatchIxn_TruthBkgROCK_isCC->Write();
    h_truthMatchIxn_TruthBkgROCK_nuPDG->Write();

    h_truthMatchIxn_TruthBkgNONROCK_isCC->Write();
    h_truthMatchIxn_TruthBkgNONROCK_nuPDG->Write();

}