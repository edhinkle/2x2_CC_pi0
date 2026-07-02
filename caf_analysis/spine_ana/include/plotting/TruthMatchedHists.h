#pragma once
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "analysis/RecoInteractionSummary.h"
#include "analysis/TruthInteractionSummary.h"
#include "analysis/MatchedInteractionSummary.h"
#include "cuts/Mx2MatchResult.h"

class TruthMatchedHists {
public:
    TruthMatchedHists(double nominalIntegratedFlux);

    //--------------------------------------------------------
    // Filling histograms
    //--------------------------------------------------------
    void FillTruthMatchMx2TrackInfo(MatchedInteractionSummary& matchSummary,
                                    RecoInteractionSummary& recoSummary);

    void FillTruthMatchDiffTruthRecoVertex(MatchedInteractionSummary& matchSummary);
    void FillTruthMatchIxnShowerMultiplicity(MatchedInteractionSummary& matchSummary,
                                             RecoInteractionSummary& recoSummary);
    void FillTruthMatchIxnIxnMode(MatchedInteractionSummary& matchSummary);
    void FillTruthMatchIxnIsRock(MatchedInteractionSummary& matchSummary);
    void FillTruthMatchIxnIsBkgCCnuPDG(MatchedInteractionSummary& matchSummary);

    //--------------------------------------------------------
    // Get methods for copying histograms for systematics
    //--------------------------------------------------------
    const TH1D* GetTruthMatchCosL() const { return h_truthMatchIxn_MuonCosL_zoomOut; }
    const TH2D* GetRespCosL() const { return h_truthMatchIxn_Reco_responseMuonCosL; }

    //--------------------------------------------------------
    // Write histograms to file
    //--------------------------------------------------------
    void Write(TDirectory* dir);

private:

    //--------------------------------------------------------
    // Histograms
    //--------------------------------------------------------
    // Mx2 Matched Track Info (Mx2 match track)
    TH1D *h_truthMatchPrim_mx2TrackPDG;
    TH1D *h_truthMatchSec_mx2TrackPDG;

    TH1D *h_truthMatchPrim_mx2TrackE;
    TH1D *h_truthMatchSec_mx2TrackE;

    TH1D *h_truthMatchPrim_mx2TrackCosL;
    TH1D *h_truthMatchSec_mx2TrackCosL;

    // Mx2 Match Info (Best match ixn muon)
    TH2D *h_truthMatchIxn_Reco_responseMx2MatchPrimCosL;
    TH2D *h_truthMatchIxn_Reco_responseMx2MatchSecCosL;

    TH1D *h_truthMatchIxn_MuonCosL;
    TH1D *h_truthMatchIxn_MuonCosL_unbinned;
    TH1D *h_truthMatchIxn_MuonCosL_zoomOut;
    TH1D *h_truthMatchIxn_MuonE;
    TH2D *h_truthMatchIxn_Reco_responseMuonCosL;

    // Mx2 Matched Track Info (Mx2 match track)
    TH1D *h_truthMatchPrim_mx2TrackRecoDiffLArStartX;
    TH1D *h_truthMatchSec_mx2TrackRecoDiffLArStartX;
    TH1D *h_truthMatchPrim_mx2TrackRecoDiffLArStartY;
    TH1D *h_truthMatchSec_mx2TrackRecoDiffLArStartY;
    TH1D *h_truthMatchPrim_mx2TrackRecoDiffLArStartZ;
    TH1D *h_truthMatchSec_mx2TrackRecoDiffLArStartZ;

    // Difference between reco and truth vertex for best-matched interaction
    TH1D *h_truthMatch_diffTruthRecoVertex;

    // Pi0, Shower, Electron, and Photon Multiplicity for Best-Matched Interaction
    // Filled for all ixns passing reco cuts
    TH1D *h_truthMatchIxn_PrimPi0Multiplicity;
    TH1D *h_truthMatchIxn_SecPi0Multiplicity;
    TH1D *h_truthMatchIxn_PrimShowerMultiplicity;
    TH1D *h_truthMatchIxn_SecShowerMultiplicity;
    TH1D *h_truthMatchIxn_PrimElectronMultiplicity;
    TH1D *h_truthMatchIxn_SecElectronMultiplicity;
    TH1D *h_truthMatchIxn_PrimPhotonMultiplicity;
    TH1D *h_truthMatchIxn_SecPhotonMultiplicity;

    // Response for truth match ixn/reco ixn shower, electron, photon multiplicity
    TH2D *h_truthMatchIxn_Reco_PrimShowerResponseMult;
    TH2D *h_truthMatchIxn_Reco_SecShowerResponseMult;
    TH2D *h_truthMatchIxn_Reco_PrimElectronResponseMult;
    TH2D *h_truthMatchIxn_Reco_SecElectronResponseMult;
    TH2D *h_truthMatchIxn_Reco_PrimPhotonResponseMult;
    TH2D *h_truthMatchIxn_Reco_SecPhotonResponseMult;

    // Reco Pi0, Shower, Electron, and Photon Multiplicity for Best-Matched Interaction
    // TRUTH MATCHED SIGNAL ONLY
    TH1D *h_truthMatchIxn_TruthSignal_RecoPrimShowerMult;
    TH1D *h_truthMatchIxn_TruthSignal_RecoSecShowerMult;
    TH1D *h_truthMatchIxn_TruthSignal_RecoPrimElectronMult;
    TH1D *h_truthMatchIxn_TruthSignal_RecoSecElectronMult;
    TH1D *h_truthMatchIxn_TruthSignal_RecoPrimPhotonMult;
    TH1D *h_truthMatchIxn_TruthSignal_RecoSecPhotonMult;

    // Reco Pi0, Shower, Electron, and Photon Multiplicity for Best-Matched Interaction
    // TRUTH MATCHED SIGNAL CC QE ONLY
    TH1D *h_truthMatchIxn_TruthSignalCCQE_RecoPrimShowerMult;
    TH1D *h_truthMatchIxn_TruthSignalCCQE_RecoSecShowerMult;
    TH1D *h_truthMatchIxn_TruthSignalCCQE_RecoPrimElectronMult;
    TH1D *h_truthMatchIxn_TruthSignalCCQE_RecoSecElectronMult;
    TH1D *h_truthMatchIxn_TruthSignalCCQE_RecoPrimPhotonMult;
    TH1D *h_truthMatchIxn_TruthSignalCCQE_RecoSecPhotonMult;

    // TRUTH MATCHED SIGNAL CC MEC ONLY
    TH1D *h_truthMatchIxn_TruthSignalCCMEC_RecoPrimShowerMult;
    TH1D *h_truthMatchIxn_TruthSignalCCMEC_RecoSecShowerMult;
    TH1D *h_truthMatchIxn_TruthSignalCCMEC_RecoPrimElectronMult;
    TH1D *h_truthMatchIxn_TruthSignalCCMEC_RecoSecElectronMult;
    TH1D *h_truthMatchIxn_TruthSignalCCMEC_RecoPrimPhotonMult;
    TH1D *h_truthMatchIxn_TruthSignalCCMEC_RecoSecPhotonMult;

    // TRUTH MATCHED SIGNAL CC DIS ONLY
    TH1D *h_truthMatchIxn_TruthSignalCCDIS_RecoPrimShowerMult;
    TH1D *h_truthMatchIxn_TruthSignalCCDIS_RecoSecShowerMult;
    TH1D *h_truthMatchIxn_TruthSignalCCDIS_RecoPrimElectronMult;
    TH1D *h_truthMatchIxn_TruthSignalCCDIS_RecoSecElectronMult;
    TH1D *h_truthMatchIxn_TruthSignalCCDIS_RecoPrimPhotonMult;
    TH1D *h_truthMatchIxn_TruthSignalCCDIS_RecoSecPhotonMult;

    // TRUTH MATCHED SIGNAL CC RES ONLY
    TH1D *h_truthMatchIxn_TruthSignalCCRES_RecoPrimShowerMult;
    TH1D *h_truthMatchIxn_TruthSignalCCRES_RecoSecShowerMult;
    TH1D *h_truthMatchIxn_TruthSignalCCRES_RecoPrimElectronMult;
    TH1D *h_truthMatchIxn_TruthSignalCCRES_RecoSecElectronMult;
    TH1D *h_truthMatchIxn_TruthSignalCCRES_RecoPrimPhotonMult;
    TH1D *h_truthMatchIxn_TruthSignalCCRES_RecoSecPhotonMult;

    // TRUTH MATCHED SIGNAL CC COH ONLY
    TH1D *h_truthMatchIxn_TruthSignalCCCOH_RecoPrimShowerMult;
    TH1D *h_truthMatchIxn_TruthSignalCCCOH_RecoSecShowerMult;
    TH1D *h_truthMatchIxn_TruthSignalCCCOH_RecoPrimElectronMult;
    TH1D *h_truthMatchIxn_TruthSignalCCCOH_RecoSecElectronMult;
    TH1D *h_truthMatchIxn_TruthSignalCCCOH_RecoPrimPhotonMult;
    TH1D *h_truthMatchIxn_TruthSignalCCCOH_RecoSecPhotonMult;

    // Pi0, Shower, Electron, and Photon Multiplicity for Best-Matched Interaction
    // TRUTH MATCHED BKG ONLY
    TH1D *h_truthMatchIxn_TruthBkg_RecoPrimShowerMult;
    TH1D *h_truthMatchIxn_TruthBkg_RecoSecShowerMult;
    TH1D *h_truthMatchIxn_TruthBkg_RecoPrimElectronMult;
    TH1D *h_truthMatchIxn_TruthBkg_RecoSecElectronMult;
    TH1D *h_truthMatchIxn_TruthBkg_RecoPrimPhotonMult;
    TH1D *h_truthMatchIxn_TruthBkg_RecoSecPhotonMult;

    // TRUTH MATCHED BKG NC ONLY
    TH1D *h_truthMatchIxn_TruthBkgNC_RecoPrimShowerMult;
    TH1D *h_truthMatchIxn_TruthBkgNC_RecoSecShowerMult;
    TH1D *h_truthMatchIxn_TruthBkgNC_RecoPrimElectronMult;
    TH1D *h_truthMatchIxn_TruthBkgNC_RecoSecElectronMult;
    TH1D *h_truthMatchIxn_TruthBkgNC_RecoPrimPhotonMult;
    TH1D *h_truthMatchIxn_TruthBkgNC_RecoSecPhotonMult;

    // TRUTH MATCHED BKG ROCK ONLY
    TH1D *h_truthMatchIxn_TruthBkgROCK_RecoPrimShowerMult;
    TH1D *h_truthMatchIxn_TruthBkgROCK_RecoSecShowerMult;
    TH1D *h_truthMatchIxn_TruthBkgROCK_RecoPrimElectronMult;
    TH1D *h_truthMatchIxn_TruthBkgROCK_RecoSecElectronMult;
    TH1D *h_truthMatchIxn_TruthBkgROCK_RecoPrimPhotonMult;
    TH1D *h_truthMatchIxn_TruthBkgROCK_RecoSecPhotonMult;

    // TRUTH MATCHED BKG CC ONLY
    TH1D *h_truthMatchIxn_TruthBkgCC_RecoPrimShowerMult;
    TH1D *h_truthMatchIxn_TruthBkgCC_RecoSecShowerMult;
    TH1D *h_truthMatchIxn_TruthBkgCC_RecoPrimElectronMult;
    TH1D *h_truthMatchIxn_TruthBkgCC_RecoSecElectronMult;
    TH1D *h_truthMatchIxn_TruthBkgCC_RecoPrimPhotonMult;
    TH1D *h_truthMatchIxn_TruthBkgCC_RecoSecPhotonMult;

    // Ixn mode of best matched interaction
    TH1D *h_truthMatchIxn_IxnMode;

    // Is best matched truth ixn rock? (0 = no, 1 = yes)
    TH1D *h_truthMatchIxn_IsRock;

    // Check if best matched truth ixn is CC/NC, nue/numu if it's rock
    TH1D *h_truthMatchIxn_TruthBkgROCK_isCC;
    TH1D *h_truthMatchIxn_TruthBkgROCK_nuPDG;

    // Check if best matched truth ixn is CC/NC, nue/numu if it's not rock
    TH1D *h_truthMatchIxn_TruthBkgNONROCK_isCC;
    TH1D *h_truthMatchIxn_TruthBkgNONROCK_nuPDG;

    double fNominalIntegratedFlux;

};