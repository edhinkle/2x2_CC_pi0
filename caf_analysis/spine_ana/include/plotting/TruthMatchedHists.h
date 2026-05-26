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
    TruthMatchedHists();

    //------------------------------------------
    // Filling histograms
    //------------------------------------------
    void FillTruthMatchMx2TrackInfo(MatchedInteractionSummary& matchSummary,
                                    RecoInteractionSummary& recoSummary);

    void FillTruthMatchDiffTruthRecoVertex(MatchedInteractionSummary& matchSummary);
    void FillTruthMatchIxnShowerMultiplicity(MatchedInteractionSummary& matchSummary);
    void Write(TDirectory* dir);

private:

    //------------------------------------------
    // Histograms
    //------------------------------------------ 
    TH1D *h_truthMatchPrim_mx2TrackPDG;
    TH1D *h_truthMatchSec_mx2TrackPDG;

    TH1D *h_truthMatchPrim_mx2TrackE;
    TH1D *h_truthMatchSec_mx2TrackE;

    TH1D *h_truthMatchPrim_mx2TrackCosL;
    TH1D *h_truthMatchSec_mx2TrackCosL;

    TH1D *h_truthMatchPrim_mx2TrackRecoDiffLArStartX;
    TH1D *h_truthMatchSec_mx2TrackRecoDiffLArStartX;
    TH1D *h_truthMatchPrim_mx2TrackRecoDiffLArStartY;
    TH1D *h_truthMatchSec_mx2TrackRecoDiffLArStartY;
    TH1D *h_truthMatchPrim_mx2TrackRecoDiffLArStartZ;
    TH1D *h_truthMatchSec_mx2TrackRecoDiffLArStartZ;

    TH1D *h_truthMatch_diffTruthRecoVertex;

    TH1D *h_truthMatchIxn_PrimPi0Multiplicity;
    TH1D *h_truthMatchIxn_SecPi0Multiplicity;
    TH1D *h_truthMatchIxn_PrimShowerMultiplicity;
    TH1D *h_truthMatchIxn_SecShowerMultiplicity;
    TH1D *h_truthMatchIxn_PrimElectronMultiplicity;
    TH1D *h_truthMatchIxn_SecElectronMultiplicity;
    TH1D *h_truthMatchIxn_PrimPhotonMultiplicity;
    TH1D *h_truthMatchIxn_SecPhotonMultiplicity;

};