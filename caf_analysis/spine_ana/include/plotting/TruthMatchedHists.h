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

};