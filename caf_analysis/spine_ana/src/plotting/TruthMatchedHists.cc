#pragma once

#include "plotting/TruthMatchedHists.h"
#include "TDirectory.h"

// Define all truth match histograms 
TruthMatchedHists::TruthMatchedHists()
{
    // Truth match mx2Track histograms
    // PDG
    h_truthMatchPrim_mx2TrackPDG = new TH1D("h_truthMatchPrim_mx2TrackPDG", "h_truthMatchPrim_mx2TrackPDG", 6000, -3000, 3000);
    h_truthMatchSec_mx2TrackPDG = new TH1D("h_truthMatchSec_mx2TrackPDG", "h_truthMatchSec_mx2TrackPDG", 6000, -3000, 3000)
    
    // Energy
    h_truthMatchPrim_mx2TrackE = new TH1D("h_truthMatchPrim_mx2TrackE", "h_truthMatchPrim_mx2TrackE", 50, 0, 20);
    h_truthMatchSec_mx2TrackE = new TH1D("h_truthMatchSec_mx2TrackE", "h_truthMatchSec_mx2TrackE", 50, 0, 20);

    // Cos Angle w/ beam
    h_truthMatchPrim_mx2TrackCosL = new TH1D("h_truthMatchPrim_mx2TrackCosL", "h_truthMatchPrim_mx2TrackCosL", 100, -1, 1);
    h_truthMatchSec_mx2TrackCosL = new TH1D("h_truthMatchSec_mx2TrackCosL", "h_truthMatchSec_mx2TrackCosL", 100, -1, 1);

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


}


//------------------------------------------
// Histogram Filling Methods
//------------------------------------------    
void TruthMatchedHists::FillTruthMatchMx2TrackInfo(MatchedInteractionSummary& matchSummary,
                                                   RecoInteractionSummary& recoSummary)
{
    if (matchSummary.isPrimary){
        h_truthMatchPrim_mx2TrackPDG->Fill(matchSummary.truthMatchMx2TrackPDG);
        h_truthMatchPrim_mx2TrackE->Fill(matchSummary.truthMatchMx2TrackE);
        h_truthMatchPrim_mx2TrackCosL->Fill(matchSummary.truthMatchMx2TrackCosL);

        h_truthMatchPrim_mx2TrackRecoDiffLArStartX->Fill(matchSummary.truthMatchMx2TrackLArStartPosX-
                                                         recoSummary.mx2MatchLArStartPosX);
        h_truthMatchPrim_mx2TrackRecoDiffLArStartY->Fill(matchSummary.truthMatchMx2TrackLArStartPosY-
                                                         recoSummary.mx2MatchLArStartPosY);
        h_truthMatchPrim_mx2TrackRecoDiffLArStartZ->Fill(matchSummary.truthMatchMx2TrackLArStartPosZ-
                                                         recoSummary.mx2MatchLArStartPosZ);
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

    h_truthMatchPrim_mx2TrackRecoDiffLArStartX->Write();
    h_truthMatchSec_mx2TrackRecoDiffLArStartX->Write();
    h_truthMatchPrim_mx2TrackRecoDiffLArStartY->Write();
    h_truthMatchSec_mx2TrackRecoDiffLArStartY->Write();
    h_truthMatchPrim_mx2TrackRecoDiffLArStartZ->Write();
    h_truthMatchSec_mx2TrackRecoDiffLArStartZ->Write();

}