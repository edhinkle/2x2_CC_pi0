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

    // Truth - Reco Diff in Vertex Position for Best Truth Match
    h_truthMatch_diffTruthRecoVertex = 
      new TH1D("h_truthMatch_diffTruthRecoVertex", "h_truthMatch_diffTruthRecoVertex", 30, 0, 15);

    // Truth match interaction pi0 multiplicity and shower multiplicity histograms
    // Bins and edges for multiplicity histograms
    int binsMult = 10;
    Double_t edgesMult[11] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    h_truthMatchIxn_PrimPi0Multiplicity = new TH1D("h_truthMatchIxn_PrimPi0Multiplicity", "h_truthMatchIxn_PrimPi0Multiplicity", binsMult, edgesMult);
    h_truthMatchIxn_SecPi0Multiplicity = new TH1D("h_truthMatchIxn_SecPi0Multiplicity", "h_truthMatchIxn_SecPi0Multiplicity", binsMult, edgesMult);
    h_truthMatchIxn_PrimShowerMultiplicity = new TH1D("h_truthMatchIxn_PrimShowerMultiplicity", "h_truthMatchIxn_PrimShowerMultiplicity", binsMult, edgesMult);
    h_truthMatchIxn_SecShowerMultiplicity = new TH1D("h_truthMatchIxn_SecShowerMultiplicity", "h_truthMatchIxn_SecShowerMultiplicity", binsMult, edgesMult);
    h_truthMatchIxn_PrimElectronMultiplicity = new TH1D("h_truthMatchIxn_PrimElectronMultiplicity", "h_truthMatchIxn_PrimElectronMultiplicity", binsMult, edgesMult);
    h_truthMatchIxn_SecElectronMultiplicity = new TH1D("h_truthMatchIxn_SecElectronMultiplicity", "h_truthMatchIxn_SecElectronMultiplicity", binsMult, edgesMult);
    h_truthMatchIxn_PrimPhotonMultiplicity = new TH1D("h_truthMatchIxn_PrimPhotonMultiplicity", "h_truthMatchIxn_PrimPhotonMultiplicity", binsMult, edgesMult);
    h_truthMatchIxn_SecPhotonMultiplicity = new TH1D("h_truthMatchIxn_SecPhotonMultiplicity", "h_truthMatchIxn_SecPhotonMultiplicity", binsMult, edgesMult);
}


//------------------------------------------
// Histogram Filling Methods
//------------------------------------------    
void TruthMatchedHists::FillTruthMatchMx2TrackInfo(MatchedInteractionSummary& matchSummary,
                                                   RecoInteractionSummary& recoSummary)
{
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

void TruthMatchedHists::FillTruthMatchDiffTruthRecoVertex(MatchedInteractionSummary& matchSummary)
{
    h_truthMatch_diffTruthRecoVertex->Fill(matchSummary.diffTruthRecoVertex);
}

void TruthMatchedHists::FillTruthMatchIxnShowerMultiplicity(MatchedInteractionSummary& matchSummary)
{
    h_truthMatchIxn_PrimPi0Multiplicity->Fill(matchSummary.truthSummaryforBestMatch.nPrimPi0);
    h_truthMatchIxn_SecPi0Multiplicity->Fill(matchSummary.truthSummaryforBestMatch.nSecPi0);

    int total_prim_showers = matchSummary.truthSummaryforBestMatch.nPrimElectron + matchSummary.truthSummaryforBestMatch.nPrimPhoton;
    int total_sec_showers = matchSummary.truthSummaryforBestMatch.nSecElectron + matchSummary.truthSummaryforBestMatch.nSecPhoton;
    h_truthMatchIxn_PrimShowerMultiplicity->Fill(total_prim_showers);
    h_truthMatchIxn_SecShowerMultiplicity->Fill(total_sec_showers);

    h_truthMatchIxn_PrimElectronMultiplicity->Fill(matchSummary.truthSummaryforBestMatch.nPrimElectron);
    h_truthMatchIxn_SecElectronMultiplicity->Fill(matchSummary.truthSummaryforBestMatch.nSecElectron);

    h_truthMatchIxn_PrimPhotonMultiplicity->Fill(matchSummary.truthSummaryforBestMatch.nPrimPhoton);
    h_truthMatchIxn_SecPhotonMultiplicity->Fill(matchSummary.truthSummaryforBestMatch.nSecPhoton);
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

    h_truthMatch_diffTruthRecoVertex->Write();

    h_truthMatchIxn_PrimPi0Multiplicity->Write();
    h_truthMatchIxn_SecPi0Multiplicity->Write();
    h_truthMatchIxn_PrimShowerMultiplicity->Write();
    h_truthMatchIxn_SecShowerMultiplicity->Write();
    h_truthMatchIxn_PrimElectronMultiplicity->Write();
    h_truthMatchIxn_SecElectronMultiplicity->Write();
    h_truthMatchIxn_PrimPhotonMultiplicity->Write();
    h_truthMatchIxn_SecPhotonMultiplicity->Write();

}