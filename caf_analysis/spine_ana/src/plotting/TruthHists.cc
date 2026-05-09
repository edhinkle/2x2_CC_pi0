#pragma once

#include "plotting/TruthHists.h"

// Define all truth histograms 
TruthHists::TruthHists()
{
    
    // True interactions above KE threshold per spill
    h_true_ixnsAboveKEThresholdPerSpill = new TH1D("true_ixnsAboveKEThresholdPerSpill", "true_ixnsAboveKEThresholdPerSpill", 10, 0, 10);

    // Edges for CosL histograms
    Double_t edges[7] = {0.91, 0.96, 0.98, 0.9887, 0.994, 0.9974, 1};

    // True muon kinematics
    h_true_CosL = new TH1D("true_CosL", "true_CosL", 6, edges);
    h_true_CosL_zoomOut = new TH1D("true_CosL_zoomOut", "true_CosL_zoomOut", 50, 0.8, 1);
    h_true_CosLNumubar_zoomOut = new TH1D(
      "trueCosLNumubar_zoomOut", "trueCosLNumubar_zoomOut", 50, 0.8, 1);
    h_true_CosLNumu_zoomOut = new TH1D(
      "trueCosLNumu_zoomOut", "trueCosLNumu_zoomOut", 50, 0.8, 1);
    h_true_Elep = new TH1D("true_Elep", "true_Elep", 50, 0, 20);

    // True neutrino energy
    h_true_Enu = new TH1D("true_Enu", "true_Enu", 100, 0, 20);

    // Bins and edges for multiplicity histograms
    int binsMult = 10;
    Double_t edgesMult[11] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

    h_true_nPi0 = new TH1D("true_nPi0", "true_nPi0", binsMult, edgesMult);

}


//------------------------------------------
// Histogram Filling Methods
//------------------------------------------    
// Fill histograms for number of interactions above KE threshold per spill
void TruthHists::FillInteractionsAboveKEThresholdPerSpill(int n)
{
    h_true_ixnsAboveKEThresholdPerSpill->Fill(n);
}

void TruthHists::FillMuonKinematics(double cosL, double Elep, bool Numubar)
{
    h_true_CosL_zoomOut->Fill(cosL);
    if(Numubar)
        h_true_CosLNumubar_zoomOut->Fill(cosL);
    if(!Numubar)
        h_true_CosLNumu_zoomOut->Fill(cosL);

    if (cosL > 0.91 && Elep > 1) {
        h_true_CosL->Fill(cosL);
        h_true_Elep->Fill(Elep);
    }
}

void TruthHists::FillEnu(double Enu)
{
    h_true_Enu->Fill(Enu);
}

void TruthHists::FillPi0Multiplicity(int nPi0)
{
    h_true_nPi0->Fill(nPi0);
}

void TruthHists::Write(TDirectory* dir)
{
    dir->cd();

    h_true_ixnsAboveKEThresholdPerSpill->Write();

    h_true_CosL->Write();
    h_true_CosL_zoomOut->Write();
    h_true_CosLNumubar_zoomOut->Write();
    h_true_CosLNumu_zoomOut->Write();
    h_true_Elep->Write();

    h_true_Enu->Write();

    h_true_nPi0->Write();

}