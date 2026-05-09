#pragma once

#include "plotting/TruthHists.h"

// Define all truth histograms 
TruthHists::TruthHists(double flux_nom)
{
    
    // True interactions above KE threshold per spill
    h_true_ixnsAboveKEThresholdPerSpill = new TH1D("true_ixnsAboveKEThresholdPerSpill", "true_ixnsAboveKEThresholdPerSpill", 10, 0, 10);

    // Cutflow histogram
    h_true_cutflow = new TH1D("h_cutflow", "Truth Selection;Selection Stage;Events", 10, 0, 10);
    h_true_cutflow->GetXaxis()->SetBinLabel(1, "Nu on Ar");
    h_true_cutflow->GetXaxis()->SetBinLabel(2, "Fiducial");


    // True neutrino energy
    h_true_Enu = new TH1D("true_Enu", "true_Enu", 100, 0, 20);

    // Bins and edges for multiplicity histograms
    int binsMult = 10;
    Double_t edgesMult[11] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

    // True particle and shower multiplicity histograms
    

    Double_t edges[7] = {0.91, 0.96, 0.98, 0.9887, 0.994, 0.9974, 1};

    h_trueCosL = new TH1D("trueCosL", Form("%0.08f", flux_nom), 6, edges);
    h_trueCosL_zoomOut = new TH1D("trueCosL_zoomOut", "", 50, 0.8, 1);
    h_trueEnu = new TH1D("trueEnu", "trueEnu", 100, 0, 20);

    h_true_mult = new TH1D("true_mult", "", binsMult, edgesMult);
    h_true_multShort = new TH1D("true_multShort", "", binsMult, edgesMult);
    h_true_multLong  = new TH1D("true_multLong", "", binsMult, edgesMult);

    h_nPi = new TH1D("nPi", "", binsMult, edgesMult);
    h_escapePi = new TH1D("escapePi", "", 20, 0, 20);
}


//------------------------------------------
// Histogram Filling Methods
//------------------------------------------    
// Fill histograms for number of interactions above KE threshold per spill
void TruthHists::FillInteractionsAboveKEThresholdPerSpill(int n)
{
    h_true_ixnsAboveKEThresholdPerSpill->Fill(n);
}

void TruthHists::FillMuonKinematics(double cosL, double Elep, double Enu)
{
    h_trueCosL_zoomOut->Fill(cosL);

    if (cosL > 0.91 && Elep > 1) {
        h_trueCosL->Fill(cosL);
        h_trueEnu->Fill(Enu);
    }
}

void TruthHists::FillMultiplicity(int nTrk, int nShort, int nLong)
{
    h_true_mult->Fill(nTrk);
    h_true_multShort->Fill(nShort);
    h_true_multLong->Fill(nLong);
}

void TruthHists::FillPionInfo(int nPi, int escapingPi)
{
    h_nPi->Fill(nPi);
    h_escapePi->Fill(escapingPi);
}

void TruthHists::Write(TDirectory* dir)
{
    dir->cd();

    h_true_ixnsAboveKEThresholdPerSpill->Write();


    h_trueCosL->Write();
    h_trueCosL_zoomOut->Write();
    h_trueEnu->Write();
    h_true_mult->Write();
    h_true_multShort->Write();
    h_true_multLong->Write();
    h_nPi->Write();
    h_escapePi->Write();
}