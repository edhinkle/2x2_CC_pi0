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
    h_true_CosL = new TH1D("true_CosL", "true_CosL", 6, edges); // Signal region (passes Mx2 cuts)
    h_true_CosL_zoomOut = new TH1D("true_CosL_zoomOut", "true_CosL_zoomOut", 50, 0.8, 1); // All events passing initial cuts
    h_true_CosLNumubar_zoomOut = new TH1D(
      "trueCosLNumubar_zoomOut", "trueCosLNumubar_zoomOut", 50, 0.8, 1); // Numubar events passing initial cuts
    h_true_CosLNumu_zoomOut = new TH1D(
      "trueCosLNumu_zoomOut", "trueCosLNumu_zoomOut", 50, 0.8, 1); // Numu events passing initial cuts
    h_true_Elep = new TH1D("true_Elep", "true_Elep", 50, 0, 20); // Signal region (passes Mx2 cuts)

    // True neutrino energy
    h_true_Enu = new TH1D("true_Enu", "true_Enu", 100, 0, 20);

    // Bins and edges for multiplicity histograms
    int binsMult = 10;
    Double_t edgesMult[11] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

    // pi0 multiplicity
    h_true_nPrimPi0 = new TH1D("true_nPrimPi0", "true_nPrimPi0", binsMult, edgesMult);
    h_true_nSecPi0_preMx2 = new TH1D("true_nSecPi0_preMx2", "true_nSecPi0_preMx2", binsMult, edgesMult);
    h_true_nSecPi0_postMx2 = new TH1D("true_nSecPi0_postMx2", "true_nSecPi0_postMx2", binsMult, edgesMult);

    // shower multiplicity
    h_true_nPrimElectron_preMx2 = new TH1D("true_nPrimElectron_preMx2", "true_nPrimElectron_preMx2", binsMult, edgesMult);
    h_true_nPrimElectron_postMx2 = new TH1D("true_nPrimElectron_postMx2", "true_nPrimElectron_postMx2", binsMult, edgesMult);
    h_true_nSecElectron_preMx2 = new TH1D("true_nSecElectron_preMx2", "true_nSecElectron_preMx2", binsMult, edgesMult);
    h_true_nSecElectron_postMx2 = new TH1D("true_nSecElectron_postMx2", "true_nSecElectron_postMx2", binsMult, edgesMult);

    h_true_nPrimPhoton_preMx2 = new TH1D("true_nPrimPhoton_preMx2", "true_nPrimPhoton_preMx2", binsMult, edgesMult);
    h_true_nPrimPhoton_postMx2 = new TH1D("true_nPrimPhoton_postMx2", "true_nPrimPhoton_postMx2", binsMult, edgesMult);
    h_true_nSecPhoton_preMx2 = new TH1D("true_nSecPhoton_preMx2", "true_nSecPhoton_preMx2", binsMult, edgesMult);
    h_true_nSecPhoton_postMx2 = new TH1D("true_nSecPhoton_postMx2", "true_nSecPhoton_postMx2", binsMult, edgesMult);

    h_true_nPrimShower_preMx2 = new TH1D("true_nPrimShower_preMx2", "true_nPrimShower_preMx2", binsMult, edgesMult);
    h_true_nPrimShower_postMx2 = new TH1D("true_nPrimShower_postMx2", "true_nPrimShower_postMx2", binsMult, edgesMult);
    h_true_nSecShower_preMx2 = new TH1D("true_nSecShower_preMx2", "true_nSecShower_preMx2", binsMult, edgesMult);
    h_true_nSecShower_postMx2 = new TH1D("true_nSecShower_postMx2", "true_nSecShower_postMx2", binsMult, edgesMult);

    h_true_IxnMode = new TH1D("true_IxnMode", "true_IxnMode");

}


//------------------------------------------
// Histogram Filling Methods
//------------------------------------------    
// Fill histograms for number of interactions above KE threshold per spill
void TruthHists::FillInteractionsAboveKEThresholdPerSpill(int n)
{
    h_true_ixnsAboveKEThresholdPerSpill->Fill(n);
}


// Fill Muon kinematics histograms
void TruthHists::FillMuonKinematics(double cosL, double Elep, bool Numubar, bool passesTruthMx2)
{
    h_true_CosL_zoomOut->Fill(cosL);
    if(Numubar)
        h_true_CosLNumubar_zoomOut->Fill(cosL);
    if(!Numubar)
        h_true_CosLNumu_zoomOut->Fill(cosL);

    if (passesTruthMx2) {
        h_true_CosL->Fill(cosL);
        h_true_Elep->Fill(Elep);
    }
}

// Fill histogram for neutrino energy
void TruthHists::FillEnu(double Enu)
{
    h_true_Enu->Fill(Enu);
}

// Fill histograms related to pi0 and shower multiplicity
void TruthHists::FillPrimPi0Multiplicity(int nPrimPi0)
{
    h_true_nPrimPi0->Fill(nPrimPi0);
}

void TruthHists::FillSecPi0MultiplicityPreMx2(int nSecPi0)
{
    h_true_nSecPi0_preMx2->Fill(nSecPi0);
}

void TruthHists::FillSecPi0MultiplicityPostMx2(int nSecPi0)
{
    h_true_nSecPi0_postMx2->Fill(nSecPi0);
}

void TruthHists::FillSecPi0MultiplicityPreMx2(int nSecPi0)
{
    h_true_nSecPi0_preMx2->Fill(nSecPi0);
}

void TruthHists::FillShowerMultiplicityPreMx2(int nPrimElectron, int nSecElectron, 
                                              int nPrimPhoton, int nSecPhoton)
{
    h_true_nPrimElectron_preMx2->Fill(nPrimElectron);
    h_true_nSecElectron_preMx2->Fill(nSecElectron);
    h_true_nPrimPhoton_preMx2->Fill(nPrimPhoton);
    h_true_nSecPhoton_preMx2->Fill(nSecPhoton);
    int total_prim_showers = nPrimElectron + nPrimPhoton;
    int total_sec_showers = nSecElectron + nSecPhoton;
    h_true_nPrimShower_preMx2->Fill(total_prim_showers);
    h_true_nSecShower_preMx2->Fill(total_sec_showers);
}

void TruthHists::FillShowerMultiplicityPostMx2(int nPrimElectron, int nSecElectron, 
                                              int nPrimPhoton, int nSecPhoton)
{
    h_true_nPrimElectron_postMx2->Fill(nPrimElectron);
    h_true_nSecElectron_postMx2->Fill(nSecElectron);
    h_true_nPrimPhoton_postMx2->Fill(nPrimPhoton);
    h_true_nSecPhoton_postMx2->Fill(nSecPhoton);
    int total_prim_showers = nPrimElectron + nPrimPhoton;
    int total_sec_showers = nSecElectron + nSecPhoton;
    h_true_nPrimShower_postMx2->Fill(total_prim_showers);
    h_true_nSecShower_postMx2->Fill(total_sec_showers);
}

// Fill histogram for ixn mode 
void TruthHists::FillIxnMode(int mode)
{
    h_true_IxnMode->Fill(mode);
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

    h_true_nPrimPi0->Write();
    h_true_nSecPi0_preMx2->Write();
    h_true_nSecPi0_postMx2->Write();

    h_true_nPrimElectron_preMx2->Write();
    h_true_nPrimElectron_postMx2->Write();
    h_true_nSecElectron_preMx2->Write();
    h_true_nSecElectron_postMx2->Write();

    h_true_nPrimPhoton_preMx2->Write();
    h_true_nPrimPhoton_postMx2->Write();
    h_true_nSecPhoton_preMx2->Write();
    h_true_nSecPhoton_postMx2->Write();

    h_true_nPrimShower_preMx2->Write();
    h_true_nPrimShower_postMx2->Write();
    h_true_nSecShower_preMx2->Write();
    h_true_nSecShower_postMx2->Write();

    h_true_IxnMode->Write();

}