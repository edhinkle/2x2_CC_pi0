#pragma once
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"

class TruthHists {
public:
    TruthHists();

    //------------------------------------------
    // Filling histograms
    //------------------------------------------

    // Fill histograms for number of interactions above KE threshold per spill
    void FillInteractionsAboveKEThresholdPerSpill(int n);

    void FillMuonKinematics(double cosL, double Elep, bool Numubar);
    void FillEnu(double Enu);
    void FillPi0Multiplicity(int nPi0);

    void Write(TDirectory* dir);

private:

    //------------------------------------------
    // Histograms
    //------------------------------------------    

    TH1D* h_true_ixnsAboveKEThresholdPerSpill;

    TH1D* h_true_CosL;
    TH1D* h_true_CosL_zoomOut;
    TH1D* h_true_CosLNumubar_zoomOut;
    TH1D* h_true_CosLNumu_zoomOut;
    TH1D* h_true_Elep;

    TH1D* h_true_Enu;

    TH1D* h_true_nPi0;

};