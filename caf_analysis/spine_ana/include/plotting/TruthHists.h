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

    void FillMuonKinematics(double cosL, double Elep, bool Numubar, bool passesTruthMx2);
    void FillEnu(double Enu);
    void FillPrimPi0Multiplicity(int nPrimPi0);
    void FillSecPi0MultiplicityPreMx2(int nSecPi0);
    void FillSecPi0MultiplicityPostMx2(int nSecPi0);
    void FillShowerMultiplicityPreMx2(int nPrimElectron, int nSecElectron, 
                                      int nPrimPhoton, int nSecPhoton);
    void FillShowerMultiplicityPostMx2(int nPrimElectron, int nSecElectron, 
                                       int nPrimPhoton, int nSecPhoton);
    void FillIxnMode(int mode);
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

    TH1D* h_true_nPrimPi0;
    TH1D* h_true_nSecPi0_preMx2;
    TH1D* h_true_nSecPi0_postMx2;

    TH1D* h_true_nPrimElectron_preMx2;
    TH1D* h_true_nPrimElectron_postMx2;
    TH1D* h_true_nSecElectron_preMx2;
    TH1D* h_true_nSecElectron_postMx2;

    TH1D* h_true_nPrimPhoton_preMx2;
    TH1D* h_true_nPrimPhoton_postMx2;
    TH1D* h_true_nSecPhoton_preMx2;
    TH1D* h_true_nSecPhoton_postMx2;

    TH1D* h_true_nPrimShower_preMx2;
    TH1D* h_true_nPrimShower_postMx2;
    TH1D* h_true_nSecShower_preMx2;
    TH1D* h_true_nSecShower_postMx2;

    TH1D* h_true_IxnMode;

};