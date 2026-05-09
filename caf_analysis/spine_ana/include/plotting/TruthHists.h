#pragma once
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"

class TruthHists {
public:
    TruthHists(double flux_nom);

    //------------------------------------------
    // Filling histograms
    //------------------------------------------

    // Fill histograms for number of interactions above KE threshold per spill
    void FillInteractionsAboveKEThresholdPerSpill(int n);

    void FillMuonKinematics(double cosL, double Elep, double Enu);
    void FillMultiplicity(int nTrk, int nShort, int nLong);
    void FillPionInfo(int nPi, int escapingPi);

    void Write(TDirectory* dir);

private:

    //------------------------------------------
    // Histograms
    //------------------------------------------    

    TH1D* h_true_ixnsAboveKEThresholdPerSpill;




    TH1D* h_true_Enu;

    TH1D* h_true_partMult;
    TH1D* h_true_showerMult;
    TH1D* h_true_ixn


    TH1D* h_trueCosL;
    TH1D* h_trueCosL_zoomOut;
    TH1D* h_trueEnu;

    TH1D* h_true_mult;
    TH1D* h_true_multShort;
    TH1D* h_true_multLong;

    TH1D* h_nPi;
    TH1D* h_escapePi;
};