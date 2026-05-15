#pragma once
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"

class RecoHists {
public:
    RecoHists();

    //------------------------------------------
    // Filling histograms
    //------------------------------------------
    void FillRecoVertexNoCuts(double x, double y, double z);

    // Fill histograms for number of interactions above KE threshold per spill

    void Write(TDirectory* dir);

private:

    //------------------------------------------
    // Histograms
    //------------------------------------------  
    TH2D *h_reco_VertexXZNoCuts;  
    TH1D *h_reco_VertexXNoCuts;
    TH1D *h_reco_VertexYNoCuts;
    TH1D *h_reco_VertexZNoCuts;

    // 


};