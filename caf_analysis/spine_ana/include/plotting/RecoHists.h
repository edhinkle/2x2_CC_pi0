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
    void FillRecoVertexWithCuts(caf::SRVector3D vertex);
    void FillRecoCosMuonAngle(double cosL)

    void Write(TDirectory* dir);

private:

    //------------------------------------------
    // Histograms
    //------------------------------------------  
    TH2D *h_reco_VertexXZNoCuts;  
    TH1D *h_reco_VertexXNoCuts;
    TH1D *h_reco_VertexYNoCuts;
    TH1D *h_reco_VertexZNoCuts;

    TH1D *h_reco_CosL;
    TH1D *h_reco_CosL_zoomOut;

};