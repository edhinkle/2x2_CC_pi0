#pragma once
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "analysis/RecoInteractionSummary.h"

class RecoHists {
public:
    RecoHists(double nominalIntegratedFlux);

    //--------------------------------------------------------
    // Filling histograms
    //--------------------------------------------------------
    void FillTotalSpillsProcessed(int n);
    void FillTotalPOT(double pot);
    void FillRecoVertexNoCuts(caf::SRVector3D vertex);
    void FillRecoVertexWithCuts(caf::SRVector3D vertex);
    void FillRecoCosMuonAngle(double cosL);
    void FillRecoShowerMultiplicity(RecoInteractionSummary& recoSummary);

    //--------------------------------------------------------
    // Get methods for copying histograms for systematics
    //--------------------------------------------------------
    const TH1D* GetRecoCosL() const { return h_reco_CosL_zoomOut; }

    //--------------------------------------------------------
    // Write histograms to file
    //--------------------------------------------------------
    void Write(TDirectory* dir);

private:

    //--------------------------------------------------------
    // Histograms
    //--------------------------------------------------------  
    TH1D *h_reco_TotalSpillsProcessed;
    TH1D *h_reco_TotalPOT;

    TH2D *h_reco_VertexZXNoCuts;  
    TH1D *h_reco_VertexXNoCuts;
    TH1D *h_reco_VertexYNoCuts;
    TH1D *h_reco_VertexZNoCuts;

    TH2D *h_reco_VertexZXWithCuts;  
    TH1D *h_reco_VertexXWithCuts;
    TH1D *h_reco_VertexYWithCuts;
    TH1D *h_reco_VertexZWithCuts;

    TH1D *h_reco_CosL;
    TH1D *h_reco_CosL_zoomOut;

    TH1D *h_reco_PrimShowerMultiplicity;
    TH1D *h_reco_SecShowerMultiplicity;
    TH1D *h_reco_PrimElectronMultiplicity;
    TH1D *h_reco_SecElectronMultiplicity;
    TH1D *h_reco_PrimPhotonMultiplicity;
    TH1D *h_reco_SecPhotonMultiplicity;

    double fNominalIntegratedFlux;

};