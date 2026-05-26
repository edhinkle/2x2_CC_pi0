#pragma once

#include "plotting/RecoHists.h"
#include "TDirectory.h"

// Define all reco histograms 
RecoHists::RecoHists()
{
    
    // All Reco interaction vertices (no cuts)
    h_reco_VertexXZNoCuts = new TH2D("h_reco_VertexXZNoCuts", "h_reco_VertexXZNoCuts", 70, -70, 70, 70, -70, 70);
    h_reco_VertexXNoCuts = new TH1D("h_reco_VertexXNoCuts", "h_reco_VertexXNoCuts", 140, -70, 70);
    h_reco_VertexYNoCuts = new TH1D("h_reco_VertexYNoCuts", "h_reco_VertexYNoCuts", 140, -70, 70);
    h_reco_VertexZNoCuts = new TH1D("h_reco_VertexZNoCuts", "h_reco_VertexZNoCuts", 140, -70, 70);

    // All Reco interaction vertices (with cuts)
    h_reco_VertexXZWithCuts = new TH2D("h_reco_VertexXZWithCuts", "h_reco_VertexXZWithCuts", 70, -70, 70, 70, -70, 70);
    h_reco_VertexXWithCuts = new TH1D("h_reco_VertexXWithCuts", "h_reco_VertexXWithCuts", 140, -70, 70);
    h_reco_VertexYWithCuts = new TH1D("h_reco_VertexYWithCuts", "h_reco_VertexYWithCuts", 140, -70, 70);
    h_reco_VertexZWithCuts = new TH1D("h_reco_VertexZWithCuts", "h_reco_VertexZWithCuts", 140, -70, 70);


    // Reco cos(lepton) wrt beam
    // Edges for CosL histograms
    Double_t edges[7] = {0.91, 0.96, 0.98, 0.9887, 0.994, 0.9974, 1};

    // Muon kinematics
    h_reco_CosL = new TH1D("h_reco_CosL", "h_reco_CosL", 6, edges); // Signal region (passes Mx2 cuts)
    h_reco_CosL_zoomOut = new TH1D("h_reco_CosL_zoomOut", "h_reco_CosL_zoomOut", 50, 0.8, 1); // All events passing initial cuts
 
    // Shower multiplicity
    // Bins and edges for multiplicity histograms
    int binsMult = 20;
    Double_t edgesMult[21] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
    h_reco_PrimShowerMultiplicity = new TH1D("h_reco_PrimShowerMultiplicity", "h_reco_PrimShowerMultiplicity", binsMult, edgesMult);
    h_reco_SecShowerMultiplicity = new TH1D("h_reco_SecShowerMultiplicity", "h_reco_SecShowerMultiplicity", binsMult, edgesMult);
    h_reco_PrimElectronMultiplicity = new TH1D("h_reco_PrimElectronMultiplicity", "h_reco_PrimElectronMultiplicity", binsMult, edgesMult);
    h_reco_SecElectronMultiplicity = new TH1D("h_reco_SecElectronMultiplicity", "h_reco_SecElectronMultiplicity", binsMult, edgesMult);
    h_reco_PrimPhotonMultiplicity = new TH1D("h_reco_PrimPhotonMultiplicity", "h_reco_PrimPhotonMultiplicity", binsMult, edgesMult);
    h_reco_SecPhotonMultiplicity = new TH1D("h_reco_SecPhotonMultiplicity", "h_reco_SecPhotonMultiplicity", binsMult, edgesMult);

}


//------------------------------------------
// Histogram Filling Methods
//------------------------------------------    
void RecoHists::FillRecoVertexNoCuts(double x, double y, double z)
{
    h_reco_VertexXZNoCuts->Fill(x, z);
    h_reco_VertexXNoCuts->Fill(x);
    h_reco_VertexYNoCuts->Fill(y);
    h_reco_VertexZNoCuts->Fill(z);
}

void FillRecoVertexWithCuts(caf::SRVector3D vertex) 
{
    h_reco_VertexXZWithCuts->Fill(vertex.x, vertex.z);
    h_reco_VertexXWithCuts->Fill(vertex.x);
    h_reco_VertexYWithCuts->Fill(vertex.y);
    h_reco_VertexZWithCuts->Fill(vertex.z);
}
  
void RecoHists::FillRecoCosMuonAngle(double cosL)
{
    h_reco_CosL_zoomOut->Fill(cosL);
    h_reco_CosL->Fill(cosL);
}

void RecoHists::FillRecoShowerMultiplicity(RecoInteractionSummary& recoSummary)
{

    int total_prim_showers = recoSummary.nPrimElectron + recoSummary.nPrimPhoton;
    int total_sec_showers = recoSummary.nSecElectron + recoSummary.nSecPhoton;
    h_reco_PrimShowerMultiplicity->Fill(total_prim_showers);
    h_reco_SecShowerMultiplicity->Fill(total_sec_showers);

    h_reco_PrimElectronMultiplicity->Fill(recoSummary.nPrimElectron);
    h_reco_SecElectronMultiplicity->Fill(recoSummary.nSecElectron);

    h_reco_PrimPhotonMultiplicity->Fill(recoSummary.nPrimPhoton);
    h_reco_SecPhotonMultiplicity->Fill(recoSummary.nSecPhoton);
}

//------------------------------------------
// Histogram Writing Method
//------------------------------------------  
void RecoHists::Write(TDirectory* dir)
{
    dir->cd();

    h_reco_VertexXZNoCuts->Write();
    h_reco_VertexXNoCuts->Write();
    h_reco_VertexYNoCuts->Write();
    h_reco_VertexZNoCuts->Write();

    h_reco_VertexXZWithCuts->Write();
    h_reco_VertexXWithCuts->Write();
    h_reco_VertexYWithCuts->Write();
    h_reco_VertexZWithCuts->Write();

    h_reco_CosL_zoomOut->Write();
    h_reco_CosL->Write();

    h_reco_PrimShowerMultiplicity->Write();
    h_reco_SecShowerMultiplicity->Write();
    h_reco_PrimElectronMultiplicity->Write();
    h_reco_SecElectronMultiplicity->Write();
    h_reco_PrimPhotonMultiplicity->Write();
    h_reco_SecPhotonMultiplicity->Write();


}