#include "plotting/RecoHists.h"
#include "TDirectory.h"

// Define all reco histograms 
RecoHists::RecoHists(double nominalIntegratedFlux)
    : fNominalIntegratedFlux(nominalIntegratedFlux)
{
    // Count spills processed
    h_reco_TotalSpillsProcessed = new TH1D("h_reco_TotalSpillsProcessed", "h_reco_TotalSpillsProcessed", 1, 0, 1);

    // Count total POT processed
    h_reco_TotalPOT = new TH1D("h_reco_TotalPOT", "h_reco_TotalPOT", 1, 0, 1);
        
    // All Reco interaction vertices (no cuts)
    h_reco_VertexZXNoCuts = new TH2D("h_reco_VertexZXNoCuts", "h_reco_VertexZXNoCuts", 70, -70, 70, 70, -70, 70);
    h_reco_VertexXNoCuts = new TH1D("h_reco_VertexXNoCuts", "h_reco_VertexXNoCuts", 140, -70, 70);
    h_reco_VertexYNoCuts = new TH1D("h_reco_VertexYNoCuts", "h_reco_VertexYNoCuts", 140, -70, 70);
    h_reco_VertexZNoCuts = new TH1D("h_reco_VertexZNoCuts", "h_reco_VertexZNoCuts", 140, -70, 70);

    // All Reco interaction vertices (with cuts)
    h_reco_VertexZXWithCuts = new TH2D("h_reco_VertexZXWithCuts", "h_reco_VertexZXWithCuts", 70, -70, 70, 70, -70, 70);
    h_reco_VertexXWithCuts = new TH1D("h_reco_VertexXWithCuts", "h_reco_VertexXWithCuts", 140, -70, 70);
    h_reco_VertexYWithCuts = new TH1D("h_reco_VertexYWithCuts", "h_reco_VertexYWithCuts", 140, -70, 70);
    h_reco_VertexZWithCuts = new TH1D("h_reco_VertexZWithCuts", "h_reco_VertexZWithCuts", 140, -70, 70);


    // Reco cos(lepton) wrt beam
    // Edges for CosL histograms
    Double_t edges[7] = {0.90, 0.96, 0.98, 0.9887, 0.994, 0.9974, 1};

    // Muon kinematics
    h_reco_CosL = new TH1D("h_reco_CosL", "h_reco_CosL", 6, edges); // Signal region (passes Mx2 cuts)
    h_reco_CosL_zoomOut = new TH1D("h_reco_CosL_zoomOut", "h_reco_CosL_zoomOut", 50, 0.8, 1); // All events passing initial cuts
    h_reco_CosL->SetTitle(Form("%0.08f", fNominalIntegratedFlux));
    h_reco_CosL_zoomOut->SetTitle(Form("%0.08f", fNominalIntegratedFlux));


    // Shower multiplicity
    // Bins and edges for multiplicity histograms
    //int binsMult = 20;
    //Double_t edgesMult[21] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
    h_reco_PrimShowerMultiplicity = new TH1D("h_reco_PrimShowerMultiplicity", "h_reco_PrimShowerMultiplicity", 20, 0, 20);
    h_reco_SecShowerMultiplicity = new TH1D("h_reco_SecShowerMultiplicity", "h_reco_SecShowerMultiplicity", 20, 0, 20);
    h_reco_PrimElectronMultiplicity = new TH1D("h_reco_PrimElectronMultiplicity", "h_reco_PrimElectronMultiplicity", 20, 0, 20);
    h_reco_SecElectronMultiplicity = new TH1D("h_reco_SecElectronMultiplicity", "h_reco_SecElectronMultiplicity", 20, 0, 20);
    h_reco_PrimPhotonMultiplicity = new TH1D("h_reco_PrimPhotonMultiplicity", "h_reco_PrimPhotonMultiplicity", 20, 0, 20);
    h_reco_SecPhotonMultiplicity = new TH1D("h_reco_SecPhotonMultiplicity", "h_reco_SecPhotonMultiplicity", 20, 0, 20);

}


//------------------------------------------
// Histogram Filling Methods
//------------------------------------------
void RecoHists::FillTotalSpillsProcessed(int n)
{
    h_reco_TotalSpillsProcessed->Fill(0.5,n);
}

void RecoHists::FillTotalPOT(double pot)
{
    h_reco_TotalPOT->Fill(0.5, pot);
}

void RecoHists::FillRecoVertexNoCuts(caf::SRVector3D vertex)
{
    h_reco_VertexZXNoCuts->Fill(vertex.z, vertex.x);
    h_reco_VertexXNoCuts->Fill(vertex.x);
    h_reco_VertexYNoCuts->Fill(vertex.y);
    h_reco_VertexZNoCuts->Fill(vertex.z);
}

void RecoHists::FillRecoVertexWithCuts(caf::SRVector3D vertex) 
{
    h_reco_VertexZXWithCuts->Fill(vertex.z, vertex.x);
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

    int total_prim_showers = recoSummary.nPrimElectrons + recoSummary.nPrimPhotons;
    int total_sec_showers = recoSummary.nSecElectrons + recoSummary.nSecPhotons;
    h_reco_PrimShowerMultiplicity->Fill(total_prim_showers);
    h_reco_SecShowerMultiplicity->Fill(total_sec_showers);

    h_reco_PrimElectronMultiplicity->Fill(recoSummary.nPrimElectrons);
    h_reco_SecElectronMultiplicity->Fill(recoSummary.nSecElectrons);

    h_reco_PrimPhotonMultiplicity->Fill(recoSummary.nPrimPhotons);
    h_reco_SecPhotonMultiplicity->Fill(recoSummary.nSecPhotons);

    // DEBUG
    if (total_prim_showers > 20) {
        std::cout << "Reco primary showers: " << total_prim_showers
                  << " (nPrimElectrons: " << recoSummary.nPrimElectrons 
                  << ", nPrimPhotons: " << recoSummary.nPrimPhotons << ")" << std::endl;
    }
    if (recoSummary.nPrimElectrons > 20) std::cout << "Reco primary electrons: " << recoSummary.nPrimElectrons << std::endl;
    if (recoSummary.nPrimPhotons > 20) std::cout << "Reco primary photons: " << recoSummary.nPrimPhotons << std::endl;

}

//------------------------------------------
// Histogram Writing Method
//------------------------------------------  
void RecoHists::Write(TDirectory* dir)
{
    dir->cd();

    h_reco_TotalSpillsProcessed->Write();
    h_reco_TotalPOT->Write();

    h_reco_VertexZXNoCuts->Write();
    h_reco_VertexXNoCuts->Write();
    h_reco_VertexYNoCuts->Write();
    h_reco_VertexZNoCuts->Write();

    h_reco_VertexZXWithCuts->Write();
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