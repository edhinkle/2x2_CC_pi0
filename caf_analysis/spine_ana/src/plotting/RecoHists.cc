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


void RecoHists::Write(TDirectory* dir)
{
    dir->cd();

    h_reco_VertexXZNoCuts->Write();
    h_reco_VertexXNoCuts->Write();
    h_reco_VertexYNoCuts->Write();
    h_reco_VertexZNoCuts->Write();

}