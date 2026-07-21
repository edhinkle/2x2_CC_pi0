#include "plotting/BeamHists.h"
#include "TDirectory.h"
#include <iostream>

// Define all beam histograms 
BeamHists::BeamHists(long NSpills)
        : fNSpills(NSpills)

{
    // Count spills processed
    h_beam_TotalSpillsProcessed = new TH1D("h_beam_TotalSpillsProcessed", "h_beam_TotalSpillsProcessed", 1, 0, 1);

    // Count total POT processed
    h_beam_TotalPOT = new TH1D("h_beam_TotalPOT", "h_beam_TotalPOT", 1, 0, 1);

    // Count POT per spill
    h_beam_POTperSpill = new TH1D("h_beam_POTperSpill", "h_beam_POTperSpill", fNSpills, 0, fNSpills);

    // Count POT by start time
    t_beam_SpillPOTAndStartTime = new TTree("t_beam_SpillPOTAndStartTime", "t_beam_SpillPOTAndStartTime");
    t_beam_SpillPOTAndStartTime->Branch("pot", &pot, "pot/D");
    t_beam_SpillPOTAndStartTime->Branch("start_time_sec", &start_time_sec, "start_time_sec/D");
    t_beam_SpillPOTAndStartTime->Branch("start_time_nsec", &start_time_nsec, "start_time_nsec/D");
    t_beam_SpillPOTAndStartTime->Branch("mx2_start_time_sec", &mx2_start_time_sec, "mx2_start_time_sec/D");
    t_beam_SpillPOTAndStartTime->Branch("mx2_start_time_nsec", &mx2_start_time_nsec, "mx2_start_time_nsec/D");

}


//------------------------------------------
// Histogram Filling Methods
//------------------------------------------
void BeamHists::FillTotalSpillsProcessed(int n)
{
    h_beam_TotalSpillsProcessed->Fill(0.5,n);
}

void BeamHists::FillTotalPOT(double pot)
{
    h_beam_TotalPOT->Fill(0.5, pot);
}

void BeamHists::FillPOTperSpill(long spill, double pot)
{
    h_beam_POTperSpill->Fill(spill, pot);
}

void BeamHists::FillPOTbyStartTime(double spillPOT, double spillStartTimeSec, 
                                   double spillStartTimeNsec, double mx2spillStartTimeSec, 
                                   double mx2spillStartTimeNsec)
{
    pot = spillPOT;
    start_time_sec = spillStartTimeSec;
    start_time_nsec = spillStartTimeNsec;
    mx2_start_time_sec = mx2spillStartTimeSec;
    mx2_start_time_nsec = mx2spillStartTimeNsec;
    //std::cout << "Filling POT by start time: " << pot << ", " << start_time_sec << ", " << start_time_nsec << std::endl;
    t_beam_SpillPOTAndStartTime->Fill();
}


//------------------------------------------
// Histogram Writing Method
//------------------------------------------  
void BeamHists::Write(TDirectory* dir)
{
    dir->cd();

    h_beam_TotalSpillsProcessed->Write();
    h_beam_TotalPOT->Write();
    h_beam_POTperSpill->Write();

    t_beam_SpillPOTAndStartTime->Write();

}