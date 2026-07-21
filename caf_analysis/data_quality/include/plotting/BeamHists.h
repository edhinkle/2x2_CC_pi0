#pragma once
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"

class BeamHists {
public:
    BeamHists(long NSpills);

    //--------------------------------------------------------
    // Filling histograms
    //--------------------------------------------------------
    void FillTotalSpillsProcessed(int n);
    void FillTotalPOT(double pot);
    void FillPOTperSpill(long spill, double pot);

    void FillPOTbyStartTime(double spillPOT, double spillStartTimeSec, 
                            double spillStartTimeNsec, double mx2spillStartTimeSec, 
                            double mx2spillStartTimeNsec);
    double pot, start_time_sec, start_time_nsec, mx2_start_time_sec, mx2_start_time_nsec;

    //--------------------------------------------------------
    // Write histograms to file
    //--------------------------------------------------------
    void Write(TDirectory* dir);

private:

    //--------------------------------------------------------
    // Histograms
    //--------------------------------------------------------  
    TH1D *h_beam_TotalSpillsProcessed;
    TH1D *h_beam_TotalPOT;
    TH1D *h_beam_POTperSpill;

    TTree *t_beam_SpillPOTAndStartTime;

    long fNSpills;
 
};