#pragma once

#include "plotting/BeamHists.h"



class HistogramManager {
public:
    explicit HistogramManager(long NSpills);

    void Write(TFile* file);

    BeamHists beam;

private:
    long fNSpills;
};