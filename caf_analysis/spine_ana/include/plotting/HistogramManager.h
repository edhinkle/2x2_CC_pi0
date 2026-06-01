#pragma once

#include "plotting/TruthHists.h"
#include "plotting/RecoHists.h"
#include "plotting/TruthMatchedHists.h"
#include "plotting/CutFlowManager.h"


class HistogramManager {
public:
    HistogramManager()
        : truth(), reco(), truthMatch(), cuts() {}

    void Write(TFile* file);

    TruthHists truth;
    RecoHists reco;
    TruthMatchedHists truthMatch;
    CutFlowManager cuts;

};