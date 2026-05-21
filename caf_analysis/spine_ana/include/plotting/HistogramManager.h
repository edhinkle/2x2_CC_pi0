#pragma once

#include "plotting/TruthHists.h"
#include "plotting/RecoHists.h"
#include "plotting/MatchedHists.h"
#include "plotting/CutFlowManager.h"


class HistogramManager {
public:
    HistogramManager()
        : truth(), reco(), truthMatch(), cuts() {}

    void Write(TFile* f);

    TruthHists truth;
    RecoHists reco;
    TruthMatchedHists truthMatch;
    CutFlowManager cuts;

};