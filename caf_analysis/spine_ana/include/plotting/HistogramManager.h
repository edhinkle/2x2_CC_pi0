#pragma once

#include "plotting/TruthHists.h"
#include "plotting/RecoHists.h"
#include "plotting/MatchedHists.h"
#include "plotting/CutFlowManager.h"


class HistogramManager {
public:
    HistogramManager()
        : truth(), reco(), match(), cuts() {}

    void Write(TFile* f);

    TruthHists truth;
    RecoHists reco;
    MatchedHists match;
    CutFlowManager cuts;
};