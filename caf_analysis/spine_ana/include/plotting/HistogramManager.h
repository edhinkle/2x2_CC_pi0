#pragma once

#include "plotting/TruthHists.h"
#include "plotting/RecoHists.h"
#include "plotting/TruthMatchedHists.h"
#include "plotting/CutFlowManager.h"
#include "plotting/FluxSystHists.h"


class HistogramManager {
public:
    explicit HistogramManager(int nThrowsFlux, double nominalIntegratedFlux);

    void Write(TFile* file);

    TruthHists truth;
    RecoHists reco;
    TruthMatchedHists truthMatch;
    CutFlowManager cuts;
    FluxSystHists fluxSyst;

private:
    int fNThrowsFlux;
    double fNominalIntegratedFlux;

};