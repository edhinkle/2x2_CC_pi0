// plotting/FluxSystHists.h
#pragma once
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include <vector>

// Owns the nThrows-deep vectors of cloned histograms.
// Initialized from the nominal histograms already in HistogramManager.
class FluxSystHists {
public:
    void Initialize(const TH1D* nomTrueCosL,
                    const TH2D* nomRespCosL,
                    const TH1D* nomRecoCosL,
                    const TH1D* nomTruthMatchCosL,
                    int nThrows);

    // Per-throw fill methods
    void FillTrue(int throwIdx, double cosL, double w);
    void FillReco(int throwIdx, double cosL, double w);
    void FillResponse(int throwIdx, double recoCosL, double trueCosL, double w);
    void FillTruthMatch(int throwIdx, double cosL, double w);

    // Set the per-throw flux value as the histogram title (matches original logic)
    void SetThrowFluxTitles(const std::vector<double>& throwFluxes);

    void Finalize(); // call after spill loop to summarize histogram info

    void Write(TDirectory* dir) const;

private:
    int fNThrows = 0;
    std::vector<TH1D> fTrueCosL;
    std::vector<TH2D> fRespCosL;
    std::vector<TH1D> fRecoCosL;
    std::vector<TH1D> fTruthMatchCosL;
    TProfile fTrueCosLProfile;
};