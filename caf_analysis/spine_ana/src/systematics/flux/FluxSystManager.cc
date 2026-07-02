// systematics/flux/FluxSystManager.cc
#include "systematics/flux/FluxSystManager.h"
#include "TFile.h"
#include "TH1D.h"
#include <stdexcept>

FluxSystManager::FluxSystManager(const FluxSystConfig& cfg)
    : fFluxSyst(cfg), fRW(fFluxSyst), fNThrows(cfg.nThrows), fEnabled(cfg.enableFluxSyst), 
      fIsRHC(cfg.isRHC), fNuSignMatters(cfg.nuSignMatters), fSignalNuPDG(cfg.signalNuPDG)
{
    if (!fEnabled) {
        fThrowFluxes.assign(fNThrows, 1.0);
        fIntegratedFluxNominal = fRW.GetIntegratedFluxNominal();
        return;
    }

    fRW.GenerateThrows(fNThrows, fFluxSyst.seed, fIsRHC, fSignalNuPDG, fNuSignMatters); // internally loads binning, covariance, flux histograms

    fIntegratedFluxNominal = fRW.GetIntegratedFluxNominal();

    if (!fEnabled) {
        fThrowFluxes.assign(fNThrows, 1.0);
        return;
    }
    else {
        fThrowFluxes = fRW.GetIntegratedFluxThrows(); // returns std::vector<double> directly
    }

    fBinRange = fRW.GetBinRangeForBeamModeAndPDG(fIsRHC, fSignalNuPDG, fNuSignMatters);
}

std::vector<double> FluxSystManager::GetWeights(int nuPDG, double nuE, bool isRHC) const
{
    std::vector<double> weights(fNThrows, 1.0);
    if (!fEnabled) return weights;

    int binIdx = fRW.GetBinIndex(nuPDG, nuE, isRHC);
    if (binIdx < 0) return weights; // out of range → all nominal

    // FillAllThrowsForBin already initializes to 1.0 for invalid binIdx,
    // so this is safe. Use a raw array to match the API.
    std::vector<Double_t> raw(fNThrows);
    if (binIdx < fBinRange.first || binIdx >= fBinRange.second) {
        // Outside the range of interest for this beam mode and PDG
        return weights; // all nominal
    }
    fRW.FillAllThrowsForBin(binIdx, raw.data());

    // Clamp unphysical weights ... not sure why they are unphysical?? -- removing logic for now
    for (int i = 0; i < fNThrows; ++i) {
        weights[i] = raw[i]; //(raw[i] > 0.01 && raw[i] < 20.0) ? raw[i] : 1.0;
    }

    return weights;
}