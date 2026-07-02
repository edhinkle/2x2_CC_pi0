// systematics/flux/FluxSystManager.h
#pragma once
#include "systematics/flux/FluxCovarianceReweight.h"
#include "config/ConfigLoader.h"
#include <vector>
#include <string>

class FluxSystManager {
public:
    explicit FluxSystManager(const FluxSystConfig& cfg);

    // Returns nThrows weights for the neutrino at this PDG+energy.
    // Returns a vector of 1.0s if flux syst is disabled or bin not found.
    std::vector<double> GetWeights(int nuPDG, double nuE, bool isRHC) const;

    // Nominal integrated flux (for reference normalization)
    double GetNominalIntegratedFlux() const { return fIntegratedFluxNominal; }

    // Per-throw integrated flux values (length = nThrows)
    const std::vector<double>& GetThrowIntegratedFluxes() const { return fThrowFluxes; }

    // Get integrated flux throws from FluxCovarianceReweight
    const std::vector<double>& GetIntegratedFluxThrows() const {
        return fRW.GetIntegratedFluxThrows();
    }

    int GetNThrows() const { return fNThrows; }
    bool IsEnabled() const { return fEnabled; }
    bool IsRHC() const { return fIsRHC; }
    bool NuSignMatters() const { return fNuSignMatters; }
    int GetSignalNuPDG() const { return fSignalNuPDG; }

private:
    FluxSystConfig fFluxSyst;
    FluxCovarianceReweight fRW;
    double fIntegratedFluxNominal;
    std::vector<double> fThrowFluxes;
    std::pair<int, int> fBinRange; // bin range for the given beam mode and PDG (start, end)
    int fNThrows;
    bool fEnabled;
    bool fIsRHC; // whether to use RHC or FHC flux
    bool fNuSignMatters; // whether to distinguish neutrinos vs antineutrinos (e.g. for RHC flux)
    int fSignalNuPDG; // PDG code of the neutrino flavor to apply flux systematics to (e.g. 14 for numu)
};