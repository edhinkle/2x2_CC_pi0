// plotting/FluxSystHists.cc
#include "plotting/FluxSystHists.h"
#include "TFile.h"
#include <stdexcept>
#include <iostream>

// =====================================================
// Initialize
// =====================================================
void FluxSystHists::Initialize(
    const TH1D* nomTrueCosL,
    const TH2D* nomRespCosL,
    const TH1D* nomRecoCosL,
    const TH1D* nomTruthMatchCosL,
    int nThrows)
{
    if (!nomTrueCosL || !nomRespCosL || 
        !nomRecoCosL || !nomTruthMatchCosL)
        throw std::invalid_argument("FluxSystHists::Initialize: null nominal histogram passed");

    fNThrows = nThrows;

    fTrueCosL.reserve(nThrows);
    fRespCosL.reserve(nThrows);
    fRecoCosL.reserve(nThrows);
    fTruthMatchCosL.reserve(nThrows);

    for (int i = 0; i < nThrows; ++i) {
        // Clone each nominal histogram with a unique name, then push by value.
        // The Clone returns a heap pointer; we dereference into the vector and
        // let ROOT's copy constructor handle the internal state. The temporary
        // pointer is deleted immediately after.
        auto* hTrueCosL   = (TH1D*)nomTrueCosL    ->Clone(Form("h_fluxSyst_true_CosL_%d",   i));
        auto* hRespCosL   = (TH2D*)nomRespCosL->Clone(Form("h_fluxSyst_truthMatchIxn_Reco_responseMuonCosL_%d", i));
        auto* hRecoCosL   = (TH1D*)nomRecoCosL->Clone(Form("h_fluxSyst_reco_CosL_%d",   i));
        auto* hTruthMatchCosL   = (TH1D*)nomTruthMatchCosL->Clone(Form("h_fluxSyst_truthMatchIxn_MuonCosL_%d", i));

        // Detach from any directory so ROOT doesn't double-delete on file close
        hTrueCosL ->SetDirectory(nullptr);
        hRespCosL ->SetDirectory(nullptr);
        hRecoCosL ->SetDirectory(nullptr);
        hTruthMatchCosL ->SetDirectory(nullptr);

        fTrueCosL.push_back(std::move(*hTrueCosL));
        fRespCosL.push_back(std::move(*hRespCosL));
        fRecoCosL.push_back(std::move(*hRecoCosL));
        fTruthMatchCosL.push_back(std::move(*hTruthMatchCosL));

        // Heap copies are now hollow shells — delete them
        delete hTrueCosL;
        delete hRespCosL;
        delete hRecoCosL;
        delete hTruthMatchCosL;
    }
}

// =====================================================
// Per-throw fill methods
// =====================================================
void FluxSystHists::FillTrue(int throwIdx, double cosL, double w)
{
    if (throwIdx < 0 || throwIdx >= fNThrows) return;
    fTrueCosL[throwIdx].Fill(cosL,  w);
}

void FluxSystHists::FillReco(int throwIdx, double cosL, double w)
{
    if (throwIdx < 0 || throwIdx >= fNThrows) return;
    fRecoCosL[throwIdx].Fill(cosL, w);
}

void FluxSystHists::FillResponse(int throwIdx,
                                      double recoCosL, double trueCosL,
                                      double w)
{
    if (throwIdx < 0 || throwIdx >= fNThrows) return;
    fRespCosL[throwIdx].Fill(recoCosL, trueCosL, w);
}

void FluxSystHists::FillTruthMatch(int throwIdx, double cosL, double w)
{
    if (throwIdx < 0 || throwIdx >= fNThrows) return;
    fTruthMatchCosL[throwIdx].Fill(cosL, w);
}


// =====================================================
// SetThrowFluxTitles
// =====================================================
void FluxSystHists::SetThrowFluxTitles(const std::vector<double>& throwFluxes)
{
    if ((int)throwFluxes.size() != fNThrows) {
        std::cerr << "FluxSystHists::SetThrowFluxTitles: "
                  << "throwFluxes size (" << throwFluxes.size()
                  << ") does not match fNThrows (" << fNThrows << ")\n";
        return;
    }

    for (int i = 0; i < fNThrows; ++i) {
        const char* title = Form("%0.08f", throwFluxes[i]);
        // Only recoCosL carries the flux title in the original script;
        // set it on all reco histograms so downstream unfolding code
        // has a consistent place to read the normalization from.
        fRecoCosL[i].SetTitle(title);
    }
}


void FluxSystHists::Finalize()
{
    // Profiles are built from the binning of the throw histograms,
    // so they can only be constructed once Initialize() has been called
    fTrueCosLProfile = TProfile("h_trueCosLProfile", "h_trueCosLProfile",
                                fTrueCosL[0].GetNbinsX(),
                                fTrueCosL[0].GetXaxis()->GetXbins()->GetArray(),
                                "S");


    for (int i = 0; i < fNThrows; ++i) {
        for (int bin = 1; bin <= fTrueCosL[i].GetNbinsX(); ++bin) {
            double binCenter = fTrueCosL[i].GetXaxis()->GetBinCenter(bin);
            fTrueCosLProfile.Fill(binCenter, fTrueCosL[i].GetBinContent(bin));
        }
    }
}


// =====================================================
// Write
// =====================================================
void FluxSystHists::Write(TDirectory* dir) const
{

    dir->cd();

    for (int i = 0; i < fNThrows; ++i) {
        fTrueCosL[i].Write();
        fRespCosL[i].Write();
        fRecoCosL[i].Write();
        fTruthMatchCosL[i].Write();
    }
    fTrueCosLProfile.Write();
}