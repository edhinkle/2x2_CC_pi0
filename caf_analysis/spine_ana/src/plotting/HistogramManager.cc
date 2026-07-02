#include "plotting/HistogramManager.h"


HistogramManager::HistogramManager(int nThrowsFlux, double nominalIntegratedFlux)
        : truth(nominalIntegratedFlux), reco(nominalIntegratedFlux), truthMatch(nominalIntegratedFlux), cuts(), fluxSyst(), fNThrowsFlux(nThrowsFlux), fNominalIntegratedFlux(nominalIntegratedFlux) 
{
  // Initialize flux syst histograms from the completed nominals
  fluxSyst.Initialize(
      truth.GetTrueCosL(),  // these are getter methods returning
      truthMatch.GetRespCosL(),  // const TH1D* to the nominal histograms
      reco.GetRecoCosL(),
      truthMatch.GetTruthMatchCosL(),
      nThrowsFlux
  );
}



void HistogramManager::Write(TFile* file)
{
  file->cd();

  TDirectory* truthDir = file->mkdir("truth");
  truth.Write(truthDir);

  TDirectory* recoDir = file->mkdir("reco");
  reco.Write(recoDir);

  TDirectory* truthMatchDir = file->mkdir("truthMatched");
  truthMatch.Write(truthMatchDir);

  TDirectory* cutDir = file->mkdir("cutflows");
  cuts.Write(cutDir);

  TDirectory* fluxSystDir = file->mkdir("fluxSyst");
  fluxSyst.Write(fluxSystDir);
}