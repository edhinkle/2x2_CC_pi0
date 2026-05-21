#pragma once

#include "plotting/HistogramManager.h"

void HistogramManager::Write(TFile* f)
{
  f->cd();

  TDirectory* truthDir = file->mkdir("truth");
  truth.Write(truthDir);

  TDirectory* recoDir = file->mkdir("reco");
  reco.Write(recoDir);

  TDirectory* truthMatchDir = file->mkdir("truthMatched");
  truthMatch.Write(truthMatchDir);

  TDirectory* cutDir = file->mkdir("cutflows");
  cuts.Write(cutDir);
}