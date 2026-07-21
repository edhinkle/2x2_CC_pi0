#include "plotting/HistogramManager.h"


HistogramManager::HistogramManager(long NSpills)
        : beam(NSpills), fNSpills(NSpills)
{
}

void HistogramManager::Write(TFile* file)
{
  file->cd();

  TDirectory* beamDir = file->mkdir("beam");
  beam.Write(beamDir);

}