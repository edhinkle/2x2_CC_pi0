#pragma once
#include "config/ConfigLoader.h"
#include "io/CAFUtils.h"
#include "analysis/TruthInteractionSummary.h"
#include "plotting/HistogramManager.h"
#include "cuts/DetectorCuts.h"
#include "systematics/flux/FluxSystManager.h"
#include "duneanaobj/StandardRecord/StandardRecord.h" 

class TruthSelection {
public:

  TruthSelection(const SelectionConfig& cfg,
                 const BeamConfig& beam,
                 const DetectorConfig& detector, 
                 const FluxSystManager& fluxSyst):
                fSelCuts(cfg), fBeam(beam), fDetector(detector), fFluxSyst(fluxSyst) {};

  void SelectTruthInteractions(const caf::StandardRecord& sr,
                               HistogramManager& hist);

  bool IxnPassesTruthLArCuts(const TruthInteractionSummary& summary) const;

  TruthInteractionSummary BuildTruthSummary(
      const caf::SRTrueInteraction& nu) const;

private:

  const SelectionConfig& fSelCuts;
  const BeamConfig& fBeam;
  const DetectorConfig& fDetector;
  const FluxSystManager& fFluxSyst;

};