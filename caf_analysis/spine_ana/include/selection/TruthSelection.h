#pragma once
#include "config/ConfigLoader.h"
#include "io/CAFUtils.h"
#include "analysis/TruthInteractionSummary.h"
#include "plotting/HistogramManager.h"
#include "cuts/DetectorCuts.h"
#include "duneanaobj/StandardRecord/StandardRecord.h" 

class TruthSelection {
public:

  TruthSelection(const SelectionConfig& cfg,
                 const BeamConfig& beam,
                 const DetectorConfig& detector):
                fSelCuts(cfg), fBeam(beam), fDetector(detector) {};

  void SelectTruthInteractions(const caf::StandardRecord& sr,
                               HistogramManager& hist);

  bool IxnPassesTruthLArCuts(const TruthInteractionSummary& summary) const;

  TruthInteractionSummary BuildTruthSummary(
      const caf::SRTrueInteraction& nu) const;

private:

  const SelectionConfig& fSelCuts;
  const BeamConfig& fBeam;
  const DetectorConfig& fDetector;

};