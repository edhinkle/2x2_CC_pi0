#pragma once
#include "config/ConfigLoader.h"
#include "io/CAFUtils.h"
#include "selection/TruthInteractionSummary.h"
#include "plotting/HistogramManager.h"
#include "duneanaobj/StandardRecord/StandardRecord.h" 

class TruthSelection {
public:

  TruthSelection(const config::SelectionConfig& cfg,
                 const config::BeamConfig& beam,
                 const config::DetectorConfig& detector);

  void SelectTruthInteractions(const caf::StandardRecord& sr,
                               HistogramManager& hist);

private:

  TruthInteractionSummary BuildSummary(
      const caf::SRTrueInteraction& nu);

  bool PassesSignal(
      const TruthInteractionSummary& summary) const;

  const config::SelectionConfig& fSelCuts;
  const config::BeamConfig& fBeam;
  const config::DetectorConfig& fDetector;
};