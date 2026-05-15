#pragma once
#include "config/ConfigLoader.h"
#include "io/CAFUtils.h"
#include "analysis/RecoInteractionSummary.h"
#include "analysis/TruthInteractionSummary.h"
#include "analysis/MatchedInteractionSummary.h"
#include "selection/TruthSelection.h"
#include "plotting/HistogramManager.h"
#include "duneanaobj/StandardRecord/StandardRecord.h" 

class RecoSelection {
public:

    RecoSelection(const config::SelectionConfig& cfg,
                const config::BeamConfig& beam,
                const config::DetectorConfig& detector,
                bool mcOnly);

    void SelectRecoInteractions(const caf::StandardRecord& sr,
                              HistogramManager& hist);

    double DiffPoints3D(const caf::SRVector3D& v1, const caf::SRVector3D& v2) const;

private:

    RecoInteractionSummary BuildRecoSummary(
        const caf::SRInteraction& dlpixn);

    MatchedInteractionSummary BuildMatchedSummary(
        const caf::StandardRecord& sr,
        const caf::SRInteraction& dlpixn);

    const config::SelectionConfig& fSelCuts;
    const config::BeamConfig& fBeam;
    const config::DetectorConfig& fDetector;
    const bool fMCOnly;
};