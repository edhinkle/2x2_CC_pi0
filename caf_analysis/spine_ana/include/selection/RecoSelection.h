#pragma once
#include "config/ConfigLoader.h"
#include "io/CAFUtils.h"
#include "analysis/RecoInteractionSummary.h"
#include "analysis/TruthInteractionSummary.h"
#include "analysis/MatchedInteractionSummary.h"
#include "cuts/Mx2Matcher.h"
#include "cuts/Mx2MatchResult.h"
#include "cuts/DetectorCuts.h"
#include "selection/TruthSelection.h"
#include "plotting/HistogramManager.h"
#include "duneanaobj/StandardRecord/StandardRecord.h" 

class RecoSelection {
public:

    RecoSelection(const SelectionConfig& cfg,
                const BeamConfig& beam,
                const DetectorConfig& detector,
                const bool mcOnly):
                fSelCuts(cfg), fBeam(beam), fDetector(detector), fMCOnly(mcOnly) {};

    void SelectRecoInteractions(const caf::StandardRecord& sr,
                                HistogramManager& hist);
    
    void FillTruthMatchedCuts(const MatchedInteractionSummary& matchSummary,
                              CutFlowManager& cuts,
                              const std::string& recoCut);

    static double DiffPoints3D(const caf::SRVector3D& v1, const caf::SRVector3D& v2);

private:

    RecoInteractionSummary BuildRecoSummary(
        const caf::SRInteraction& dlpixn,
        const Mx2MatchResult& mx2MatchResult) const;

    MatchedInteractionSummary BuildMatchedIxnSummary(
        const caf::StandardRecord& sr,
        const caf::SRInteraction& dlpixn) const;

    void FillParticleTruthMatching(const caf::StandardRecord& sr,
                                   const caf::SRInteraction& dlpixn,
                                   MatchedInteractionSummary& matchSummary,
                                   RecoInteractionSummary& recoSummary,
                                   const Mx2MatchResult& mx2MatchResult);

    const SelectionConfig& fSelCuts;
    const BeamConfig& fBeam;
    const DetectorConfig& fDetector;
    const bool fMCOnly;

};