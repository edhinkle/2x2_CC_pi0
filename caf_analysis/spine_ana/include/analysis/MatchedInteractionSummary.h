#pragma once

#include "selection/TruthInteractionSummary.h"
#include "duneanaobj/StandardRecord/StandardRecord.h"

struct MatchedInteractionSummary {

    // ixn-level info about best-matched truth interaction
    double bestMatchOverlap = -999;
    int bestMatchIndex = -999;

    // Check if rock interaction by looking at truth id
    bool isRockIxn = false;

    TruthInteractionSummary truthSummary;

    double diffVertex = 999; // distance between reco and truth vertex for best match

    bool passesLArCuts = false; // check if best-matched truth interaction passes truth cuts for signal definition
    bool passesMx2 = false; // check if best-matched truth interaction passes Mx2 signal definition cuts

    // -----------------------------
    // Truth Match Particle Info
    // -----------------------------
    bool isPrimary = false;
    int truthMatchMx2TrackPDG = 0;
    double truthMatchMx2TrackE = -1;
    double truthMatchMx2TrackCosL = -999;
    double truthMatchMx2TrackLArStartPosX = -999;
    double truthMatchMx2TrackLArStartPosY = -999;
    double truthMatchMx2TrackLArStartPosZ = -999;
};