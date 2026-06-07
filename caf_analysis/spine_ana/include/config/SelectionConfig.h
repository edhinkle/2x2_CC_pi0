// SelectionConfig.h
#pragma once

struct SelectionConfig {
    
    // Truth Cuts to select muons passing through Mx2
    double cosThetaCut;
    double muonEnergyCut;

    // Minimum overlap to consider Reco/Truth Ixn match successful
    double minTruthIxnOverlap;
    double maxRockIxnTruthId;
    double maxTruthRecoVertexDiff;

    // Mx2 Matching cuts (Allowances vertex and track Match)
    double maxMx2TrackVertexDiff;
    double maxTrkMatchMx2MatchDiff;

    // Values for Mx2 matching (misalignment in physical detectors -- NDCAFMaker copied script)
    double mx2MatchExtrapLocZ; // cm, Z location to which tracks in LAr and 2x2 are extrapolated
    double maxMx2MatchExtrapDiffX; // cm, allowance for difference in X of extrapolated LAr and Mx2 tracks
    double maxMx2MatchExtrapDiffY;  // cm, allowance for difference in Y of extrapolated LAr and Mx2 tracks
    double maxMx2MatchArctanDiffXZ;  // rad, allowance for differnce in arctangent in X/Z plane for LAr and Mx2 tracks
    double maxMx2MatchArctanDiffYZ;  // rad, allowance for differnce in arctangent in Y/Z plane for LAr and Mx2 tracks

    // Values for Mx2 by-hand matching tracks
    double minMx2TrackLength; // cm, minimum length of track in Mx2 for matching
    double maxMx2MatchTrackStartZ; // cm, maximum value for Mx2 track start point in Z which allows for assumption that track starts outside Mx2
    double minMx2MatchTrackEndZ; // cm, minimum value for Mx2 track end point in Z which allows for assumption that track is muon

    // Mx2 Matched Track Criteria
    double mx2TrackMatchDotProdThreshold; // Minimum dot product between direction of Mx2 and 2x2 track to count as match
};