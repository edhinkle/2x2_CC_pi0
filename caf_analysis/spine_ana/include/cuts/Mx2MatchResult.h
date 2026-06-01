#pragma once

#include "TVector3.h"

struct Mx2MatchResult {
  bool trackIsBackward = false;

  double trkMatchDotProdCAF = -999; // dot product of track match in CAF track match branch (for comparison with Mx2 matching score)
  double NDCAFMatchDotProd = -999; // dot product using NDCAFMaker Mx2 Logic to get Mx2 track info
  double mx2TrackEndZ = -999;

  bool isGoodMatch = false;

  int mx2IxnIdx = -1;
  int mx2TrackIdx = -1;

  int LArTrackIdx = -1;

  TVector3 LArTrackDir;       // matched LAr direction
  double LArTrackLength; // matched LAr length

  // Optional truth info (only if MC)
  int truthIxnMx2PartIdx = -1;
  int truthIxnMx2PartType = -1; // primary = 1, secondary = 3
  int truthIxnMx2IxnIdx = -1;
};