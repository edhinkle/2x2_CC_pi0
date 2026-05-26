#pragma once

#include "duneanaobj/StandardRecord/StandardRecord.h"

struct RecoInteractionSummary {


  // -----------------------------
  // Detector-level classifications
  // -----------------------------

  // -----------------------------
  // Interaction-level info
  // -----------------------------

  caf::SRVector3D vertex;

  // -----------------------------
  // Muon kinematics
  // -----------------------------

  double muonCosL   = -999.0;
  double mx2MatchLArStartPosX = -999.0;
  double mx2MatchLArStartPosY = -999.0;
  double mx2MatchLArStartPosZ = -999.0;

  // -----------------------------
  // Shower Info
  // -----------------------------

  // -----------------------------
  // Shower Truth Match Info
  // -----------------------------


  //int nPi0        = 0;

    //TODO: Add pi0 accounting
  
};