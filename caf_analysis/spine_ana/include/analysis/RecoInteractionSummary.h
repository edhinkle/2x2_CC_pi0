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

  // -----------------------------
  // Hadron content
  // -----------------------------

  //int nPi0        = 0;

    //TODO: Add pi0 accounting
  // -----------------------------
  // Signal flags
  // -----------------------------

  
};