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

  double muonEnergy = -999.0;
  double muonCosL   = -999.0;

  // -----------------------------
  // Hadron content
  // -----------------------------

  int nPi0        = 0;


  // -----------------------------
  // Signal flags
  // -----------------------------

  bool passesRecoMx2 = false;
};