#pragma once

#include "duneanaobj/StandardRecord/StandardRecord.h"

struct TruthInteractionSummary {


  // -----------------------------
  // Detector-level classifications
  // -----------------------------
  bool KEOverThreshold = false;

  // -----------------------------
  // Interaction-level info
  // -----------------------------

  int nuPDG = 0;
  double nuE = -999.0;
  bool Numubar = false;

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

  bool passesMx2 = false;
};