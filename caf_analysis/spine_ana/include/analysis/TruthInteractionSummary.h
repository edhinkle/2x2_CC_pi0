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

  caf::SRVector3D vertex;

  // -----------------------------
  // Muon kinematics
  // -----------------------------

  double muonEnergy = -999.0;
  double muonCosL   = -999.0;

  // -----------------------------
  // Multiplicity info
  // -----------------------------

  int nVisibleTracks = 0;

  int nShortTracks = 0;
  int nLongTracks  = 0;

  // -----------------------------
  // Hadron content
  // -----------------------------

  int nPions      = 0;
  int nProtons    = 0;
  int nNeutrons   = 0;

  int nPiPlus     = 0;
  int nPiMinus    = 0;
  int nPi0        = 0;

  int nEscapingPi = 0;

  // -----------------------------
  // Signal flags
  // -----------------------------

  bool passesSignal = false;
};