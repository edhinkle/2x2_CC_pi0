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
  bool iscc = false;
  int targetPDG = 0;
  int ixnMode = -999; 

  caf::SRVector3D vertex;

  // -----------------------------
  // Muon kinematics
  // -----------------------------

  double muonEnergy = -999.0;
  double muonCosL   = -999.0;

  // -----------------------------
  // Pi0/Shower counting
  // -----------------------------

  int nPrimPi0        = 0;
  int nSecPi0        = 0;
  int nPrimElectrons = 0;
  int nPrimPhotons = 0;
  int nSecElectrons = 0;
  int nSecPhotons = 0;

  int nPrimPions = 0;

  // -----------------------------
  // Signal flags
  // -----------------------------

  bool passesMx2 = false;
};