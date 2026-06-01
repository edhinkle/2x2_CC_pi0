// src/cuts/DetectorCuts.cxx

#include "cuts/DetectorCuts.h"

#include <cmath>

// Check if point is in module boundaries
bool DetectorCuts::InModuleVolumes(
    const caf::SRVector3D& pt,
    const DetectorConfig& cfg) 
{
  if (abs(pt.x) > cfg.absModuleXMax ||
      abs(pt.x) < cfg.absModuleXMin)
    return false;

  if (abs(pt.z) > cfg.absModuleZMax ||
      abs(pt.z) < cfg.absModuleZMin)
    return false;

  if (abs(pt.y) > cfg.absModuleYMax)
    return false;

  return true;
}

// Check if point is in fiducial volume of detector
bool DetectorCuts::InFiducialVolume(
    const caf::SRVector3D& pt,
    const DetectorConfig& cfg) 
{
  if (abs(pt.x) > cfg.absFiducialXMax ||
      abs(pt.x) < cfg.absFiducialXMin)
    return false;

  if (abs(pt.z) > cfg.absFiducialZMax ||
      abs(pt.z) < cfg.absFiducialZMin)
    return false;

  if (abs(pt.y) > cfg.absFiducialYMax)
    return false;

  return true;
}

// Check if true interaction is above KE threshold for detector
bool DetectorCuts::TrueIxnAboveKEThreshold(
    const caf::SRTrueInteraction& nu,
    const DetectorConfig& cfg) 
{
  // Check primaries
  for (const auto& prim : nu.prim) {
    if (TrueParticleAboveKEThreshold(prim, cfg))
      return true;
  }
  // Check secondaries
  for (const auto& sec : nu.sec) {
    if (TrueParticleAboveKEThreshold(sec, cfg))
      return true;
  }
  // If we get here, no particles above threshold
  return false;
}

// Check if true particle is over KE threshold for detector
bool DetectorCuts::TrueParticleAboveKEThreshold(
    const caf::SRTrueParticle& part,
    const DetectorConfig& cfg)
{

  auto start_pos = part.start_pos;

  if (!(InFiducialVolume(start_pos, cfg)))
    return false;

  const double px = part.p.px;
  const double py = part.p.py;
  const double pz = part.p.pz;

  const double momentum =
      std::sqrt(px*px + py*py + pz*pz);

  const double energy = part.p.E;

  const double mass = std::sqrt(energy*energy - momentum*momentum);

  const double KE = energy - mass;

  if (KE > cfg.keThreshold) {
    return true;
  }
  else {
    return false;
  }

}