// include/cuts/DetectorCuts.h

#pragma once

#include "duneanaobj/StandardRecord/StandardRecord.h" 
#include "config/DetectorConfig.h"

bool InModuleVolumes(
    const caf::SRVector3D& pt,
    const DetectorConfig& cfg);

bool InFiducialVolume(
    const caf::SRVector3D& pt,
    const DetectorConfig& cfg);

bool TrueIxnAboveKEThreshold(
    const caf::SRTrueInteraction& nu,
    const DetectorConfig& cfg);

bool TrueParticleAboveKEThreshold(
    const caf::SRTrueParticle& part,
    const DetectorConfig& cfg);