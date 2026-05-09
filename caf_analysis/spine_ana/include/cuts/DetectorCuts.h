// include/cuts/DetectorCuts.h

#pragma once

#include "duneanaobj/StandardRecord/StandardRecord.h" 
#include "config/DetectorConfig.h"

bool InFiducialVolume(
    const caf::SRVector3D& pt,
    const DetectorConfig& cfg);

bool IxnAboveKEThreshold(
    const caf::SRTrueInteraction& nu,
    const DetectorConfig& cfg);