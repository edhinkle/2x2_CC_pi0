// include/cuts/DetectorCuts.h

#pragma once

#include "duneanaobj/StandardRecord/StandardRecord.h" 
#include "config/DetectorConfig.h"

class DetectorCuts {
public:
    DetectorCuts();
    
    static bool InModuleVolumes(
        const caf::SRVector3D& pt,
        const DetectorConfig& cfg);
    
    static bool InFiducialVolume(
        const caf::SRVector3D& pt,
        const DetectorConfig& cfg);
    
    static bool TrueIxnAboveKEThreshold(
        const caf::SRTrueInteraction& nu,
        const DetectorConfig& cfg);
    
    static bool TrueParticleAboveKEThreshold(
        const caf::SRTrueParticle& part,
        const DetectorConfig& cfg);
};