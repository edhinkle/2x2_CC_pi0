// DetectorConfig.h
#pragma once

struct DetectorConfig {

    // Can we see the particle?
    double keThreshold; // GeV, kinetic energy threshold for seeing particle
    double minTrackLength;

    // Mx2 offsets for matching (misalignment in physical detectors)
    double mx2OffsetXTruth; // cm, offset in X for Mx2 matching for truth particles
    double mx2OffsetYTruth; // cm, offset in Y for Mx2 matching for truth particles
    double mx2OffsetXReco;  // cm, offset in X for Mx2 matching for reco particles
    double mx2OffsetYReco;  // cm, offset in Y for Mx2 matching for reco particles

    // Detector boundaries
    double absModuleXMin; // cm
    double absModuleXMax; // cm
    double absModuleYMin; // cm
    double absModuleYMax; // cm
    double absModuleZMin; // cm
    double absModuleZMax; // cm

    // Module boundaries for FV cuts
    double absFiducialXMin; // cm
    double absFiducialXMax; // cm
    double absFiducialYMin; // cm
    double absFiducialYMax; // cm
    double absFiducialZMin; // cm
    double absFiducialZMax; // cm


};