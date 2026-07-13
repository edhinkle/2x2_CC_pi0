// DetectorConfig.h
#pragma once

struct DetectorConfig {

    // Can we see the particle?
    double keThreshold; // GeV, kinetic energy threshold for seeing particle
    double minTrackLength;

    // Detector boundaries
    double absModuleXMin; // cm
    double absModuleXMax; // cm
    double absModuleYMin; // cm
    double absModuleYMax; // cm
    double absModuleZMin; // cm
    double absModuleZMax; // cm
    double absCathodeLocX; // cm, absolute X position of cathode planes in modules

    // Module boundaries for FV cuts
    double absFiducialXMin; // cm
    double absFiducialXMax; // cm
    double absFiducialYMin; // cm
    double absFiducialYMax; // cm
    double absFiducialZMin; // cm
    double absFiducialZMax; // cm

    // Cathode Buffer min/max x
    double absCathodeBufferXMin; // cm, minimum X position of cathode buffer region
    double absCathodeBufferXMax; // cm, maximum X position of cathode buffer region

    // FV cut for Mx2 exiting track requirement
    double absFiducialZMaxExiting; // cm, max Z for vertex for Mx2 exiting track requirement

    // Mx2 boundaries from DUNE DocDB 32440
    double mx2DownstreamZMin; // cm
    double mx2DownstreamZMax; // cm
    double mx2UpstreamZMin; // cm
    double mx2UpstreamZMax; // cm
};