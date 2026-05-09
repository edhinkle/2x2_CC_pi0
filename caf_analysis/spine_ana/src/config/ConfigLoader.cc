// ConfigLoader.cc
#include "config/ConfigLoader.h"
#include <yaml-cpp/yaml.h>

namespace config {

    SelectionConfig LoadSelectionConfig(const std::string& path) {
        SelectionConfig cfg;
    
        // Load config yaml
        YAML::Node config = YAML::LoadFile("config_data/CC1pi0.yaml");
    
        cfg.cosThetaCut = config["selection"]["cosThetaCut"].as<double>();
        cfg.muonEnergyCut = config["selection"]["muonEnergyCut"].as<double>();
    
        return cfg;
    }
    
    BeamConfig LoadBeamConfig(const std::string& path) {
        BeamConfig cfg;
    
        // Load config yaml
        YAML::Node config = YAML::LoadFile("config_data/CC1pi0.yaml");
    
        // Load from config yaml
        cfg.beam_x = config["beam"]["beam_x"].as<double>();
        cfg.beam_y = config["beam"]["beam_y"].as<double>();
        cfg.beam_z = config["beam"]["beam_z"].as<double>();
        cfg.dir = TVector3(cfg.beam_x, cfg.beam_y, cfg.beam_z);
    
        // Get unit direction vector for XYZ
        auto beam_mag = TMath::Sqrt(cfg.beam_x * cfg.beam_x + cfg.beam_y * cfg.beam_y + cfg.beam_z * cfg.beam_z);
        cfg.beam_x = cfg.beam_x / beam_mag;
        cfg.beam_y = cfg.beam_y / beam_mag;
        cfg.beam_z = cfg.beam_z / beam_mag;
    
        return cfg;
    }

    DetectorConfig LoadDetectorConfig(const std::string& path) {
        DetectorConfig cfg;
    
        // Load config yaml
        YAML::Node config = YAML::LoadFile("config_data/CC1pi0.yaml");
    
        // Can we see the particle?
        cfg.keThreshold = config["detector"]["keThreshold"].as<double>();
        cfg.minTrackLength = config["detector"]["minTrackLength"].as<double>();

        // Mx2 offsets for matching (misalignment in physical detectors)
        cfg.mx2OffsetXTruth = config["detector"]["mx2OffsetTruth"].as<double>();
        cfg.mx2OffsetYTruth = config["detector"]["mx2OffsetTruth"].as<double>();
        cfg.mx2OffsetXReco = config["detector"]["mx2OffsetReco"].as<double>();
        cfg.mx2OffsetYReco = config["detector"]["mx2OffsetReco"].as<double>();

        // Module boundaries for FV cuts
        // First, get module boundaries (abs value; assumes symmetric in + and - directions for XYZ)
        cfg.absModuleXMin = config["detector"]["absModuleMinX"].as<double>();
        cfg.absModuleXMax = config["detector"]["absModuleMaxX"].as<double>();
        cfg.absModuleYMin = config["detector"]["absModuleMinY"].as<double>();
        cfg.absModuleYMax = config["detector"]["absModuleMaxY"].as<double>();
        cfg.absModuleZMin = config["detector"]["absModuleMinZ"].as<double>();
        cfg.absModuleZMax = config["detector"]["absModuleMaxZ"].as<double>();

        // Then, get FV cuts
        detFVCutAbsMinX = config["detector"]["detFVCutAbsMinX"].as<double>();
        detFVCutAbsMaxX = config["detector"]["detFVCutAbsMaxX"].as<double>();
        detFVCutAbsMinY = config["detector"]["detFVCutAbsMinY"].as<double>();
        detFVCutAbsMaxY = config["detector"]["detFVCutAbsMaxY"].as<double>();
        detFVCutAbsMinZ = config["detector"]["detFVCutAbsMinZ"].as<double>();
        detFVCutAbsMaxZ = config["detector"]["detFVCutAbsMaxZ"].as<double>();

        // Finally, apply FV cuts
        cfg.absFiducialXMin = cfg.absModuleXMin + detFVCutAbsMinX;
        cfg.absFiducialXMax = cfg.absModuleXMax - detFVCutAbsMaxX;
        cfg.absFiducialYMin = cfg.absModuleYMin - detFVCutAbsMinY; // symmetry in y (+ no gaps)
        cfg.absFiducialYMax = cfg.absModuleYMax - detFVCutAbsMaxY; // symmetry in y (+ no gaps)
        cfg.absFiducialZMin = cfg.absModuleZMin + detFVCutAbsMinZ;
        cfg.absFiducialZMax = cfg.absModuleZMax - detFVCutAbsMaxZ;

    
        return cfg;
    }

} // namespace config