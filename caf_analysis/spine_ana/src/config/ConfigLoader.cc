// ConfigLoader.cc
#include "config/ConfigLoader.h"
#include "yaml-cpp/yaml.h"

namespace config {

    SelectionConfig LoadSelectionConfig(const std::string& path) {
        SelectionConfig cfg;
    
        // Load config yaml
        YAML::Node config = YAML::LoadFile(path);
    
        cfg.cosThetaCut = config["selection"]["cosThetaCut"].as<double>();
        cfg.muonEnergyCut = config["selection"]["muonEnergyCut"].as<double>();

        cfg.minTruthIxnOverlap = config["selection"]["minTruthIxnOverlap"].as<double>();
        cfg.maxRockIxnTruthId = config["selection"]["maxRockIxnTruthId"].as<double>();
        cfg.maxTruthRecoVertexDiff = config["selection"]["maxTruthRecoVertexDiff"].as<double>();
    
        cfg.maxMx2TrackVertexDiff = config["selection"]["maxMx2TrackVertexDiff"].as<double>();
        cfg.maxTrkMatchMx2MatchDiff = config["selection"]["maxTrkMatchMx2MatchDiff"].as<double>();

        cfg.mx2MatchExtrapLocZ = config["selection"]["mx2MatchExtrapLocZ"].as<double>();
        cfg.maxMx2MatchExtrapDiffX = config["selection"]["maxMx2MatchExtrapDiffX"].as<double>();
        cfg.maxMx2MatchExtrapDiffY = config["selection"]["maxMx2MatchExtrapDiffY"].as<double>();
        cfg.maxMx2MatchArctanDiffXZ = config["selection"]["maxMx2MatchArctanDiffXZ"].as<double>();
        cfg.maxMx2MatchArctanDiffYZ = config["selection"]["maxMx2MatchArctanDiffYZ"].as<double>();

        cfg.minMx2TrackLength = config["selection"]["minMx2TrackLength"].as<double>();

        // Mx2 Track Start Limit
        double mx2EnterDownstreamAllowance = config["selection"]["mx2EnterDownstreamAllowance"].as<double>();
        double mx2DownstreamZMin = config["detector"]["mx2DownstreamZMin"].as<double>(); // cm

        cfg.maxMx2MatchTrackStartZ = mx2DownstreamZMin + mx2EnterDownstreamAllowance;

        // Mx2 Track End Limit (for muon)
        double mx2DistFromEndDownstreamZ = config["selection"]["mx2DistFromEndDownstreamZ"].as<double>();
        double mx2DownstreamZMax = config["detector"]["mx2DownstreamZMax"].as<double>(); // cm

        cfg.minMx2MatchTrackEndZ = mx2DownstreamZMax - mx2DistFromEndDownstreamZ;

        // Mx2 Matching Threshold
        cfg.mx2TrackMatchDotProdThreshold = config["selection"]["mx2TrackMatchDotProdThreshold"].as<double>();

        return cfg;
    }
    
    BeamConfig LoadBeamConfig(const std::string& path) {
        BeamConfig cfg;
    
        // Load config yaml
        YAML::Node config = YAML::LoadFile(path);
    
        // Load from config yaml
        cfg.beam_x = config["beam"]["beam_x"].as<double>();
        cfg.beam_y = config["beam"]["beam_y"].as<double>();
        cfg.beam_z = config["beam"]["beam_z"].as<double>();
    
        // Get unit direction vector for XYZ
        auto beam_mag = TMath::Sqrt(cfg.beam_x * cfg.beam_x + cfg.beam_y * cfg.beam_y + cfg.beam_z * cfg.beam_z);
        cfg.beam_x = cfg.beam_x / beam_mag;
        cfg.beam_y = cfg.beam_y / beam_mag;
        cfg.beam_z = cfg.beam_z / beam_mag;
        cfg.beam_dir = TVector3(cfg.beam_x, cfg.beam_y, cfg.beam_z);
    
        return cfg;
    }

    DetectorConfig LoadDetectorConfig(const std::string& path) {
        DetectorConfig cfg;
    
        // Load config yaml
        YAML::Node config = YAML::LoadFile(path);
    
        // Can we see the particle?
        cfg.keThreshold = config["detector"]["keThreshold"].as<double>();
        cfg.minTrackLength = config["detector"]["minTrackLength"].as<double>();

        // Module boundaries for FV cuts
        // First, get module boundaries (abs value; assumes symmetric in + and - directions for XYZ)
        cfg.absModuleXMin = config["detector"]["absModuleMinX"].as<double>();
        cfg.absModuleXMax = config["detector"]["absModuleMaxX"].as<double>();
        cfg.absModuleYMin = config["detector"]["absModuleMinY"].as<double>();
        cfg.absModuleYMax = config["detector"]["absModuleMaxY"].as<double>();
        cfg.absModuleZMin = config["detector"]["absModuleMinZ"].as<double>();
        cfg.absModuleZMax = config["detector"]["absModuleMaxZ"].as<double>();

        // Then, get FV cuts
        double detFVCutAbsMinX = config["selection"]["detFVCutAbsMinX"].as<double>();
        double detFVCutAbsMaxX = config["selection"]["detFVCutAbsMaxX"].as<double>();
        double detFVCutAbsMinY = config["selection"]["detFVCutAbsMinY"].as<double>();
        double detFVCutAbsMaxY = config["selection"]["detFVCutAbsMaxY"].as<double>();
        double detFVCutAbsMinZ = config["selection"]["detFVCutAbsMinZ"].as<double>();
        double detFVCutAbsMaxZ = config["selection"]["detFVCutAbsMaxZ"].as<double>();
        double detFVCutAbsMaxZExiting = config["selection"]["detFVCutAbsMaxZExiting"].as<double>(); // Mx2 exiting track requirement

        // Apply FV cuts
        cfg.absFiducialXMin = cfg.absModuleXMin + detFVCutAbsMinX;
        cfg.absFiducialXMax = cfg.absModuleXMax - detFVCutAbsMaxX;
        cfg.absFiducialYMin = cfg.absModuleYMin - detFVCutAbsMinY; // symmetry in y (+ no gaps)
        cfg.absFiducialYMax = cfg.absModuleYMax - detFVCutAbsMaxY; // symmetry in y (+ no gaps)
        cfg.absFiducialZMin = cfg.absModuleZMin + detFVCutAbsMinZ;
        cfg.absFiducialZMax = cfg.absModuleZMax - detFVCutAbsMaxZ;

        // Mx2 exiting track requirement:
        cfg.absFiducialZMaxExiting = cfg.absModuleZMax - detFVCutAbsMaxZExiting;

        // Mx2 boundaries from DUNE DocDB 32440
        cfg.mx2DownstreamZMin = config["detector"]["mx2DownstreamZMin"].as<double>(); // cm
        cfg.mx2DownstreamZMax = config["detector"]["mx2DownstreamZMax"].as<double>(); // cm
        cfg.mx2UpstreamZMin = config["detector"]["mx2UpstreamZMin"].as<double>(); // cm
        cfg.mx2UpstreamZMax = config["detector"]["mx2UpstreamZMax"].as<double>(); // cm

        return cfg;
    }

} // namespace config