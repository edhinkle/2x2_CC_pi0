// ConfigLoader.cc
#include "config/ConfigLoader.h"
#include <yaml-cpp/yaml.h>

namespace config {

    SelectionConfig LoadSelectionConfig(const std::string& path) {
        SelectionConfig cfg;
    
        // Pseudocode
        YAML::Node config = YAML::LoadFile("config_data/CC1pi0.yaml");
    
        cfg.minTrackLength = config["selection"]["minTrackLength"].as<double>();
        cfg.cosThetaCut = config["selection"]["cosThetaCut"].as<double>();
        cfg.muonEnergyCut = config["selection"]["muonEnergyCut"].as<double>();
    
        return cfg;
    }
    
    BeamConfig LoadBeamConfig(const std::string& path) {
        BeamConfig cfg;
    
        /// Pseudocode
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

} // namespace config