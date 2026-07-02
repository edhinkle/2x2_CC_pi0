// FluxSystConfig.h
#pragma once

struct FluxSystConfig {
    
    // Flux systematics related file paths
    std::string binFilePath;
    std::string fluxFilePath;

    // Enable flux systematics? 
    bool enableFluxSyst;
    int  nThrows;
    int  seed;

    // Flux Systematics info/settings
    bool isRHC; // whether to use RHC or FHC flux
    int signalNuPDG; // PDG code of the neutrino flavor to apply flux systematics to (e.g. 14 for numu)
    bool nuSignMatters; // whether to distinguish neutrinos vs antineutrinos (e.g. for RHC flux)

};