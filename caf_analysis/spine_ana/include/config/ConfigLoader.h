// ConfigLoader.h
#pragma once

#include <string>
#include "config/SelectionConfig.h"
#include "config/BeamConfig.h"
#include "config/DetectorConfig.h"
#include "config/FluxSystConfig.h"

namespace config {

    SelectionConfig LoadSelectionConfig(const std::string& path);
    BeamConfig LoadBeamConfig(const std::string& path);
    DetectorConfig LoadDetectorConfig(const std::string& path);

    FluxSystConfig LoadFluxSystConfig(const std::string& path);

} // namespace config