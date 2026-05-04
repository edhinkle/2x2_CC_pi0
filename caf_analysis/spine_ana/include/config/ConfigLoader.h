// ConfigLoader.h
#pragma once

#include <string>
#include "SelectionConfig.h"
#include "BeamConfig.h"

namespace config {

    SelectionConfig LoadSelectionConfig(const std::string& path);
    BeamConfig LoadBeamConfig(const std::string& path);

} // namespace config