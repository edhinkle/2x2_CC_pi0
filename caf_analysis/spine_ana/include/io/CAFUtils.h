// CAFUtils.h
#pragma once

#include <string>

class TChain;

namespace io {

    TChain* BuildCAFChain(const std::string& file_list);

}