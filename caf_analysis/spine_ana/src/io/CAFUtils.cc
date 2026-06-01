// CAFUtils.cc
#include "io/CAFUtils.h"

#include <fstream>
#include <iostream>
#include <string>
#include <TChain.h>

namespace io {

    TChain* BuildCAFChain(const std::string& file_list) {

        // Give input file list
        std::ifstream in(file_list.c_str());

        // Check if input list is present
        if (!in.is_open()) {
            std::cerr << Form("File %s not found", file_list.c_str())
                      << std::endl;
            return nullptr;
        }
    
        // Create new TChain
        TChain* caf_chain = new TChain("cafTree");
    
        // Add files to CAF chain from input list
        std::string file;
        while (in >> file) {
            caf_chain->Add(file.c_str());
        }

        return caf_chain;
    }

}