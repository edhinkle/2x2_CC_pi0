#pragma once
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include <map>
#include <string>
#include <vector>

struct CutFlow {
  std::vector<std::string> order;      // preserves cut order
  std::map<std::string, double> values; // fast lookup of counts
};


class CutFlowManager {
public:
    void Count(const std::string& flow,
             const std::string& cut,
             double weight = 1.0);

    void Write(TDirectory* dir) const;

private:
    std::map<std::string, CutFlow> cutFlows; // flow name -> cut flow

};