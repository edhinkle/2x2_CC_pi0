
#pragma once

#include "plotting/CutFlowManager.h"
#include "TDirectory.h"


void CutFlowManager::Count(const std::string& flow,
                           const std::string& cut,
                           double weight)
{
  auto& cf = flows[flow];

  // If this cut has never been seen before, record its order
  if (cf.values.find(cut) == cf.values.end()) {
    cf.order.push_back(cut);
    cf.values[cut] = 0.0;
  }

  // Increment count
  cf.values[cut] += weight;
}

void CutFlowManager::Write(TDirectory* dir) const
{
  dir->cd();

  for (const auto& [flowName, cf] : flows) {

    int nCuts = cf.order.size();
    if (nCuts == 0) continue;

    TH1D* h = new TH1D(flowName.c_str(),
                       flowName.c_str(),
                       nCuts, 0, nCuts);

    for (int i = 0; i < nCuts; ++i) {
      const std::string& cutName = cf.order[i];
      double value = cf.values.at(cutName);

      int bin = i + 1;
      h->SetBinContent(bin, value);
      h->GetXaxis()->SetBinLabel(bin, cutName.c_str());
    }

    h->Write();
  }
}