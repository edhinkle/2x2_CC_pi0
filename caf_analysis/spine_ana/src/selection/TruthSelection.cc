#include "analysis/TruthInteractionSummary.h"
#include "cuts/DetectorCuts.h"
#include "selection/TruthSelection.h"

void TruthSelection::SelectTruthInteractions(const caf::StandardRecord& sr,
                                             HistogramManager& hist)
{
  int ixnsOverKEThreshold = 0;
  //----------------------------------------------------------------------
  // Loop over truth neutrino interactions
  //----------------------------------------------------------------------
  for (const auto& nu : sr.mc.nu) {

    // Check if ixn is above KE threshold for detector
    if (TrueIxnAboveKEThreshold(nu, fDetector))
      ++ixnsOverKEThreshold;
    
    //--------------------------------------------------------------------
    // Basic interaction bookkeeping
    //--------------------------------------------------------------------

    hist.cuts.Count("Truth", "All");

    // Argon target only
    if (nu.targetPDG != 1000180400)
      continue;
    hist.cuts.Count("Truth", "Argon Target");

    // Active Volume Cut
    if (!InModuleVolumes(nu.vtx, fDetector))
      continue;
    hist.cuts.Count("Truth", "Active Volume");

    // Numu Cut
    if (std::abs(nu.pdg) != 14)
      continue;
    hist.cuts.Count("Truth", "NuMu");

    // Fiducial Volume Cut
    if (!InFiducialVolume(nu.vtx, fDetector))
      continue;
    hist.cuts.Count("Truth", "Fiducial Volume");

    // CC only
    if (!nu.iscc)
      continue;
    hist.cuts.Count("Truth", "CC");


    //--------------------------------------------------------------------
    // Event-level truth summary quantities
    //--------------------------------------------------------------------

    TruthInteractionSummary summary = BuildTruthSummary(nu);
    
    // Fill histograms for number of pi0s
    hist.truth.FillPrimPi0Multiplicity(summary.nPrimPi0);

    // Only one pi0
    if (summary.nPrimPi0 != 1)
      continue;
    hist.cuts.Count("Truth", "1 Pi0");

    // Fill histograms for muon kinematics and secondary pi0 multiplicity before Mx2 cuts
    hist.truth.FillMuonKinematics(summary.muonCosL, summary.muonEnergy, summary.Numubar, summary.passesMx2);
    hist.truth.FillSecPi0MultiplicityPreMx2(summary.nSecPi0);
    hist.truth.FillShowerMultiplicityPreMx2(summary.nPrimElectrons, summary.nSecElectrons,
                                            summary.nPrimPhotons, summary.nSecPhotons);

    // Passes Mx2 signal definition
    if (!summary.passesMx2)
      continue;
    hist.cuts.Count("Truth", "Muon Through Mx2");

    // Fill histograms for neutrino energy
    hist.truth.FillEnu(summary.nuE);

    // Fill post-Mx2 multiplicity histograms
    hist.truth.FillSecPi0MultiplicityPostMx2(summary.nSecPi0);
    hist.truth.FillShowerMultiplicityPostMx2(summary.nPrimElectrons, summary.nSecElectrons,
                                             summary.nPrimPhotons, summary.nSecPhotons);

  }
  // Fill histograms for number of interactions above KE threshold for detector
  hist.truth.FillInteractionsAboveKEThresholdPerSpill(ixnsOverKEThreshold);
}

// Helper method to build TruthInteractionSummary
TruthInteractionSummary TruthSelection::BuildTruthSummary(
    const caf::SRTrueInteraction& nu) const
{
  TruthInteractionSummary summary;

  //--------------------------------------------------
  // Detector-level classifications
  //--------------------------------------------------
  summary.KEOverThreshold = TrueIxnAboveKEThreshold(nu, fDetector);

  //--------------------------------------------------
  // Interaction-level info
  //--------------------------------------------------

  summary.nuPDG    = nu.pdg;
  summary.nuE = nu.E;
  summary.vertex   = nu.vtx;
  summary.Numubar = nu.pdg < 0;
  summary.iscc    = nu.iscc;
  summary.targetPDG = nu.targetPDG;
  summary.ixnMode = nu.mode;

  //--------------------------------------------------
  // Loop over primaries
  //--------------------------------------------------

  for (const auto& prim : nu.prim) {

    const int pdg = std::abs(prim.pdg);

    //----------------------------------------------
    // Muon
    //----------------------------------------------

    if (pdg == 13) {

      const auto& p4 = prim.p;

      const double p =
          std::sqrt(p4.px*p4.px +
                    p4.py*p4.py +
                    p4.pz*p4.pz);

      summary.muonEnergy = p4.E;

      const double dirX = p4.px / p;
      const double dirY = p4.py / p;
      const double dirZ = p4.pz / p;

      summary.muonCosL =
          dirX * fBeam.beamX +
          dirY * fBeam.beamY +
          dirZ * fBeam.beamZ;
    }

    //----------------------------------------------
    // Pi0 counting
    //----------------------------------------------

    if (abs(pdg) == 111)
      ++summary.nPrimPi0;

    if (abs(pdg) == 11)
      ++summary.nPrimElectrons;

    if (abs(pdg) == 22)
      ++summary.nPrimPhotons;

  }

  //--------------------------------------------------
  // Loop over secondaries
  //--------------------------------------------------
  for (const auto& sec : nu.sec) {

    const int pdg = std::abs(sec.pdg);

    if (abs(pdg) == 111)
      ++summary.nSecPi0;

    if (abs(pdg) == 11)
      ++summary.nSecElectrons;

    if (abs(pdg) == 22)
      ++summary.nSecPhotons;

  }

  //--------------------------------------------------
  // Signal definition
  //--------------------------------------------------

  summary.passesMx2 =
      (summary.muonCosL > fSelCuts.minMuonCosL &&
       summary.muonEnergy > fSelCuts.minMuonEnergy);

  return summary;
}

// Check if truth interaction summary passes truth cuts for signal definition
bool TruthSelection::IxnPassesTruthLArCuts(const TruthInteractionSummary& summary) const
{
  // Argon target only
  if (summary.targetPDG != 1000180400)
    return false;

  // Active Volume Cut
  if (!InModuleVolumes(summary.vertex, fDetector))
    return false;

  // Numu Cut
  if (std::abs(summary.nuPDG) != 14)
    return false;

  // Fiducial Volume Cut
  if (!InFiducialVolume(summary.vertex, fDetector))
    return false;

  // CC only
  if (!summary.iscc)
    return false;

  // Only one pi0
  if (summary.nPrimPi0 != 1)
    return false;

  return true;
}
