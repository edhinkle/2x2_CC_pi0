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
    if (IxnAboveKEThreshold(nu, fDetector))
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

    TruthInteractionSummary summary = BuildSummary(nu);

    // Fill histograms for number of pi0s
    hist.truth.FillPi0Multiplicity(summary.nPi0);

    // Only one pi0
    if (summary.nPi0 != 1)
      continue;
    hist.cuts.Count("Truth", "1 Pi0");

    // Fill histograms for muon kinematics
    hist.truth.FillMuonKinematics(summary.muonCosL, summary.muonEnergy, summary.Numubar);

    // Passes Mx2 signal definition
    if (!summary.passesMx2)
      continue;
    hist.cuts.Count("Truth", "Muon Through Mx2");

    // Fill histograms for neutrino energy
    hist.truth.FillEnu(summary.nuE);

  }
  // Fill histograms for number of interactions above KE threshold for detector
  hist.truth.FillInteractionsAboveKEThresholdPerSpill(ixnsOverKEThreshold);
}

// Helper method to build TruthInteractionSummary
TruthInteractionSummary TruthSelection::BuildSummary(
    const caf::SRTrueInteraction& nu) const
{
  TruthInteractionSummary summary;

  //--------------------------------------------------
  // Detector-level classifications
  //--------------------------------------------------
  summary.KEOverThreshold = IxnAboveKEThreshold(nu, fDetector);

  //--------------------------------------------------
  // Interaction-level info
  //--------------------------------------------------

  summary.nuPDG    = nu.pdg;
  summary.nuE = nu.E;
  summary.vertex   = nu.vtx;
  summary.Numubar = nu.pdg < 0;

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

    if (pdg == 111)
      ++summary.nPi0;

  }

  //--------------------------------------------------
  // Signal definition
  //--------------------------------------------------

  summary.passesMx2 =
      (summary.muonCosL > fSelCuts.minMuonCosL &&
       summary.muonEnergy > fSelCuts.minMuonEnergy);

  return summary;
}