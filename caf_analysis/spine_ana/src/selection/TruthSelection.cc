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

    TruthInteractionSummary summary;

    summary.nuPDG = nu.pdg;
    summary.nuE = nu.E;
    summary.vertex = nu.vtx;

    //--------------------------------------------------------------------
    // Primary hadron counting
    //--------------------------------------------------------------------

    int nPiPlus   = 0;
    int nPiMinus  = 0;
    int nPi0      = 0;
    int nProton   = 0;
    int nNeutron  = 0;

    for (const auto& prim : nu.prim) {

      switch (prim.pdg) {

        case 211:
          ++nPiPlus;
          break;

        case -211:
          ++nPiMinus;
          break;

        case 111:
          ++nPi0;
          break;

        case 2212:
          ++nProton;
          break;

        case 2112:
          ++nNeutron;
          break;
      }
    }

    summary.nPiPlus  = nPiPlus;
    summary.nPiMinus = nPiMinus;
    summary.nPi0     = nPi0;
    summary.nProton  = nProton;
    summary.nNeutron = nNeutron;

    //--------------------------------------------------------------------
    // Final-state particle loop
    //--------------------------------------------------------------------

    int nVisibleTracks = 0;
    int nShortTracks   = 0;
    int nLongTracks    = 0;
    int nPions         = 0;
    int nContainedPi   = 0;
    int nEscapingPi    = 0;

    double muonEnergy = -999.0;
    double muonCosL   = -999.0;

    for (const auto& prim : nu.prim) {

      const int pdg = std::abs(prim.pdg);

      //------------------------------------------------------------------
      // Keep only reconstructable charged particles
      //------------------------------------------------------------------

      if (pdg != 13 &&
          pdg != 2212 &&
          pdg != 211 &&
          pdg != 321)
        continue;

      //------------------------------------------------------------------
      // Kinematics
      //------------------------------------------------------------------

      const auto& start = prim.start_pos;
      const auto& end   = prim.end_pos;
      const auto& p4    = prim.p;

      const double px = p4.px;
      const double py = p4.py;
      const double pz = p4.pz;

      const double momentum =
          std::sqrt(px*px + py*py + pz*pz);

      const double dx = end.x - start.x;
      const double dy = end.y - start.y;
      const double dz = end.z - start.z;

      const double length =
          std::sqrt(dx*dx + dy*dy + dz*dz);

      //------------------------------------------------------------------
      // Muon kinematics
      //------------------------------------------------------------------

      if (pdg == 13) {

        muonEnergy = p4.E;

        const double dirX = px / momentum;
        const double dirY = py / momentum;
        const double dirZ = pz / momentum;

        muonCosL =
            dirX * beam.beamX +
            dirY * beam.beamY +
            dirZ * beam.beamZ;

        hist.FillTrueMuonKinematics(muonEnergy, muonCosL);
      }

      //------------------------------------------------------------------
      // Pion containment
      //------------------------------------------------------------------

      if (pdg == 211) {

        ++nPions;

        const bool escapes =
            (std::abs(end.x) > cfg.detectorXMax ||
             std::abs(end.z) > cfg.detectorZMax);

        if (escapes) {
          ++nEscapingPi;
        }
        else {
          ++nContainedPi;
          hist.FillContainedPionLength(length);
        }
      }

      //------------------------------------------------------------------
      // Proton containment
      //------------------------------------------------------------------

      if (pdg == 2212) {

        const bool contained =
            (std::abs(end.x) < cfg.detectorXMax &&
             std::abs(end.z) < cfg.detectorZMax);

        if (contained) {
          hist.FillContainedProtonLength(length);
        }
      }

      //------------------------------------------------------------------
      // Track multiplicity
      //------------------------------------------------------------------

      if (length > cfg.minTrackLength) {

        ++nVisibleTracks;

        if (length < 10.0)
          ++nShortTracks;
        else
          ++nLongTracks;
      }
    }

    //--------------------------------------------------------------------
    // Signal definition
    //--------------------------------------------------------------------

    const bool passesSignal =
        (muonCosL > cfg.minMuonCosL &&
         muonEnergy > cfg.minMuonEnergy);

    if (!passesSignal)
      continue;

    //--------------------------------------------------------------------
    // Store summary quantities
    //--------------------------------------------------------------------

    summary.muonEnergy    = muonEnergy;
    summary.muonCosL      = muonCosL;

    summary.nVisibleTracks = nVisibleTracks;
    summary.nShortTracks   = nShortTracks;
    summary.nLongTracks    = nLongTracks;

    summary.nPions         = nPions;
    summary.nEscapingPi    = nEscapingPi;

    //--------------------------------------------------------------------
    // Fill histograms
    //--------------------------------------------------------------------

    hist.FillTruth(summary);

    hist.FillTruthMultiplicity(nVisibleTracks);
    hist.FillTruthMultiplicity2D(nShortTracks, nLongTracks);

    hist.FillPionMultiplicity(nPions);
    hist.FillEscapingPions(nEscapingPi);

    hist.FillMuonCosTheta(muonCosL);
    hist.FillNeutrinoEnergy(nu.E);
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
    // Proton counting
    //----------------------------------------------

    if (pdg == 2212)
      ++summary.nProtons;

    //----------------------------------------------
    // Pion counting
    //----------------------------------------------

    if (pdg == 211)
      ++summary.nPions;

    //----------------------------------------------
    // Track multiplicity
    //----------------------------------------------

    const double length =
        ComputeTrackLength(prim);

    if (length > fConfig.minTrackLength) {

      ++summary.nVisibleTracks;

      if (length < 10.0)
        ++summary.nShortTracks;
      else
        ++summary.nLongTracks;
    }
  }

  //--------------------------------------------------
  // Signal definition
  //--------------------------------------------------

  summary.passesSignal =
      (summary.muonCosL > fConfig.minMuonCosL &&
       summary.muonEnergy > fConfig.minMuonEnergy);

  return summary;
}