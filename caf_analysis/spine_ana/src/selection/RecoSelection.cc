#pragma once

#include "analysis/RecoInteractionSummary.h"
#include "analysis/TruthInteractionSummary.h"
#include "cuts/DetectorCuts.h"
#include "selection/RecoSelection.h"
#include "cuts/Mx2Matcher.h"
#include "cuts/Mx2MatchResult.h"

void RecoSelection::SelectRecoInteractions(const caf::StandardRecord& sr,
                                           HistogramManager& hist)
{
  //----------------------------------------------------------------------
  // Loop over reco neutrino interactions
  //----------------------------------------------------------------------
  for (long unsigned nixn = 0; nixn < sr->common.ixn.dlp.size(); nixn++) {

    // Interaction index is necessary for matching tracks with Mx2
    dlpixn = sr->common.ixn.dlp[nixn];
    
    //--------------------------------------------------------------------
    // Basic interaction bookkeeping
    //--------------------------------------------------------------------

    hist.cuts.Count("Reco", "All");

    reco_vtx = dlpixn.vtx;

    // Fill histogram for reco vertex distribution (no cuts)
    hist.reco.FillRecoVertexXZNoCuts(reco_vtx.x, reco_vtx.z);

    // Initialize MatchedInteractionSummary for this interaction
    MatchedInteractionSummary matchSummary;
    // Match to Truth if MC And Check if signal
    if (fMcOnly) {
        matchSummary = BuildMatchedSummary(sr, dlpixn);
    }
    // Fill cutflow for matched interactions passing vertex cut and LAr cuts for signal definition
    if (matchSummary.diffVertex < fSelCuts.maxTruthRecoVertexDiff) {
        hist.cuts.Count("RecoMatchedValidVtx", "All");
        if (matchSummary.passesLArCuts){
            hist.cuts.Count("RecoMatchedSignal", "All");
            if (matchSummary.passesMx2) {
                hist.cuts.Count("RecoMatchedSignalwithMx2", "All");
            }
        }
    }

    // Active Volume Cut
    if (!InModuleVolumes(reco_vtx.vtx, fDetector))
      continue;
    hist.cuts.Count("Reco", "Active Volume");
    if (matchSummary.diffVertex < fSelCuts.maxTruthRecoVertexDiff) {
        hist.cuts.Count("RecoMatchedValidVtx", "Active Volume");
        if (matchSummary.passesLArCuts){
            hist.cuts.Count("RecoMatchedSignal", "Active Volume");
            if (matchSummary.passesMx2) {
                hist.cuts.Count("RecoMatchedSignalwithMx2", "Active Volume");
            }
        }
    }

    // Fiducial Volume Cut
    if (!InFiducialVolume(reco_vtx.vtx, fDetector))
      continue;
    hist.cuts.Count("Reco", "Fiducial Volume");
    if (matchSummary.diffVertex < fSelCuts.maxTruthRecoVertexDiff) {
        hist.cuts.Count("RecoMatchedValidVtx", "Fiducial Volume");
        if (matchSummary.passesLArCuts){
            hist.cuts.Count("RecoMatchedSignal", "Fiducial Volume");
            if (matchSummary.passesMx2) {
                hist.cuts.Count("RecoMatchedSignalwithMx2", "Fiducial Volume");
            }
        }
    } 

    // Mx2 Matching Muon Cut
    Mx2MatchResult Mx2Matcher::MatchInteraction(nixn, dlpixn, sr);
    if (!Mx2MatchResult.isGoodMatch)
        continue;
    hist.cuts.Count("Reco", "Mx2 Muon Match");
    if (matchSummary.diffVertex < fSelCuts.maxTruthRecoVertexDiff) {
        hist.cuts.Count("RecoMatchedValidVtx", "Mx2 Muon Match");
        if (matchSummary.passesLArCuts){
            hist.cuts.Count("RecoMatchedSignal", "Mx2 Muon Match");
            if (matchSummary.passesMx2) {
                hist.cuts.Count("RecoMatchedSignalwithMx2", "Mx2 Muon Match");
            }
        }
    } 

    //--------------------------------------------------------------------
    // Event-level reco summary quantities
    //--------------------------------------------------------------------

    RecoInteractionSummary summary = BuildRecoSummary(dlpixn);

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

}

// Helper method to build RecoInteractionSummary
RecoInteractionSummary RecoSelection::BuildRecoSummary(
    const caf::SRInteraction& dlpixn) const
{
  RecoInteractionSummary summary;

  //--------------------------------------------------
  // Detector-level classifications
  //--------------------------------------------------

  //--------------------------------------------------
  // Interaction-level info
  //--------------------------------------------------
  summary.vertex   = dlpixn.vtx;

  //--------------------------------------------------
  // Loop over primaries
  //--------------------------------------------------

  for (const auto& part : dlpixn.part.dlp) {

    const int pdg = std::abs(part.pdg);

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

// Helper method to build MatchedInteractionSummary
MatchedInteractionSummary RecoSelection::BuildMatchedSummary(
    const caf::StandardRecord& sr,
    const caf::SRInteraction& dlpixn) const
{
    MatchedInteractionSummary summary;

    // Match reco interaction to a truth interaction 
    for (long unsigned int ntruth = 0; ntruth < dlpixn.truth.size(); ntruth++) {

          if (summary.bestMatchOverlap < dlpixn.truthOverlap.at(ntruth) && dlpixn.truthOverlap.at(ntruth) > fSelCuts.minTruthIxnOverlap) {
            summary.bestMatchOverlap = dlpixn.truthOverlap.at(ntruth);
            summary.bestMatchIndex = dlpixn.truth.at(ntruth);
          }
    }

    // Check if truth match is a rock interaction by looking at truth id
    if(sr.mc.nu[summary.bestMatchIndex].id > fSelCuts.minRockIxnTruthId) {
      summary.isRockIxn = true;
    }

    // Fill in truth summary for best-matched interaction if it exists
    if (summary.bestMatchIndex != -999) {
        const auto& bestTruthIxn = sr.mc.nu[summary.bestMatchIndex];
        summary.truthSummary = TruthSelection::BuildTruthSummary(bestTruthIxn);
        
        // Get difference in truth and reco vertex positions for best-matched interaction
        summary.diffVertex = DiffPoints3D(summary.truthSummary.vertex, dlpixn.vtx);
        
        summary.passesLArCuts = TruthSelection::IxnPassesTruthLArCuts(summary.truthSummary);
        summary.passesMx2 = summary.truthSummary.passesMx2;
    }

    return summary;
}

double RecoSelection::DiffPoints3D(const caf::SRVector3D& v1, const caf::SRVector3D& v2) const
{
    return TMath::Sqrt((v1.x - v2.x) * (v1.x - v2.x) +
                       (v1.y - v2.y) * (v1.y - v2.y) +
                       (v1.z - v2.z) * (v1.z - v2.z));
}