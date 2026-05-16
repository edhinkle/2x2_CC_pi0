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

    // Initialize MatchedInteractionSummary for this interaction -- only used if MC
    MatchedInteractionSummary matchSummary;
    // Match to Truth if MC And Check if signal
    if (fMcOnly) {
        matchSummary = BuildMatchedSummary(sr, dlpixn);
    }
    // Fill cutflow for matched interactions passing vertex cut and LAr cuts for signal definition
    if (fMcOnly) FillTruthMatchedCuts(matchSummary, hist.cuts, "All");

    // Active Volume Cut
    if (!InModuleVolumes(reco_vtx.vtx, fDetector))
      continue;
    hist.cuts.Count("Reco", "Active Volume");
    if (fMcOnly) FillTruthMatchedCuts(matchSummary, hist.cuts, "Active Volume");

    // Fiducial Volume Cut
    if (!InFiducialVolume(reco_vtx.vtx, fDetector))
      continue;
    hist.cuts.Count("Reco", "Fiducial Volume");
    if (fMcOnly) FillTruthMatchedCuts(matchSummary, hist.cuts, "Fiducial Volume");

    // Mx2 Matching Muon Cut
    Mx2MatchResult Mx2Matcher::MatchInteraction(nixn, dlpixn, sr);
    if (!Mx2MatchResult.isGoodMatch)
        continue;
    hist.cuts.Count("Reco", "Mx2 Muon Match");
    if (fMcOnly) FillTruthMatchedCuts(matchSummary, hist.cuts, "Mx2 Muon Match");

    //--------------------------------------------------------------------
    // Event-level reco summary quantities
    //--------------------------------------------------------------------

    RecoInteractionSummary recoSummary = BuildRecoSummary(dlpixn, Mx2MatchResult);

    // TODO: Add Cut on Pi0s

    // Fill reco histograms
    hist.reco.FillRecoCosMuonAngle(recoSummary.muonCosL);

  }

}

// Helper method to build RecoInteractionSummary
RecoInteractionSummary RecoSelection::BuildRecoSummary(
    const caf::SRInteraction& dlpixn,
    const Mx2MatchResult&) const
{
  RecoInteractionSummary summary;

  //--------------------------------------------------
  // Detector-level classifications
  //--------------------------------------------------

  //--------------------------------------------------
  // Interaction-level info
  //--------------------------------------------------
  summary.vertex = dlpixn.vtx;

  //--------------------------------------------------
  // Mx2 Match info
  //--------------------------------------------------
  summary.muonCosL = Mx2MatchResult.LArTrackDir.Dot(fBeam.beam_dir);

  //--------------------------------------------------
  // Loop over primaries
  //--------------------------------------------------

  for (const auto& part : dlpixn.part.dlp) {

    const int pdg = std::abs(part.pdg);
    //TODO: Add pi0 accounting

    //----------------------------------------------
    // Showers
    //----------------------------------------------

    //if (pdg == 13) {
//
    //  const auto& p4 = prim.p;
//
    //  const double p =
    //      std::sqrt(p4.px*p4.px +
    //                p4.py*p4.py +
    //                p4.pz*p4.pz);
//
    //  summary.muonEnergy = p4.E;
//
    //  const double dirX = p4.px / p;
    //  const double dirY = p4.py / p;
    //  const double dirZ = p4.pz / p;
//
    //  summary.muonCosL =
    //      dirX * fBeam.beamX +
    //      dirY * fBeam.beamY +
    //      dirZ * fBeam.beamZ;
    //}

  }

  //--------------------------------------------------
  // Signal definition
  //--------------------------------------------------

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

void RecoSelection::FillTruthMatchedCuts(const MatchedInteractionSummary& matchSummary,
                     CutFlowManager& cuts,
                     const std::string& recoCut)
{
    if (matchSummary.diffVertex >= fSelCuts.maxTruthRecoVertexDiff)
        return;

    cuts.Count("RecoMatchedValidVtx", recoCut);

    if (!matchSummary.passesLArCuts)
        return;

    cuts.Count("RecoMatchedSignal", recoCut);

    if (!matchSummary.passesMx2)
        return;

    cuts.Count("RecoMatchedSignalwithMx2", recoCut);
}

double RecoSelection::DiffPoints3D(const caf::SRVector3D& v1, const caf::SRVector3D& v2) const
{
    return TMath::Sqrt((v1.x - v2.x) * (v1.x - v2.x) +
                       (v1.y - v2.y) * (v1.y - v2.y) +
                       (v1.z - v2.z) * (v1.z - v2.z));
}