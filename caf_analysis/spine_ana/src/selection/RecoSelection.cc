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
  // Initialize Object for Mx2Matcher
  //----------------------------------------------------------------------
  Mx2Matcher mx2Matcher(fSelCuts, fBeam, fDetector, fMCOnly);
  //----------------------------------------------------------------------
  // Loop over reco neutrino interactions
  //----------------------------------------------------------------------
  for (long unsigned nixn = 0; nixn < sr.common.ixn.dlp.size(); nixn++) {

    // Interaction index is necessary for matching tracks with Mx2
    auto dlpixn = sr.common.ixn.dlp[nixn];
    
    //--------------------------------------------------------------------
    // Basic interaction bookkeeping
    //--------------------------------------------------------------------

    hist.cuts.Count("Reco", "All");

    auto reco_vtx = dlpixn.vtx;

    // Fill histogram for reco vertex distribution (no cuts)
    hist.reco.FillRecoVertexNoCuts(reco_vtx);

    // Initialize MatchedInteractionSummary for this interaction -- only used if MC
    MatchedInteractionSummary matchSummary;
    // Match to Truth if MC And Check if signal
    if (fMCOnly) {
        matchSummary = BuildMatchedIxnSummary(sr, dlpixn);
    }
    // Fill cutflow for matched interactions passing vertex cut and LAr cuts for signal definition
    if (fMCOnly) FillTruthMatchedCuts(matchSummary, hist.cuts, "All");

    // Active Volume Cut
    if (!DetectorCuts::InModuleVolumes(reco_vtx, fDetector))
      continue;
    hist.cuts.Count("Reco", "Active Volume");
    if (fMCOnly) FillTruthMatchedCuts(matchSummary, hist.cuts, "Active Volume");

    // Fiducial Volume Cut
    if (!DetectorCuts::InFiducialVolume(reco_vtx, fDetector))
      continue;
    hist.cuts.Count("Reco", "Fiducial Volume");
    if (fMCOnly) FillTruthMatchedCuts(matchSummary, hist.cuts, "Fiducial Volume");

    // Mx2 Matching Muon Cut
    Mx2MatchResult mx2RecoMatch = mx2Matcher.MatchInteraction(nixn, dlpixn, sr);
    if (!mx2RecoMatch.isGoodMatch)
        continue;
    hist.cuts.Count("Reco", "Mx2 Muon Match");
    if (fMCOnly) FillTruthMatchedCuts(matchSummary, hist.cuts, "Mx2 Muon Match");

    //--------------------------------------------------------------------
    // Event-level reco summary quantities
    //--------------------------------------------------------------------

    RecoInteractionSummary recoSummary = BuildRecoSummary(dlpixn, mx2RecoMatch);

    // Updated MatchedInteractionSummary with Particle Truth Matching from Mx2, etc.
    if (fMCOnly) {
        FillParticleTruthMatching(sr, dlpixn, matchSummary, recoSummary, mx2RecoMatch);
        // Fill truth match Mx2 Track (and muon from truth summary)
        hist.truthMatch.FillTruthMatchMx2TrackInfo(matchSummary,recoSummary);
        hist.truthMatch.FillTruthMatchDiffTruthRecoVertex(matchSummary);
        // Fill histograms for shower multiplicity for best matched interaction
        hist.truthMatch.FillTruthMatchIxnShowerMultiplicity(matchSummary,recoSummary);
        // Fill histograms for counting ixn mode
        hist.truthMatch.FillTruthMatchIxnIxnMode(matchSummary);
        // Fill histograms for counting selected ixns where best match truth is rock
        hist.truthMatch.FillTruthMatchIxnIsRock(matchSummary);
        // Fill histograms for counting CC/NC, nu PDG for rock and non-rock truth-matched BKG
        hist.truthMatch.FillTruthMatchIxnIsBkgCCnuPDG(matchSummary);
    }
    // TODO: Add Cut on Pi0s

    // Fill reco histograms
    hist.reco.FillRecoVertexWithCuts(recoSummary.vertex);
    hist.reco.FillRecoCosMuonAngle(recoSummary.muonCosL);
    hist.reco.FillRecoShowerMultiplicity(recoSummary);

  }

}

// Helper method to build RecoInteractionSummary
RecoInteractionSummary RecoSelection::BuildRecoSummary(
    const caf::SRInteraction& dlpixn,
    const Mx2MatchResult& mx2MatchResult) const
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
  summary.muonCosL = mx2MatchResult.LArTrackDir.Dot(fBeam.beam_dir);
  const auto& mx2MatchLArTrack = dlpixn.part.dlp[mx2MatchResult.LArTrackIdx];
  summary.mx2MatchLArStartPosX = mx2MatchLArTrack.start.x;
  summary.mx2MatchLArStartPosY = mx2MatchLArTrack.start.y;
  summary.mx2MatchLArStartPosZ = mx2MatchLArTrack.start.z;

  //--------------------------------------------------
  // Loop over primaries
  //--------------------------------------------------

  for (const auto& part : dlpixn.part.dlp) {

    const int pdg = std::abs(part.pdg);
    //TODO: Add pi0 accounting

    //----------------------------------------------
    // Showers
    //----------------------------------------------

    if (abs(pdg) == 11) {
        if (part.primary == true)
            ++summary.nPrimElectrons;
        else
            ++summary.nSecElectrons;
    }

    if (abs(pdg) == 22) {
        if (part.primary == true)
            ++summary.nPrimPhotons;
        else
            ++summary.nSecPhotons;
    }
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
MatchedInteractionSummary RecoSelection::BuildMatchedIxnSummary(
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
    if(sr.mc.nu[summary.bestMatchIndex].id < fSelCuts.maxRockIxnTruthId) {
      summary.isRockIxn = true;
    }

    // Fill in truth summary for best-matched interaction if it exists
    if (summary.bestMatchIndex != -999) {
        const auto& bestTruthIxn = sr.mc.nu[summary.bestMatchIndex];
        TruthSelection truthSel(fSelCuts, fBeam, fDetector);
        summary.truthSummaryforBestMatch = truthSel.BuildTruthSummary(bestTruthIxn);
        
        // Get difference in truth and reco vertex positions for best-matched interaction
        summary.diffTruthRecoVertex = DiffPoints3D(summary.truthSummaryforBestMatch.vertex, dlpixn.vtx);
        
        summary.passesLArCuts = truthSel.IxnPassesTruthLArCuts(summary.truthSummaryforBestMatch);
        summary.passesMx2 = summary.truthSummaryforBestMatch.passesMx2;
    }

    return summary;
}

void RecoSelection::FillParticleTruthMatching(const caf::StandardRecord& sr,
                               const caf::SRInteraction& dlpixn,
                               MatchedInteractionSummary& matchSummary,
                               RecoInteractionSummary& recoSummary,
                               const Mx2MatchResult& mx2MatchResult)
{
    // --------------------------------------
    // Mx2 Track Matching Truth Information
    // --------------------------------------
    auto mx2TrackMatchP = caf::SRLorentzVector();
    auto mx2TrackLArStartPos = caf::SRVector3D();
    if (mx2MatchResult.truthIxnMx2PartType == 1) { // Primary track
        matchSummary.truthMatchMx2TrackisPrimary=true;
        matchSummary.truthMatchMx2TrackPDG = sr.mc.nu[mx2MatchResult.truthIxnMx2IxnIdx].prim[mx2MatchResult.truthIxnMx2PartIdx].pdg;
        std::cout<<"[DEBUG] Truth Match Track PDG "<<matchSummary.truthMatchMx2TrackPDG<<std::endl;
        mx2TrackMatchP = sr.mc.nu[mx2MatchResult.truthIxnMx2IxnIdx].prim[mx2MatchResult.truthIxnMx2PartIdx].p;
        mx2TrackLArStartPos = sr.mc.nu[mx2MatchResult.truthIxnMx2IxnIdx].prim[mx2MatchResult.truthIxnMx2PartIdx].start_pos;
    }
    if (mx2MatchResult.truthIxnMx2PartType == 3) { // Secondary track
        matchSummary.truthMatchMx2TrackisPrimary=false;
        matchSummary.truthMatchMx2TrackPDG = sr.mc.nu[mx2MatchResult.truthIxnMx2IxnIdx].sec[mx2MatchResult.truthIxnMx2PartIdx].pdg;
        mx2TrackMatchP = sr.mc.nu[mx2MatchResult.truthIxnMx2IxnIdx].sec[mx2MatchResult.truthIxnMx2PartIdx].p;
        mx2TrackLArStartPos = sr.mc.nu[mx2MatchResult.truthIxnMx2IxnIdx].sec[mx2MatchResult.truthIxnMx2PartIdx].start_pos;
    }
    // Fill energy/cosL info
    matchSummary.truthMatchMx2TrackE = mx2TrackMatchP.E;
    double mx2TrackMatchPMag = TMath::Sqrt(mx2TrackMatchP.px * mx2TrackMatchP.px + 
                                   mx2TrackMatchP.py * mx2TrackMatchP.py + 
                                   mx2TrackMatchP.pz * mx2TrackMatchP.pz);
    TVector3 mx2TrackMatchDir(mx2TrackMatchP.px / mx2TrackMatchPMag, 
                              mx2TrackMatchP.py / mx2TrackMatchPMag,
                              mx2TrackMatchP.pz / mx2TrackMatchPMag);
    matchSummary.truthMatchMx2TrackCosL = mx2TrackMatchDir.Dot(fBeam.beam_dir);
    // Fill Truth Match Start Points for LAr Track
    matchSummary.truthMatchMx2TrackLArStartPosX = mx2TrackLArStartPos.x;
    matchSummary.truthMatchMx2TrackLArStartPosY = mx2TrackLArStartPos.y;
    matchSummary.truthMatchMx2TrackLArStartPosZ = mx2TrackLArStartPos.z;

}

void RecoSelection::FillTruthMatchedCuts(const MatchedInteractionSummary& matchSummary,
                     CutFlowManager& cuts,
                     const std::string& recoCut)
{
    if (matchSummary.diffTruthRecoVertex >= fSelCuts.maxTruthRecoVertexDiff)
        return;

    cuts.Count("RecoMatchedValidVtx", recoCut);

    if (!matchSummary.passesLArCuts)
        return;

    cuts.Count("RecoMatchedSignal", recoCut);

    if (!matchSummary.passesMx2)
        return;

    cuts.Count("RecoMatchedSignalwithMx2", recoCut);
}

double RecoSelection::DiffPoints3D(const caf::SRVector3D& v1, const caf::SRVector3D& v2) 
{
    return TMath::Sqrt((v1.x - v2.x) * (v1.x - v2.x) +
                       (v1.y - v2.y) * (v1.y - v2.y) +
                       (v1.z - v2.z) * (v1.z - v2.z));
}