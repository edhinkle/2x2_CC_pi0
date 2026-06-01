#pragma once

#include "config/ConfigLoader.h"
#include "cuts/Mx2MatchResult.h"
#include "duneanaobj/StandardRecord/StandardRecord.h"

class Mx2Matcher {
public:
  Mx2Matcher(const SelectionConfig& cfg,
             const BeamConfig& beam,
             const DetectorConfig& detector,
             const bool mcOnly):
             fSelCuts(cfg), fBeam(beam), fDetector(detector), fMCOnly(mcOnly) {};

  Mx2MatchResult MatchInteraction(
      const int LArIxnIdx,
      const caf::SRInteraction& dlpixn,
      const caf::StandardRecord& sr);

private:

  Mx2MatchResult MatchTrack(
      int LArTrackIdx,
      const int LArIxnIdx,
      const caf::SRInteraction& ixn,
      const caf::StandardRecord& sr);

  bool IsPrimaryTrack(const caf::SRRecoParticle& part, const caf::SRVector3D& vertex) const;
  bool IsVisibleTrack(const double length) const;
  bool IsExitingTrack(const caf::SRRecoParticle& part) const;

  double ComputeLength(const caf::SRRecoParticle& part) const;
  TVector3 ComputeDirection(const caf::SRRecoParticle& part) const;
  bool IsBackwardTrack(const caf::SRRecoParticle& part) const;
  bool Mx2TrackIsEnteringDownstream(const caf::SRTrack& trackMx2) const;
  bool Mx2TrackIsMuonCandidate(const caf::SRTrack& trackMx2) const;

  double GetNDCAFTrkMatchDotProd(
      const caf::StandardRecord& sr, 
      const caf::SRRecoParticle& LArPart,
      const int LArIxnIdx);

  bool PassesMx2NDCAFMakerCut(
      const caf::SRTrack& trackMx2, 
      const int LArIxnIdx,
      const caf::SRVector3D& LArStart,
      const caf::SRVector3D& LArEnd, 
      double &NDCAFMatchDotProd, 
      double &residual);

  const SelectionConfig& fSelCuts;
  const BeamConfig& fBeam;
  const DetectorConfig& fDetector;
  const bool fMCOnly;
};