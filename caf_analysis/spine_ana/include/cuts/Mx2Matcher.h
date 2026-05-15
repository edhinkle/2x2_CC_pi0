#pragma once

#include "cuts/Mx2MatchResult.h"
#include "duneanaobj/StandardRecord/StandardRecord.h"

class Mx2Matcher {
public:
  Mx2Matcher(const config::SelectionConfig& cfg,
             const config::BeamConfig& beam,
             const config::DetectorConfig& detector,
             bool mcOnly);

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

  bool IsPrimaryTrack(const caf::SRParticle& part, const caf::SRVector3D& vertex) const;
  bool IsVisibleTrack(const double length) const;
  bool IsExitingTrack(const caf::SRParticle& part) const;

  double ComputeLength(const caf::SRParticle& part) const;
  TVector3 ComputeDirection(const caf::SRParticle& part) const;
  bool IsBackwardTrack(const caf::SRParticle& part) const;
  bool Mx2TrackIsEnteringDownstream(const caf::SRTrack& trackMx2) const;
  bool Mx2TrackIsMuonCandidate(const caf::SRTrack& trackMx2) const;

  double GetNDCAFTrkMatchDotProd(
      const caf::StandardRecord& sr, 
      const caf::SRParticle& LArPart,
      const int LArIxnIdx);

  bool PassesMx2NDCAFMakerCut(
      const caf::SRTrack& trackMx2, 
      const int LArIxnIdx,
      const caf::SRVector3D& LArStart,
      const caf::SRVector3D& LArEnd, 
      double &NDCAFMatchDotProd, 
      double &residual);

  const config::SelectionConfig& fSelCuts;
  const config::BeamConfig& fBeam;
  const config::DetectorConfig& fDetector;
  const bool fMCOnly;
};