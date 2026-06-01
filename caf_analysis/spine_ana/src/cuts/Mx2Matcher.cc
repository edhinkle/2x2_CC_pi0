#include "cuts/Mx2Matcher.h"
#include "selection/RecoSelection.h"

// Constructor

Mx2MatchResult Mx2Matcher::MatchTrack(
    int LArTrackIdx,
    const int LArIxnIdx,
    const caf::SRInteraction& dlpixn,
    const caf::StandardRecord& sr)
{
  Mx2MatchResult result;
  result.LArTrackIdx = LArTrackIdx;

  const auto& part = dlpixn.part.dlp[LArTrackIdx];

  result.LArTrackDir = ComputeDirection(part);
  result.LArTrackLength = ComputeLength(part);
  result.trackIsBackward = IsBackwardTrack(part);

  // Check Track match in CAF track match branch (look at dot product between matches -- want close to 1)
  // Returns -999. if no match, in which case Mx2 matching will not be attempted for this track
  result.trkMatchDotProdCAF = GetNDCAFTrkMatchDotProd(sr, part, LArIxnIdx); 
  if (result.trkMatchDotProdCAF < 0.)
    return result;

  double bestMatchDotProd = -99.;

  for (size_t i = 0; i < sr.nd.minerva.ixn.size(); ++i) {

    const auto& mx2Ixn = sr.nd.minerva.ixn[i];

    for (size_t j = 0; j < mx2Ixn.ntracks; ++j) {

      const auto& mx2Track = mx2Ixn.tracks[j];
      result.mx2TrackEndZ = mx2Track.end.z;

      // Check if track passes Mx2 NDCAFMaker matching criteria (necessary due to offsets between Mx2 and 2x2 in data)
      caf::SRVector3D LArTrackStart;
      caf::SRVector3D LArTrackEnd;
      
      if (!result.trackIsBackward) {
        LArTrackStart = part.start;
        LArTrackEnd = part.end;
      }
      else {
        LArTrackStart = part.end;
        LArTrackEnd = part.start;
      }
      double NDCAFMatchDispl = 0; // only used as placeholder ptr for following method
      bool PassesMx2LArTrackMatch = PassesMx2NDCAFMakerCut(mx2Track, 
                                                           LArIxnIdx, 
                                                           LArTrackStart, 
                                                           LArTrackEnd, 
                                                           result.NDCAFMatchDotProd, 
                                                           NDCAFMatchDispl);
      if (!PassesMx2LArTrackMatch) continue;
      // Set Mx2LArTrackMatchDispl to abs. value to protect vs. direction differences w/ trkmatch
      result.NDCAFMatchDotProd = abs(result.NDCAFMatchDotProd);
      // TODO: Update when not MC only (because track match cuts)
      if (fMCOnly && abs(result.NDCAFMatchDotProd - result.trkMatchDotProdCAF)>fSelCuts.maxTrkMatchMx2MatchDiff)
        continue;

      // Require forward-going track into Mx2
      // TODO: Consider adding check to allow for vetoing interactions where
      //       there's a throughgoing muon (goes through Mx2 upstream and downstream)
      //       which could cause confusion with the muon from the neutrino
      if (!Mx2TrackIsEnteringDownstream(mx2Track)) continue;
      if (!Mx2TrackIsMuonCandidate(mx2Track)) continue; // check that it goes far enough in Mx2 DS in Z

      if (result.NDCAFMatchDotProd > bestMatchDotProd) {
        bestMatchDotProd = result.NDCAFMatchDotProd;

        result.mx2IxnIdx = i;
        result.mx2TrackIdx = j;

        // Optional truth
        if (fMCOnly) {
          result.truthIxnMx2PartIdx = mx2Track.truth[0].part;
          result.truthIxnMx2PartType = mx2Track.truth[0].type; // prim = 1, sec = 3
          result.truthIxnMx2IxnIdx  = mx2Track.truth[0].ixn;
        }
      }
    }
  }

  return result;
}

// Check if Track is a primary track of the interaction, starting near the vertex
bool Mx2Matcher::IsPrimaryTrack(const caf::SRRecoParticle& part, const caf::SRVector3D& vertex) const
{
    int pdg = part.pdg;

    if (!(abs(pdg) == 13 || abs(pdg) == 2212 || abs(pdg) == 211 || abs(pdg) == 321))
        return false;

    if (!(part.primary == true))
        return false;

    // Check that track starts near the vertex (to exclude secondaries from interactions in detector)
    double diffVertexStart = RecoSelection::DiffPoints3D(part.start, vertex);
    double diffVertexEnd = RecoSelection::DiffPoints3D(part.end, vertex);

    if (diffVertexStart < fSelCuts.maxMx2TrackVertexDiff || diffVertexEnd < fSelCuts.maxMx2TrackVertexDiff)
        return true;
    
    return false;
}

// Check if track is long enough to be visible in 2x2
bool Mx2Matcher::IsVisibleTrack(const double length) const
{
    // Require track to be over minimum track threshold for visibility in Mx2
    return (length > fDetector.minTrackLength);
}

// Check if track is exiting the detector (for Mx2 matching requirement)
bool Mx2Matcher::IsExitingTrack(const caf::SRRecoParticle& part) const
{
  return (abs(part.start.z) > fDetector.absFiducialZMaxExiting || abs(part.end.z) > fDetector.absFiducialZMaxExiting);
}

bool Mx2Matcher::IsBackwardTrack(const caf::SRRecoParticle& part) const
{
    // Check if track is backward-going based on start and end z positions
    return (part.end.z < part.start.z);
}

// Helper to check if Mx2 track enters downstream 
bool Mx2Matcher::Mx2TrackIsEnteringDownstream(const caf::SRTrack& trackMx2) const
{
    // Check if Mx2 track starts downstream of 2x2 within allowance region of 
    // start of downstream Mx2
    return (trackMx2.start.z > fDetector.absModuleZMax 
            && trackMx2.start.z < fSelCuts.maxMx2MatchTrackStartZ
            && trackMx2.end.z > fSelCuts.maxMx2MatchTrackStartZ);
}

// Helper to check if Mx2 Track travels far enough into Mx2 DS Z to be considered for muon candidate
bool Mx2Matcher::Mx2TrackIsMuonCandidate(const caf::SRTrack& trackMx2) const
{
    // Check if Mx2 track starts downstream of 2x2 within allowance region of 
    // start of downstream Mx2
    return (trackMx2.end.z > fSelCuts.minMx2MatchTrackEndZ);
}

// Helper method to compute length of track
double Mx2Matcher::ComputeLength(const caf::SRRecoParticle& part) const
{
  auto dX = part.end.x - part.start.x;
  auto dY = part.end.y - part.start.y;
  auto dZ = part.end.z - part.start.z;

  return std::sqrt(dX*dX + dY*dY + dZ*dZ);
}

// Helper method to compute direction of track
TVector3 Mx2Matcher::ComputeDirection(const caf::SRRecoParticle& part) const
{
  double L = ComputeLength(part);
  if (L <= 0) return TVector3(0,0,0);

  // Check direction using dZ (forward if positive, backward if negative)
  double dirZ = (part.end.z - part.start.z)/L;
  if (dirZ < 0) {
    // If backward-going, flip direction to point into Mx2
    return TVector3(
      (part.start.x - part.end.x)/L,
      (part.start.y - part.end.y)/L,
      (part.start.z - part.end.z)/L
    );
  }
  else {
    // If forward-going, direction is as-is
    return TVector3(
      (part.end.x - part.start.x)/L,
      (part.end.y - part.start.y)/L,
      (part.end.z - part.start.z)/L
    );
  }
}


// Helper method to get dot product of track extrapolation matching from trkmatch CAF branch
double Mx2Matcher::GetNDCAFTrkMatchDotProd(const caf::StandardRecord& sr, 
  const caf::SRRecoParticle& LArPart,
  const int LArIxnIdx)
{
  double trkMatchDotProdCAF=-999;
  for(size_t i=0; i<sr.nd.trkmatch.extrap.size(); i++)
  {
    // Check that you you are looking at the right interaction and that the track is reco (not truth)
    if (sr.nd.trkmatch.extrap[i].larid.ixn==LArIxnIdx && sr.nd.trkmatch.extrap[i].larid.reco==1)
      {
        int trkMatchIdx=sr.nd.trkmatch.extrap[i].larid.idx;
        // Check that the track matched in trkmatch is the same track as the one being considered in Mx2 matching (by comparing start and end z positions)
        if (sr.nd.lar.dlp[LArIxnIdx].tracks[trkMatchIdx].start.z==LArPart.start.z && sr.nd.lar.dlp[LArIxnIdx].tracks[trkMatchIdx].end.z==LArPart.end.z)
          {
            double cosAnglDispl=abs(sr.nd.trkmatch.extrap[i].angdispl);
            if (trkMatchDotProdCAF<cosAnglDispl) trkMatchDotProdCAF=cosAnglDispl;
          }
      }
  }
  return trkMatchDotProdCAF;
}


// Track Matching logic from NDCAFMaker to check if trkmatch branch logic is consistent
// Only change from CC numu inclusive selection code (from NDCAFMaker) is inputting LAr
// track start and end points as SRVector3D instead of individual points
bool Mx2Matcher::PassesMx2NDCAFMakerCut(
      const caf::SRTrack& trackMx2, 
      const int LArIxnIdx,
      const caf::SRVector3D& LArStart,
      const caf::SRVector3D& LArEnd, 
      double &NDCAFMatchDotProd, 
      double &residual)
{

    double z_extr = fSelCuts.mx2MatchExtrapLocZ;
    double d_x = fSelCuts.maxMx2MatchExtrapDiffX;
    double d_y = fSelCuts.maxMx2MatchExtrapDiffY;
    double d_thetax = fSelCuts.maxMx2MatchArctanDiffXZ;
    double d_thetay = fSelCuts.maxMx2MatchArctanDiffYZ;

    // Get LAr Track points
    double x1_lar = LArStart.x;
    double x2_lar = LArEnd.x;
    double y1_lar = LArStart.y;
    double y2_lar = LArEnd.y;
    double z1_lar = LArStart.z;
    double z2_lar = LArEnd.z;

    // Get Mx2 Track points
    double x1_minerva = trackMx2.start.x;
    double x2_minerva = trackMx2.end.x;
    double y1_minerva = trackMx2.start.y;
    double y2_minerva = trackMx2.end.y;
    double z1_minerva = trackMx2.start.z;
    double z2_minerva = trackMx2.end.z;

    double dX=x2_lar-x1_lar;
    double dY=y2_lar-y1_lar;
    double dZ=z2_lar-z1_lar;
    double track_LarLen=TMath::Sqrt(dX*dX+dY*dY+dZ*dZ);

    /*
    The experimental setup: Liquid Argon Detector is placed betbeen two MINERvA planes.
    To define matching criteria it is needed to find angles between LAr and MINERvA tracks.
    For LAr detector resolution is diffeent in X and Y direction, therefore it is needed to find angles between tracks
    as finction of the angle in X direction and as the function of an angle in Y direction. Distances between tracks
    will be calculated as distancec between extrapolated points - points of intersection of LAr and MINERvA tracls with
    the plane (parallel to plane XY) of LAr detector.
    */

    double tg_theta_mn_x = (x2_minerva - x1_minerva) / (z2_minerva - z1_minerva); // tangent of an angle between minerva track and X-axis
    double tg_theta_mn_y = (y2_minerva - y1_minerva) / (z2_minerva - z1_minerva); // tangent of an angle between minerva track and Y-axis
    double theta_mn_x = atan(tg_theta_mn_x);                                      // angle between minerva track and X-axis
    double theta_mn_y = atan(tg_theta_mn_y);                                      // angle between minerva track and Y-axis

    double tg_theta_nd_x = (x2_lar - x1_lar) / (z2_lar - z1_lar); // tangent of the angle between LAr track and X-axis
    double tg_theta_nd_y = (y2_lar - y1_lar) / (z2_lar - z1_lar); // tangent of the angle between LAr track and Y-axis
    double theta_nd_x = atan(tg_theta_nd_x);                      // angle between LAr track and X-axis
    double theta_nd_y = atan(tg_theta_nd_y);                      // angle between LAr track and Y-axis

    double delta_theta_x = theta_mn_y - theta_nd_y;
    double delta_theta_y = theta_mn_x - theta_nd_x;

    // Extrapolating Both tracks to the same point z = zextr (here it's the front of Lar)
    double t_mn = (z_extr - z1_minerva) / (z2_minerva - z1_minerva);
    double x_mn = t_mn * (x2_minerva - x1_minerva) + x1_minerva; // X-coordinate of extrapolated point of LAr track
    double y_mn = t_mn * (y2_minerva - y1_minerva) + y1_minerva; // Y-coordinate of extrapolated point of LAr track

    double t_nd = (z_extr - z1_lar) / (z2_lar - z1_lar); // parametr of the equation of the line (LAr track)
    double x_nd = t_nd * (x2_lar - x1_lar) + x1_lar;     // X-coordinate of extrapolated point of LAr track
    double y_nd = t_nd * (y2_lar - y1_lar) + y1_lar;     // Y-coordinate of extrapolated point of LAr track

    double dist_x = (x_mn - x_nd); // distance between X-coordinates of extrapolated points of minerva and LAr tracks
    double dist_y = (y_mn - y_nd); // distance between Y-coordinates of extrapolated points of minerva and LAr tracks

    residual = sqrt(pow(dist_x, 2) + pow(dist_y, 2));
    NDCAFMatchDotProd = ((x2_minerva - x1_minerva) * (x2_lar - x1_lar) +
                (y2_minerva - y1_minerva) * (y2_lar - y1_lar) +
                (z2_minerva - z1_minerva) * (z2_lar - z1_lar)) /
               (track_LarLen * trackMx2.len_cm); // angle between minerva and Lar tracks

    return (abs(delta_theta_x) < d_thetax && abs(delta_theta_y) < d_thetay && abs(dist_x) < d_x && abs(dist_y) < d_y);
}

// Main method to check if interaction has a track that matches Mx2 cuts and compute matching score
Mx2MatchResult Mx2Matcher::MatchInteraction(
    const int LArIxnIdx,
    const caf::SRInteraction& dlpixn,
    const caf::StandardRecord& sr)
{
  Mx2MatchResult bestMatch;

  for (size_t i = 0; i < dlpixn.part.dlp.size(); ++i) {

    const auto& part = dlpixn.part.dlp[i];

    if (!IsPrimaryTrack(part, dlpixn.vtx)) continue;
    if (!IsVisibleTrack(ComputeLength(part))) continue;
    if (!IsExitingTrack(part)) continue;

    auto result = MatchTrack(i, LArIxnIdx, dlpixn, sr);

    if (result.NDCAFMatchDotProd > bestMatch.NDCAFMatchDotProd) {
      bestMatch = result;
    }
  }

  bestMatch.isGoodMatch = (bestMatch.NDCAFMatchDotProd > fSelCuts.mx2TrackMatchDotProdThreshold);
  return bestMatch;
}