#ifndef TauAnalysis_SVfit_SVfitEventVertexRefitter_h
#define TauAnalysis_SVfit_SVfitEventVertexRefitter_h

/** \class SVfitEventVertexRefitter
 *
 * Class to refit position of primary event vertex,
 * excluding the tracks associated to tau decay products
 *
 * \author Evan Friis, Christian Veelken; UC Davis
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

class SVfitEventVertexRefitter
{
 public:
  SVfitEventVertexRefitter(const edm::ParameterSet&);
  ~SVfitEventVertexRefitter();

  void beginEvent(const edm::Event&, const edm::EventSetup&);

  /// Refit the primary vertex after subtracting the tau associated tracks.
  TransientVertex refit(const reco::Vertex*, const std::vector<const reco::Track*>* = 0) const;

 private:
  edm::InputTag srcBeamSpot_;

  const reco::BeamSpot* beamSpot_;
  VertexState beamSpotState_;
  bool beamSpotIsValid_;

  const TransientTrackBuilder* trackBuilder_;

  const VertexFitter<5>* vertexFitAlgorithm_;

  unsigned minNumTracksRefit_;

  bool applyBeamSpotConstraint_;

  int verbosity_;
};

#endif
