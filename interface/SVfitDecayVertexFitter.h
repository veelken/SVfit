#ifndef TauAnalysis_SVfit_SVfitDecayVertexFitter_h
#define TauAnalysis_SVfit_SVfitDecayVertexFitter_h

/** \class SVfitDecayVertexFitter
 *
 * Class to fit position of tau decay vertex
 * (speficic to hadronic 3-prong tau decays)
 *
 * \author Evan Friis, Christian Veelken; UC Davis
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

class SVfitDecayVertexFitter
{
 public:
  SVfitDecayVertexFitter(const edm::ParameterSet&);
  ~SVfitDecayVertexFitter();

  void beginEvent(const edm::Event&, const edm::EventSetup&);

  /// Fit decay vertex of three-prong tau
  TransientVertex fitSecondaryVertex(const std::vector<const reco::Track*>&) const;

 private:
  const TransientTrackBuilder* trackBuilder_;

  const KalmanVertexFitter* vertexFitAlgorithm_;

  unsigned minNumTracksFit_;
};

#endif
