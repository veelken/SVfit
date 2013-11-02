#include "TauAnalysis/SVfit/interface/SVfitDecayVertexFitter.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <TMath.h>

SVfitDecayVertexFitter::SVfitDecayVertexFitter(const edm::ParameterSet& cfg)
  : trackBuilder_(0),
    vertexFitAlgorithm_(0)
{
  vertexFitAlgorithm_ = new KalmanVertexFitter(true);

  minNumTracksFit_ = ( cfg.exists("minNumTracksFit") ) ? 
    cfg.getParameter<unsigned>("minNumTracksFit") : 2;
}

SVfitDecayVertexFitter::~SVfitDecayVertexFitter()
{
  delete vertexFitAlgorithm_;
}

void SVfitDecayVertexFitter::beginEvent(const edm::Event& evt, const edm::EventSetup& es)
{
//--- get pointer to TransientTrackBuilder
  edm::ESHandle<TransientTrackBuilder> trackBuilderHandle;
  es.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilderHandle);
  trackBuilder_ = trackBuilderHandle.product();
  if ( !trackBuilder_ ) 
    throw cms::Exception("SVfitDecayVertexFitter::beginEvent")
      << " Failed to access TransientTrackBuilder !!\n";
}

TransientVertex SVfitDecayVertexFitter::fitSecondaryVertex(const std::vector<const reco::Track*>& tracks) const 
{
//-- return dummy vertex for one-prong tau decays
//   for which tau decay vertex cannot be determined

  if ( tracks.size() < minNumTracksFit_ ) return TransientVertex();
  assert(trackBuilder_);

//--- build transient tracks
  std::vector<reco::TransientTrack> track_trajectories;
  for ( std::vector<const reco::Track*>::const_iterator track = tracks.begin();
	track != tracks.end(); ++track ) {
    track_trajectories.push_back(trackBuilder_->build(*track));
  }
  
//--- fit tau decay vertex
  return vertexFitAlgorithm_->vertex(track_trajectories);
}
