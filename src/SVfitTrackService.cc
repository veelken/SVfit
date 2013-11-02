#include "TauAnalysis/SVfit/interface/SVfitTrackService.h"

#include "FWCore/ServiceRegistry/interface/ActivityRegistry.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

SVfitTrackService::SVfitTrackService(const edm::ParameterSet& pset, edm::ActivityRegistry& reg)
  : builder_(NULL) 
{
  isValid_ = false;
  // Register the reset() function to be called after all modules have
  // processed the event.  This make sure the state is always consistent for
  // caller.
  reg.watchPostProcessEvent(this, &SVfitTrackService::reset);
}

void SVfitTrackService::setup(const edm::EventSetup& es, const reco::Candidate::Point& refPoint) 
{
  if ( !isValid_ ) {
    es.get<TransientTrackRecord>().get("TransientTrackBuilder", builder_);
    isValid_ = true;
    refPoint_ = AlgebraicVector3(refPoint.x(), refPoint.y(), refPoint.z());
  }
}

reco::TransientTrack
SVfitTrackService::transientTrack(const reco::Track* track) const 
{
  // Check if we've already made it.
  TransTrackCache::const_iterator lookup = cacheTransientTrack_.find(track);
  if ( lookup != cacheTransientTrack_.end() ) return lookup->second;

  // Build the transient track
  if ( !isValid_ || !builder_.isValid() ) {
    throw cms::Exception("NoTransTrackBuilder")
      << "<SVfitTrackService> The handle to the transient track builder is"
      << " not valid.  Please ensure that it is declared in the cfg file"
      << " and that setEventSetup(es) has been called on the track service"
      << " before use.";
  }
  reco::TransientTrack result = builder_->build(track);
  std::pair<TransTrackCache::const_iterator, bool> insertResult =
    cacheTransientTrack_.insert(std::make_pair(track, result));
  assert(insertResult.second);
  return insertResult.first->second;
}

void SVfitTrackService::reset(const edm::Event& evt, const edm::EventSetup&) 
{
  cacheTransientTrack_.clear();
  cacheTrackExtrapolations_.clear();
  isValid_ = false;
}
