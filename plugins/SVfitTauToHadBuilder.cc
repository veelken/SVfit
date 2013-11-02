
/** \class SVfitTauToHadBuilder
 *
 * Auxiliary class reconstructing tau --> had decays and
 * building SVfitTauToHadHypothesis objects;
 * used by SVfit algorithm
 *
 * \author Evan Friis, Christian Veelken; UC Davis
 *
 */

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TauAnalysis/SVfit/interface/SVfitTrackService.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

#include "TauAnalysis/SVfit/interface/SVfitTauDecayBuilder.h"
#include "TauAnalysis/SVfit/interface/SVfitSingleParticleTrackExtractor.h"
#include "TauAnalysis/SVfit/interface/SVfitParameter.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitTauToHadHypothesis.h"

using namespace svFit_namespace;

class SVfitTauToHadBuilder : public SVfitTauDecayBuilder
{
 public:
  SVfitTauToHadBuilder(const edm::ParameterSet& cfg)
    : SVfitTauDecayBuilder(cfg)
  {}
  
  void beginJob(SVfitAlgorithmBase* algorithm) 
  {
    SVfitTauDecayBuilder::beginJob(algorithm);

    //if ( verbosity_ ) print(std::cout);
  }

  SVfitSingleParticleHypothesis* build(const SVfitTauDecayBuilder::inputParticleMap& inputParticles) const 
  {
    inputParticleMap::const_iterator particlePtr = inputParticles.find(prodParticleLabel_);
    assert(particlePtr != inputParticles.end());

    SVfitTauToHadHypothesis* hypothesis = new SVfitTauToHadHypothesis(particlePtr->second, prodParticleLabel_, barcodeCounter_);
    if ( verbosity_ ) {
      std::cout << "<SVfitTauToHadBuilder::build>:" << std::endl;
      std::cout << " hypothesis #" << barcodeCounter_ << ": " << hypothesis << std::endl;
    }
    ++barcodeCounter_;

    SVfitTauDecayBuilder::initialize(hypothesis, particlePtr->second.get());

    return hypothesis;
  }

  // The neutrino system in hadronic decays is massless
  bool nuSystemIsMassless() const { return true; }

  virtual int getDecayMode(const reco::Candidate* candidate) const
  {
    const pat::Tau* tauPtr = dynamic_cast<const pat::Tau*>(candidate);
    assert(tauPtr);
    return tauPtr->decayMode();
  }

  virtual std::vector<const reco::Track*> extractTracks(const reco::Candidate* candidate) const
  {
    const pat::Tau* tauPtr = dynamic_cast<const pat::Tau*>(candidate);
    if ( !tauPtr ) 
      throw cms::Exception("SVfitTauToHadBuilder")
	<< "Failed to extract tracks. Please check if correct Builder plugin type is used !!\n";
    return trackExtractor_(*tauPtr);
  }

  void print(std::ostream& stream) const 
  {
    SVfitTauDecayBuilder::print(stream);
  }

 private:
  edm::Service<SVfitTrackService> trackService_;
  SVfitSingleParticleTrackExtractor<pat::Tau> trackExtractor_;
};

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(SVfitSingleParticleBuilderPluginFactory, SVfitTauToHadBuilder, "SVfitTauToHadBuilder");


