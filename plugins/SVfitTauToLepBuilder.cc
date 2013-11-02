
/** \class SVfitTauToLepBuilder
 *
 * Auxiliary class reconstructing tau --> e/mu decays and
 * building SVfitTauToLepHypothesis objects;
 * used by SVfit algorithm
 *
 * \author Evan Friis, Christian Veelken; UC Davis
 *
 */

#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/TauReco/interface/PFTauDecayMode.h"

#include "TauAnalysis/CandidateTools/interface/SVfitTauDecayBuilder.h"
#include "TauAnalysis/CandidateTools/interface/SVfitSingleParticleTrackExtractor.h"
#include "TauAnalysis/CandidateTools/interface/SVfitParameter.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitTauToDaughtersHypothesisBaseT1T2.h"

template<typename T>
class SVfitTauToLepBuilder : public SVfitTauDecayBuilder 
{
 public:
  SVfitTauToLepBuilder(const edm::ParameterSet& cfg)
    : SVfitTauDecayBuilder(cfg)
  {}

  void beginJob(SVfitAlgorithmBase* algorithm) 
  {
    SVfitTauDecayBuilder::beginJob(algorithm);

    idxFitParameter_nuInvMass_ = getFitParameterIdx(algorithm, prodParticleLabel_, nSVfit_namespace::kTau_nuInvMass);

    //if ( verbosity_ ) print(std::cout);
  }

  SVfitSingleParticleHypothesis* build(const SVfitTauDecayBuilder::inputParticleMap& inputParticles) const 
  {
    inputParticleMap::const_iterator particlePtr = inputParticles.find(prodParticleLabel_);
    assert(particlePtr != inputParticles.end());

    SVfitTauToDaughtersHypothesisBaseT1T2<SVfitTauDecayHypothesis, T>* hypothesis = 
      new SVfitTauToDaughtersHypothesisBaseT1T2<SVfitTauDecayHypothesis, T>(particlePtr->second, prodParticleLabel_, barcodeCounter_);
    if ( verbosity_ ) {
      std::cout << "<SVfitTauToLepBuilder::build>:" << std::endl;
      std::cout << " hypothesis #" << barcodeCounter_ << ": " << hypothesis << std::endl;
    }
    ++barcodeCounter_;

    SVfitTauDecayBuilder::initialize(hypothesis, particlePtr->second.get());

    return hypothesis;
  }

  // The two neutrino system in leptonic tau decays can have non-zero mass
  bool nuSystemIsMassless() const { return false; }
  
  virtual int getDecayMode(const reco::Candidate* candidate) const 
  {
    assert(0); // force template specializations for pat::Electrons/pat::Muons to be used
  }

  virtual std::vector<const reco::Track*> extractTracks(const reco::Candidate* candidate) const 
  {
    const T* objPtr = dynamic_cast<const T*>(candidate);
    if ( !objPtr ) 
      throw cms::Exception("SVfitTauToLepBuilder")
	<< "Failed to extract tracks. Please check if correct Builder plugin type is used !!\n";
    return trackExtractor_(*objPtr);
  }

 private:
  SVfitSingleParticleTrackExtractor<T> trackExtractor_;
};

template<>
int SVfitTauToLepBuilder<pat::Electron>::getDecayMode(const reco::Candidate* candidate) const 
{
  return reco::PFTauDecayMode::tauDecaysElectron;
}

template <>
int SVfitTauToLepBuilder<pat::Muon>::getDecayMode(const reco::Candidate* candidate) const
{
  return reco::PFTauDecayMode::tauDecayMuon;
}

typedef SVfitTauToLepBuilder<pat::Electron> SVfitTauToElecBuilder;
typedef SVfitTauToLepBuilder<pat::Muon> SVfitTauToMuBuilder;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(SVfitSingleParticleBuilderPluginFactory, SVfitTauToElecBuilder, "SVfitTauToElecBuilder");
DEFINE_EDM_PLUGIN(SVfitSingleParticleBuilderPluginFactory, SVfitTauToMuBuilder, "SVfitTauToMuBuilder");

