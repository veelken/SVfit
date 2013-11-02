#include "TauAnalysis/SVfit/plugins/SVfitTauDecayLikelihoodDummy.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include "TauAnalysis/SVfit/interface/SVfitAlgorithmBase.h"
#include "TauAnalysis/SVfit/interface/SVfitParameter.h"
#include "TauAnalysis/SVfit/interface/svFitAuxFunctions.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitTauDecayHypothesis.h"

using namespace svFit_namespace;

template <typename T>
SVfitTauDecayLikelihoodDummy<T>::SVfitTauDecayLikelihoodDummy(const edm::ParameterSet& cfg)
  : SVfitSingleParticleLikelihood(cfg),
    algorithm_(0)
{
  if ( this->verbosity_ ) std::cout << "<SVfitTauDecayLikelihoodDummy::SVfitTauDecayLikelihoodDummy>:" << std::endl;
}

template <typename T>
SVfitTauDecayLikelihoodDummy<T>::~SVfitTauDecayLikelihoodDummy()
{
// nothing to be done yet...
}

template <typename T>
void SVfitTauDecayLikelihoodDummy<T>::beginJob(SVfitAlgorithmBase* algorithm)
{
  assert(0); // force template specializations for pat::Electrons/pat::Muons/pat::Taus to be used
}

template <>
void SVfitTauDecayLikelihoodDummy<pat::Electron>::beginJob(SVfitAlgorithmBase* algorithm)
{
  algorithm_ = algorithm;
  
  algorithm->requestFitParameter(prodParticleLabel_, nSVfit_namespace::kTau_visEnFracX, pluginName_);
  algorithm->requestFitParameter(prodParticleLabel_, nSVfit_namespace::kTau_phi_lab,    pluginName_);
  algorithm->requestFitParameter(prodParticleLabel_, nSVfit_namespace::kTau_nuInvMass,  pluginName_);
}

template <>
void SVfitTauDecayLikelihoodDummy<pat::Muon>::beginJob(SVfitAlgorithmBase* algorithm)
{
  algorithm_ = algorithm;

  algorithm->requestFitParameter(prodParticleLabel_, nSVfit_namespace::kTau_visEnFracX, pluginName_);
  algorithm->requestFitParameter(prodParticleLabel_, nSVfit_namespace::kTau_phi_lab,    pluginName_);
  algorithm->requestFitParameter(prodParticleLabel_, nSVfit_namespace::kTau_nuInvMass,  pluginName_);
}

template <>
void SVfitTauDecayLikelihoodDummy<pat::Tau>::beginJob(SVfitAlgorithmBase* algorithm)
{
  algorithm_ = algorithm;

  algorithm->requestFitParameter(prodParticleLabel_, nSVfit_namespace::kTau_visEnFracX, pluginName_);
  algorithm->requestFitParameter(prodParticleLabel_, nSVfit_namespace::kTau_phi_lab,    pluginName_);
}

//
//-------------------------------------------------------------------------------
//

template <typename T>
double SVfitTauDecayLikelihoodDummy<T>::operator()(const SVfitSingleParticleHypothesis* hypothesis, int polSign) const
{
  const SVfitTauDecayHypothesis* hypothesis_T = dynamic_cast<const SVfitTauDecayHypothesis*>(hypothesis);
  assert(hypothesis_T != 0);

  double visEnFracX = hypothesis_T->visEnFracX();

  //if ( this->verbosity_ ) {
  //  std::cout << "<SVfitTauDecayLikelihoodDummy::operator()>:" << std::endl;
  //  std::cout << " visEnFracX = " << visEnFracX << std::endl;
  //}

  double prob = 1.0;

  if ( algorithm_->applyJacobiFactors() && visEnFracX > 0. ) {
    double jacobiFactor = 1./(cube(visEnFracX));
    //if ( this->verbosity_ ) std::cout << "jacobiFactor = " << jacobiFactor << std::endl;
    prob *= jacobiFactor;
  }

  //if ( this->verbosity_ ) std::cout << "--> prob = " << prob << std::endl;

  return prob;
}

typedef SVfitTauDecayLikelihoodDummy<pat::Electron> SVfitTauToElecLikelihoodDummy;
typedef SVfitTauDecayLikelihoodDummy<pat::Muon> SVfitTauToMuLikelihoodDummy;
typedef SVfitTauDecayLikelihoodDummy<pat::Tau> SVfitTauToHadLikelihoodDummy;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(SVfitSingleParticleLikelihoodPluginFactory, SVfitTauToElecLikelihoodDummy, "SVfitTauToElecLikelihoodDummy");
DEFINE_EDM_PLUGIN(SVfitSingleParticleLikelihoodPluginFactory, SVfitTauToMuLikelihoodDummy, "SVfitTauToMuLikelihoodDummy");
DEFINE_EDM_PLUGIN(SVfitSingleParticleLikelihoodPluginFactory, SVfitTauToHadLikelihoodDummy, "SVfitTauToHadLikelihoodDummy");

