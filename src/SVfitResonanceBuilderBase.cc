#include "TauAnalysis/SVfit/interface/SVfitResonanceBuilderBase.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitResonanceHypothesis.h"
#include "AnalysisDataFormats/SVfit/interface/SVfitSingleParticleHypothesis.h"

#include "TauAnalysis/SVfit/interface/svFitAuxFunctions.h"

SVfitResonanceBuilderBase::SVfitResonanceBuilderBase(const edm::ParameterSet& cfg)
  : SVfitBuilderBase(cfg),
    prodResonanceLabel_(cfg.getParameter<std::string>("prodResonanceLabel")),
    numDaughterBuilders_(0)
{
// nothing to be done yet... 
}

SVfitResonanceBuilderBase::~SVfitResonanceBuilderBase()
{
  for ( std::vector<SVfitSingleParticleBuilderBase*>::iterator it = daughterBuilders_.begin();
	it != daughterBuilders_.end(); ++it ) {
    delete (*it);
  }
}

void SVfitResonanceBuilderBase::beginJob(SVfitAlgorithmBase* algorithm)
{
  for ( std::vector<SVfitSingleParticleBuilderBase*>::iterator daughterBuilder = daughterBuilders_.begin();
	daughterBuilder != daughterBuilders_.end(); ++daughterBuilder ) {
    (*daughterBuilder)->beginJob(algorithm);
  }
}

void SVfitResonanceBuilderBase::beginEvent(const edm::Event& evt, const edm::EventSetup& es)
{
  SVfitBuilderBase::beginEvent(evt, es);
  for ( std::vector<SVfitSingleParticleBuilderBase*>::iterator daughterBuilder = daughterBuilders_.begin();
	daughterBuilder != daughterBuilders_.end(); ++daughterBuilder ) {
    (*daughterBuilder)->beginEvent(evt, es);
  }
}

SVfitResonanceHypothesis* SVfitResonanceBuilderBase::build(const inputParticleMap& inputParticles) const
{
  SVfitResonanceHypothesis* resonance = new SVfitResonanceHypothesis();

  reco::Candidate::LorentzVector p4(0,0,0,0);

  for ( std::vector<SVfitSingleParticleBuilderBase*>::const_iterator daughterBuilder = daughterBuilders_.begin();
	daughterBuilder != daughterBuilders_.end(); ++daughterBuilder ) {
    SVfitSingleParticleHypothesis* daughter = (*daughterBuilder)->build(inputParticles);
    daughter->setMother(resonance);

    p4 += daughter->p4();

    lastBuiltDaughters_[*daughterBuilder] = resonance->numDaughters();

    resonance->daughters_.push_back(daughter);
  }

  resonance->p4_ = p4;

  if ( resonance->numDaughters() == 2 ) {
    const SVfitSingleParticleHypothesis* daughter1 = resonance->daughter(0);
    const SVfitSingleParticleHypothesis* daughter2 = resonance->daughter(1);
    resonance->dPhiVis_ = TMath::ACos(TMath::Cos(daughter1->p4().phi() - daughter2->p4().phi()));
  }

  resonance->name_ = prodResonanceLabel_;

  resonance->barcode_ = barcodeCounter_;
  ++barcodeCounter_;

  return resonance;
}

void SVfitResonanceBuilderBase::finalize(SVfitResonanceHypothesis* resonance) const
{
  for ( std::vector<SVfitSingleParticleBuilderBase*>::const_iterator daughterBuilder = daughterBuilders_.begin();
	daughterBuilder != daughterBuilders_.end(); ++daughterBuilder ) {
    int idx = lastBuiltDaughters_[*daughterBuilder];
    SVfitSingleParticleHypothesis* daughter = dynamic_cast<SVfitSingleParticleHypothesis*>(&resonance->daughters_[idx]);
    assert(daughter);
    (*daughterBuilder)->finalize(daughter);
  }
}

bool SVfitResonanceBuilderBase::applyFitParameter(SVfitResonanceHypothesis* resonance, const double* params) const
{
  bool isValidSolution = true;

  for ( unsigned iDaughterBuilder = 0; iDaughterBuilder < numDaughterBuilders_; ++iDaughterBuilder ) {
    isValidSolution &= daughterBuilders_[iDaughterBuilder]->applyFitParameter(resonance->daughter(iDaughterBuilder), params);
  }

  reco::Candidate::LorentzVector dp4(0,0,0,0);

  size_t numDaughters = resonance->numDaughters();
  for ( size_t iDaughter = 0; iDaughter < numDaughters; ++iDaughter ) {
    const SVfitSingleParticleHypothesis* daughter = resonance->daughter(iDaughter);
    dp4 += daughter->dp4_fitted();
  }
  
  resonance->dp4_ = dp4;

  if ( numDaughters == 2 ) {
    const SVfitSingleParticleHypothesis* daughter1 = resonance->daughter(0);
    resonance->prod_angle_rf_ = svFit_namespace::gjAngleFromLabMomenta(resonance->p4_fitted(), daughter1->p4_fitted());
  }

  resonance->isValidSolution_ = isValidSolution;

  return isValidSolution;
}

void SVfitResonanceBuilderBase::print(std::ostream& stream) const
{
  stream << "<SVfitResonanceBuilderBase::print>:" << std::endl;
  stream << " pluginName = " << pluginName_ << std::endl;
  stream << " pluginType = " << pluginType_ << std::endl;
  stream << " prodResonanceLabel = " << prodResonanceLabel_ << std::endl;
  for ( std::vector<SVfitSingleParticleBuilderBase*>::const_iterator daughterBuilder = daughterBuilders_.begin();
	daughterBuilder != daughterBuilders_.end(); ++daughterBuilder ) {
    (*daughterBuilder)->print(stream);
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

EDM_REGISTER_PLUGINFACTORY(SVfitResonanceBuilderPluginFactory, "SVfitResonanceBuilderPluginFactory");


