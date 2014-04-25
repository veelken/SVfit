#include "TauAnalysis/SVfit/interface/SVfitAlgorithmBase.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TauAnalysis/SVfit/interface/svFitAuxFunctions.h"

#include <TMath.h>

using namespace svFit_namespace;

const SVfitAlgorithmBase* SVfitAlgorithmBase::gSVfitAlgorithm = 0;

SVfitAlgorithmBase::SVfitAlgorithmBase(const edm::ParameterSet& cfg)
  : currentEventHypothesis_(0),
    fitParameterCounter_(0)
{
  pluginName_ = cfg.getParameter<std::string>("pluginName");
  pluginType_ = cfg.getParameter<std::string>("pluginType");

  edm::ParameterSet cfgEvent = cfg.getParameter<edm::ParameterSet>("event");
  eventModel_ = new eventModelType(cfgEvent, allLikelihoods_);

  verbosity_ = cfg.exists("verbosity") ?
    cfg.getParameter<int>("verbosity") : 0;

  //std::cout << "<SVfitAlgorithmBase::SVfitAlgorithmBase (pluginName = " << pluginName_ << ")>:" << std::endl;
  //std::cout << " verbosity = " << verbosity_ << std::endl;
}

SVfitAlgorithmBase::~SVfitAlgorithmBase()
{
  delete eventModel_;
}

void SVfitAlgorithmBase::beginJob()
{
  for ( std::vector<SVfitLikelihoodBase*>::iterator likelihood = allLikelihoods_.begin();
	likelihood != allLikelihoods_.end(); ++likelihood ) {
    (*likelihood)->beginJob(this);
  }

  eventModel_->builder_->beginJob(this);
}

void SVfitAlgorithmBase::beginEvent(const edm::Event& evt, const edm::EventSetup& es)
{
  for ( std::vector<SVfitLikelihoodBase*>::iterator likelihood = allLikelihoods_.begin();
	likelihood != allLikelihoods_.end(); ++likelihood ) {
    (*likelihood)->beginEvent(evt, es);
  }

  eventModel_->builder_->beginEvent(evt, es);

  currentEventSetup_ = &es;
}

void SVfitAlgorithmBase::requestFitParameter(const std::string& name, int type, const std::string& requester)
{
  if ( name == "allTauDecays" ||
       name == "allLeptons"   ||
       name == "allNeutrinos" ) {
    //edm::LogWarning ("SVfitAlgorithmBase::requestFitParameter")
    //  << " Value = " << name << " not supported yet"
    //  << " --> relying on SingleParticleLikelihood plugins to initialize fitParameter for now.";
    return;
  }

  SVfitParameter* fitParameter = getFitParameter(name, type);

  if ( !fitParameter ) {
    assert(type >= 0 && type < svFit_namespace::kNumFitParameter);
    SVfitParameter newFitParameter(fitParameterCounter_, name, type);
    fitParameters_.push_back(newFitParameter);
    fitParameter = &fitParameters_.back();
    ++fitParameterCounter_;
  }

  fitParameter->regUsedBy(requester);
}

void SVfitAlgorithmBase::fixFitParameter(int fitParameterIdx)
{
  if ( fitParameterIdx != -1 ) {
    SVfitParameter* fitParameter = this->getFitParameter(fitParameterIdx);
    assert(fitParameter);
    fitParameter->setIsFixed(true);
  }
}

unsigned SVfitAlgorithmBase::getNumFitParameter(const std::string& name) const
{
  unsigned retVal = 0;

  for ( std::vector<SVfitParameter>::iterator fitParameter = fitParameters_.begin();
	fitParameter != fitParameters_.end(); ++fitParameter ) {
    if ( fitParameter->Name() == name ) ++retVal;
  }

  return retVal;
}

SVfitParameter* SVfitAlgorithmBase::getFitParameter(const std::string& name, int type) const
{
  SVfitParameter* retVal = 0;

  for ( std::vector<SVfitParameter>::iterator fitParameter = fitParameters_.begin();
	fitParameter != fitParameters_.end(); ++fitParameter ) {
    if ( fitParameter->Name() == name && fitParameter->Type() == type ) retVal = &(*fitParameter);
  }
  
  return retVal;
}

SVfitParameter* SVfitAlgorithmBase::getFitParameter(int idx) const
{
  assert(idx >= 0 && idx < (int)fitParameters_.size());

  SVfitParameter* retVal = &fitParameters_[idx];
  assert(retVal && retVal->index() == idx);

  return retVal;
}

void SVfitAlgorithmBase::setFitParameterInitialValue(int fitParameterIdx, double initialValue)
{
  if ( fitParameterIdx != -1 ) {
    SVfitParameter* fitParameter = this->getFitParameter(fitParameterIdx);
    assert(fitParameter);
    fitParameter->setInitialValue(initialValue);
  }
}

void SVfitAlgorithmBase::setFitParameterLimit(int fitParameterIdx, double lowerLimit, double upperLimit)
{
  if ( fitParameterIdx != -1 ) {
    SVfitParameter* fitParameter = this->getFitParameter(fitParameterIdx);
    assert(fitParameter);
    fitParameter->setLowerLimit(lowerLimit);
    fitParameter->setUpperLimit(upperLimit);
  }
}

void SVfitAlgorithmBase::setFitParameterStepSize(int fitParameterIdx, double stepSize)
{
  if ( fitParameterIdx != -1 ) {
    SVfitParameter* fitParameter = this->getFitParameter(fitParameterIdx);
    assert(fitParameter);
    fitParameter->setStepSize(stepSize);
  }
}

SVfitEventHypothesisBase* SVfitAlgorithmBase::fit(const inputParticleMap& inputParticles, const reco::Vertex* eventVertex) const
{
  // beginEvent should always be called before fit(...)
  assert(currentEventSetup_);

  // Setup the track service
  reco::Candidate::Point eventVertexPosition(0,0,0);
  if ( eventVertex ) eventVertexPosition = eventVertex->position();
  trackService_->setup(*currentEventSetup_, eventVertexPosition);

  currentEventHypothesis_ = eventModel_->builder_->build(inputParticles, eventVertex);
  currentEventHypothesis_->name_ = pluginName_;
  currentEventHypothesis_isValidSolution_ = true;

  eventModel_->beginCandidate(currentEventHypothesis_);

  gSVfitAlgorithm = this;

  if ( verbosity_ >= 1 ) {
    std::cout << "<SVfitAlgorithmBase::fit (pluginName = " << pluginName_ << ")>:" << std::endl;
    for ( std::vector<SVfitParameter>::const_iterator fitParameter = fitParameters_.begin();
	  fitParameter != fitParameters_.end(); ++fitParameter ) {
      fitParameter->dump(std::cout);
    }
  }

  fitImp();
  fittedEventHypothesis_->nll_ = fittedEventHypothesis_nll_;
  if ( verbosity_ >= 2 ) fittedEventHypothesis_->print(std::cout);

  return fittedEventHypothesis_;
}

bool SVfitAlgorithmBase::update(const double* x, const double* param) const
{
  currentEventHypothesis_isValidSolution_ = eventModel_->builder_->applyFitParameter(currentEventHypothesis_, x);
  if ( verbosity_ >= 2 ) {
    currentEventHypothesis_->print(std::cout);
    std::cout << "isValidSolution = " << currentEventHypothesis_isValidSolution_ << std::endl;
  }
  return currentEventHypothesis_isValidSolution_;
}

double SVfitAlgorithmBase::nll(const double* x, const double* param) const
{
  update(x, param);

  double nll = eventModel_->nll(currentEventHypothesis_);
  if ( TMath::IsNaN(nll) ) nll = std::numeric_limits<float>::max();

  return nll;
}

void SVfitAlgorithmBase::setMassResults(SVfitResonanceHypothesisBase* resonance, double value, double errUp, double errDown) const
{
  resonance->mass_ = value;
  resonance->massErrUp_ = errUp;
  resonance->massErrDown_ = errDown;
}

#include "FWCore/Framework/interface/MakerMacros.h"

EDM_REGISTER_PLUGINFACTORY(SVfitAlgorithmPluginFactory, "SVfitAlgorithmPluginFactory");


