#include "TauAnalysis/SVfit/plugins/SVfitResonanceLikelihoodPolarization.h"

#include "TauAnalysis/SVfit/interface/svFitAuxFunctions.h"

#include <TMath.h>

#include <string>

using namespace SVfit_namespace;

SVfitResonanceLikelihoodPolarization::SVfitResonanceLikelihoodPolarization(const edm::ParameterSet& cfg)
  : SVfitResonanceLikelihood(cfg)
{
  initializePolProbFunction(SVfitResonanceHypothesis::kPolLR, "LR", cfg);
  initializePolProbFunction(SVfitResonanceHypothesis::kPolRL, "RL", cfg);
  initializePolProbFunction(SVfitResonanceHypothesis::kPolLL, "LL", cfg);
  initializePolProbFunction(SVfitResonanceHypothesis::kPolRR, "RR", cfg);
  
  power_ = cfg.getParameter<double>("power");
}

SVfitResonanceLikelihoodPolarization::~SVfitResonanceLikelihoodPolarization()
{
  for ( std::map<int, polProbFunctionType*>::iterator it = polProbFunctions_.begin();
	it != polProbFunctions_.end(); ++it ) {
    delete it->second;
  }
}

void SVfitResonanceLikelihoodPolarization::initializePolProbFunction(int polHandedness, 
							       const std::string& cfgPolProbFunctionName, const edm::ParameterSet& cfg)
{
  if ( cfg.exists(cfgPolProbFunctionName) ) {
    edm::ParameterSet cfgPolProbFunction = cfg.getParameter<edm::ParameterSet>(cfgPolProbFunctionName);
    polProbFunctionType* polProbFunction = new polProbFunctionType(pluginName_, cfgPolProbFunction);
    polProbFunctions_.insert(std::pair<int, polProbFunctionType*>(polHandedness, polProbFunction));
  }
}

double SVfitResonanceLikelihoodPolarization::operator()(const SVfitResonanceHypothesis* resonance, int polHandedness) const 
{
  //if ( verbosity_ ) {
  //  std::cout << "<SVfitResonanceLikelihoodPolarization::operator()>:" << std::endl;
  //  std::cout << " mass = " << resonance->p4_fitted().mass() << std::endl;
  //  std::string polHandedness_string = "undefined";
  //  if      ( polHandedness == SVfitResonanceHypothesis::kPolLR ) polHandedness_string = "LR";
  //  else if ( polHandedness == SVfitResonanceHypothesis::kPolRL ) polHandedness_string = "RL";
  //  else if ( polHandedness == SVfitResonanceHypothesis::kPolLL ) polHandedness_string = "LL";
  //  else if ( polHandedness == SVfitResonanceHypothesis::kPolRR ) polHandedness_string = "RR";
  //  std::cout << " polHandedness = " << polHandedness_string << std::endl;
  //}

  assert(resonance);

  double prob = 0.;
  std::map<int, polProbFunctionType*>::const_iterator polProbFunction = polProbFunctions_.find(polHandedness);
  if ( polProbFunction != polProbFunctions_.end() ) {
    double x = resonance->p4_fitted().mass();
    prob = polProbFunction->second->eval(x);
  }
  
  //if ( this->verbosity_ ) std::cout << "--> prob = " << prob << std::endl;

  return prob;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(SVfitResonanceLikelihoodPluginFactory, SVfitResonanceLikelihoodPolarization, "SVfitResonanceLikelihoodPolarization");
