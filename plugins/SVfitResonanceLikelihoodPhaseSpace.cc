#include "TauAnalysis/SVfit/plugins/SVfitResonanceLikelihoodPhaseSpace.h"

#include <TMath.h>

SVfitResonanceLikelihoodPhaseSpace::SVfitResonanceLikelihoodPhaseSpace(const edm::ParameterSet& cfg)
  : SVfitResonanceLikelihood(cfg)
{
  applySinThetaFactor_ = cfg.exists("applySinThetaFactor") ?
    cfg.getParameter<bool>("applySinThetaFactor") : false;

  power_ = cfg.getParameter<double>("power");
}

SVfitResonanceLikelihoodPhaseSpace::~SVfitResonanceLikelihoodPhaseSpace()
{
// nothing to be done yet...
}

double SVfitResonanceLikelihoodPhaseSpace::operator()(const SVfitResonanceHypothesis* resonance, int polHandedness) const 
{
  //if ( this->verbosity_ ) std::cout << "<SVfitResonanceLikelihoodPhaseSpace::operator()>:" << std::endl;

  double prodAngle_rf = resonance->prod_angle_rf();
  //if ( this->verbosity_ ) std::cout << " prodAngle_rf = " << prodAngle_rf << std::endl;

  double prob = 1.;
  if ( applySinThetaFactor_ ) prob *= (0.5*TMath::Sin(prodAngle_rf)); // phase-space factor 
				                                      // (to be used only in "fit", **not** in integration mode)
  
  //if ( this->verbosity_ ) std::cout << "--> prob = " << prob << std::endl;

  return prob;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(SVfitResonanceLikelihoodPluginFactory, SVfitResonanceLikelihoodPhaseSpace, "SVfitResonanceLikelihoodPhaseSpace");
