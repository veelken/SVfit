#include "TauAnalysis/SVfit/plugins/SVfitResonanceLikelihoodRegularization.h"

#include <TMath.h>
#include <TString.h>

#include <string>

SVfitResonanceLikelihoodRegularization::SVfitResonanceLikelihoodRegularization(const edm::ParameterSet& cfg)
  : SVfitResonanceLikelihood(cfg)
{
  TString nll_formula_string = cfg.getParameter<std::string>("nll").data();
  nll_formula_string.ReplaceAll("mass", "x");
  nll_formula_string.ReplaceAll("pt",   "y");
  
  nll_formula_ = new TFormula(TString(pluginName_.data()).Append("_nll").Data(), nll_formula_string.Data());

  power_ = cfg.getParameter<double>("power");
}

SVfitResonanceLikelihoodRegularization::~SVfitResonanceLikelihoodRegularization()
{
  delete nll_formula_;
}

double SVfitResonanceLikelihoodRegularization::operator()(const SVfitResonanceHypothesis* resonance, int polHandedness) const 
{
  assert(resonance);

  double x = resonance->p4_fitted().mass();
  double y = resonance->p4_fitted().pt();

  double nll = nll_formula_->Eval(x, y);
  //std::cout << "<SVfitResonanceLikelihoodRegularization::operator()>:" << std::endl;
  //std::cout << " formula = " << nll_formula_->GetTitle() << ":" 
  //          << " mass = " << x << ", pt = " << y << " --> nll = " << nll << std::endl;

  double prob = TMath::Exp(-power_*nll);
  
  return prob;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(SVfitResonanceLikelihoodPluginFactory, SVfitResonanceLikelihoodRegularization, "SVfitResonanceLikelihoodRegularization");
