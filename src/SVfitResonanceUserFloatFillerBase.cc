#include "TauAnalysis/SVfit/interface/SVfitResonanceUserFloatFillerBase.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitResonanceHypothesis.h"

#include "TauAnalysis/SVfit/interface/SVfitAlgorithmBase.h"
#include "TauAnalysis/SVfit/interface/svFitAuxFunctions.h"

using namespace svFit_namespace;

enum { kMax, kMedian, kMean3sigmaWithinMax, kMean5sigmaWithinMax };

SVfitResonanceUserFloatFillerBase::SVfitResonanceUserFloatFillerBase(const edm::ParameterSet& cfg)
  : pluginName_(cfg.getParameter<std::string>("pluginName")),
    pluginType_(cfg.getParameter<std::string>("pluginType"))
{
  std::string max_or_median_string = cfg.getParameter<std::string>("max_or_median");
  if      ( max_or_median_string == "max"                 ) max_or_median_ = kMax;
  else if ( max_or_median_string == "median"              ) max_or_median_ = kMedian;
  else if ( max_or_median_string == "mean3sigmaWithinMax" ) max_or_median_ = kMean3sigmaWithinMax;
  else if ( max_or_median_string == "mean5sigmaWithinMax" ) max_or_median_ = kMean5sigmaWithinMax;
  else throw cms::Exception("SVfitAlgorithmByIntegration2")
    << " Invalid Configuration Parameter 'max_or_median' = " << max_or_median_string << " !!\n";

  verbosity_ = cfg.exists("verbosity") ?
    cfg.getParameter<int>("verbosity") : 0;
}

SVfitResonanceUserFloatFillerBase::~SVfitResonanceUserFloatFillerBase()
{
// nothing to be done yet... 
}

/*
void SVfitResonanceUserFloatFillerBase::beginEvent(const edm::Event& evt, const edm::EventSetup& es)
{
// nothing to be done yet... 
}
 */

void SVfitResonanceUserFloatFillerBase::beginCandidate(const SVfitResonanceHypothesis*)
{
// nothing to be done yet... 
}

void SVfitResonanceUserFloatFillerBase::addUserFloat(SVfitResonanceHypothesis* resonance, const std::string& key, const TH1* histogram) const
{
  TH1* histogram_density = compHistogramDensity(histogram);

  double valueMaximum, valueMaximum_interpol, valueMean, valueQuantile016, valueQuantile050, valueQuantile084, valueMean3sigmaWithinMax, valueMean5sigmaWithinMax;
  extractHistogramProperties(
    histogram, histogram,
    valueMaximum, valueMaximum_interpol, valueMean, valueQuantile016, valueQuantile050, valueQuantile084, valueMean3sigmaWithinMax, valueMean5sigmaWithinMax,
    0);
  double value = 0.;
  if      ( max_or_median_ == kMax                 ) value = valueMaximum_interpol;
  else if ( max_or_median_ == kMedian              ) value = valueQuantile050;
  else if ( max_or_median_ == kMean3sigmaWithinMax ) value = valueMean3sigmaWithinMax;
  else if ( max_or_median_ == kMean5sigmaWithinMax ) value = valueMean5sigmaWithinMax;

  resonance->addUserFloat(key, value);

  delete histogram_density;
}

#include "FWCore/Framework/interface/MakerMacros.h"

EDM_REGISTER_PLUGINFACTORY(SVfitResonanceUserFloatFillerPluginFactory, "SVfitResonanceUserFloatFillerPluginFactory");


