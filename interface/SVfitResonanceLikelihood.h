#ifndef TauAnalysis_SVfit_SVfitResonanceLikelihoodBase_h
#define TauAnalysis_SVfit_SVfitResonanceLikelihoodBase_h

/** \class SVfitResonanceLikelihoodBase
 *
 * Abstract base-class for plugins computing likelihood of resonances;
 * used by nSVfit algorithm
 *
 * \author Christian Veelken, UC Davis
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/SVfit/interface/SVfitLikelihoodBase.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitResonanceHypothesis.h"

#include <string>
#include <iostream>

class SVfitResonanceLikelihood : public SVfitLikelihoodBase
{
 public:
  SVfitResonanceLikelihood(const edm::ParameterSet& cfg)
    : SVfitLikelihoodBase(cfg)
  {}
  virtual ~SVfitResonanceLikelihood() {}

  virtual void beginCandidate(const SVfitResonanceHypothesis*) const {}

  virtual double operator()(const SVfitResonanceHypothesis*, int) const = 0;
};

#include "FWCore/PluginManager/interface/PluginFactory.h"

typedef edmplugin::PluginFactory<SVfitResonanceLikelihood* (const edm::ParameterSet&)> SVfitResonanceLikelihoodPluginFactory;

#endif
