#ifndef TauAnalysis_SVfit_SVfitEventLikelihoodBase_h
#define TauAnalysis_SVfit_SVfitEventLikelihoodBase_h

/** \class SVfitEventLikelihoodBase
 *
 * Abstract base-class for plugins computing likelihood of combinations of resonances (an "event");
 * used by nSVfit algorithm
 *
 * \author Christian Veelken, UC Davis
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/SVfit/interface/SVfitLikelihoodBase.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitEventHypothesis.h"

#include <string>
#include <iostream>

class SVfitEventLikelihood : public SVfitLikelihoodBase
{
 public:
  SVfitEventLikelihood(const edm::ParameterSet& cfg)
    : SVfitLikelihoodBase(cfg)
  {}
  virtual ~SVfitEventLikelihood() {}

  virtual void beginCandidate(const SVfitEventHypothesis*) const {}

  virtual double operator()(const SVfitEventHypothesis*) const = 0;
};

#include "FWCore/PluginManager/interface/PluginFactory.h"

typedef edmplugin::PluginFactory<SVfitEventLikelihood* (const edm::ParameterSet&)> SVfitEventLikelihoodPluginFactory;

#endif
