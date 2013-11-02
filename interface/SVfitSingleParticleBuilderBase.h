#ifndef TauAnalysis_SVfit_SVfitSingleParticleBuilderBase_h
#define TauAnalysis_SVfit_SVfitSingleParticleBuilderBase_h

/** \class SVfitSingleParticleBuilderBase
 *
 * Base-class for building objects derrived from SVfitSingleParticleHypothesis class;
 * used by SVfit algorithm
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.8 $
 *
 * $Id: SVfitSingleParticleBuilderBase.h,v 1.8 2012/08/28 15:00:20 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/SVfit/interface/SVfitBuilderBase.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitSingleParticleHypothesis.h"
#include "AnalysisDataFormats/SVfit/interface/SVfitResonanceHypothesis.h"

#include <string>
#include <iostream>

class SVfitResonanceHypothesis;

class SVfitSingleParticleBuilderBase : public SVfitBuilderBase
{
 public:
  SVfitSingleParticleBuilderBase(const edm::ParameterSet& cfg)
    : SVfitBuilderBase(cfg),
      prodParticleLabel_(cfg.getParameter<std::string>("prodParticleLabel"))
  {}
  virtual ~SVfitSingleParticleBuilderBase() {}

  typedef edm::Ptr<reco::Candidate> CandidatePtr;
  typedef std::map<std::string, CandidatePtr> inputParticleMap;
  virtual SVfitSingleParticleHypothesis* build(const inputParticleMap&) const = 0;

  virtual void finalize(SVfitSingleParticleHypothesis*) const = 0;

  virtual bool applyFitParameter(SVfitSingleParticleHypothesis*, const double*) const = 0;

  virtual void print(std::ostream&) const {}

protected:
  std::string prodParticleLabel_;
};

#include "FWCore/PluginManager/interface/PluginFactory.h"

typedef edmplugin::PluginFactory<SVfitSingleParticleBuilderBase* (const edm::ParameterSet&)> SVfitSingleParticleBuilderPluginFactory;

#endif



