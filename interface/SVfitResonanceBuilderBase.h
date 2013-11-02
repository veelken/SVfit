#ifndef TauAnalysis_SVfit_SVfitResonanceBuilderBase_h
#define TauAnalysis_SVfit_SVfitResonanceBuilderBase_h

/** \class SVfitResonanceBuilderBase
 *
 * Base-class for building SVfitResonanceHypothesis objects;
 * used by SVfit algorithm
 *
 * \author Christian Veelken, UC Davis
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/SVfit/interface/SVfitBuilderBase.h"
#include "TauAnalysis/SVfit/interface/SVfitSingleParticleBuilderBase.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitResonanceHypothesis.h"
#include "AnalysisDataFormats/SVfit/interface/SVfitEventHypothesis.h"

#include <string>
#include <iostream>

class SVfitResonanceBuilderBase : public SVfitBuilderBase
{
 public:
  SVfitResonanceBuilderBase(const edm::ParameterSet&);
  virtual ~SVfitResonanceBuilderBase();

  virtual void beginJob(SVfitAlgorithmBase*);
  virtual void beginEvent(const edm::Event&, const edm::EventSetup&);

  typedef edm::Ptr<reco::Candidate> CandidatePtr;
  typedef std::map<std::string, CandidatePtr> inputParticleMap;
  virtual SVfitResonanceHypothesis* build(const inputParticleMap&) const;

  virtual void finalize(SVfitResonanceHypothesis*) const;

  virtual bool applyFitParameter(SVfitResonanceHypothesis*, const double*) const;

  virtual void print(std::ostream&) const;

 protected:
  std::string prodResonanceLabel_;

  std::vector<SVfitSingleParticleBuilderBase*> daughterBuilders_;
  unsigned numDaughterBuilders_;

  mutable std::map<SVfitSingleParticleBuilderBase*, int> lastBuiltDaughters_;
};

#include "FWCore/PluginManager/interface/PluginFactory.h"

typedef edmplugin::PluginFactory<SVfitResonanceBuilderBase* (const edm::ParameterSet&)> SVfitResonanceBuilderPluginFactory;

#endif
