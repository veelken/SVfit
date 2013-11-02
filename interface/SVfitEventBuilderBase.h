#ifndef TauAnalysis_SVfit_SVfitEventBuilderBase_h
#define TauAnalysis_SVfit_SVfitEventBuilderBase_h

/** \class SVfitEventBuilderBase
 *
 * Base-class for building SVfitEventHypothesis objects;
 * used by SVfit algorithm
 *
 * \author Christian Veelken, UC Davis
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h" 

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/Ptr.h"

#include "TauAnalysis/SVfit/interface/SVfitBuilderBase.h"
#include "TauAnalysis/SVfit/interface/SVfitResonanceBuilderBase.h"
#include "TauAnalysis/SVfit/interface/SVfitEventVertexRefitter.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitEventHypothesis.h"

#include <string>
#include <iostream>

class SVfitEventBuilderBase : public SVfitBuilderBase
{
 public:
  SVfitEventBuilderBase(const edm::ParameterSet&);
  virtual ~SVfitEventBuilderBase();

  virtual void beginJob(SVfitAlgorithmBase*);
  virtual void beginEvent(const edm::Event&, const edm::EventSetup&);

  typedef edm::Ptr<reco::Candidate> CandidatePtr;
  typedef std::map<std::string, CandidatePtr> inputParticleMap;
  virtual SVfitEventHypothesis* build(const inputParticleMap&, const reco::Vertex*) const;

  virtual bool applyFitParameter(SVfitEventHypothesis*, const double*) const;

  virtual void print(std::ostream&) const;

 protected:
  std::vector<SVfitResonanceBuilderBase*> resonanceBuilders_;
  unsigned numResonanceBuilders_;

  mutable std::map<SVfitResonanceBuilderBase*, int> lastBuiltResonances_;

  edm::InputTag srcBeamSpot_;

  reco::BeamSpot beamSpot_;

  SVfitEventVertexRefitter* eventVertexRefitAlgorithm_;
  bool doEventVertexRefit_;

  int idxFitParameter_pvShiftX_;
  int idxFitParameter_pvShiftY_;
  int idxFitParameter_pvShiftZ_;
  bool doFitParameter_pvShift_;
};

#include "FWCore/PluginManager/interface/PluginFactory.h"

typedef edmplugin::PluginFactory<SVfitEventBuilderBase* (const edm::ParameterSet&)> SVfitEventBuilderPluginFactory;

#endif
