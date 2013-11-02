#ifndef TauAnalysis_SVfit_SVfitEventBuilder_h
#define TauAnalysis_SVfit_SVfitEventBuilder_h

/** \class SVfitEventBuilder
 *
 * Auxiliary class for building SVfitEventHypothesis objects;
 * used by SVfit algorithm
 *
 * \author Christian Veelken, UC Davis
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h" 

#include "DataFormats/VertexReco/interface/Vertex.h"

#include "TauAnalysis/SVfit/interface/SVfitEventBuilderBase.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitEventHypothesisBase.h"

class SVfitEventBuilder : public SVfitEventBuilderBase
{
 public:

  SVfitEventBuilder(const edm::ParameterSet&);
  ~SVfitEventBuilder() {}

  void beginJob(SVfitAlgorithmBase*);
  void beginEvent(const edm::Event&, const edm::EventSetup&);

  typedef edm::Ptr<reco::Candidate> CandidatePtr;
  typedef std::map<std::string, CandidatePtr> inputParticleMap;
  SVfitEventHypothesis* build(const inputParticleMap&, const reco::Vertex*) const;
  
 private:
  /// different possible polarization states of W bosons
  std::vector<int> polHandedness_;
  unsigned numPolStates_;

  SVfitAlgorithmBase* algorithm_;

  /// optional parameters for setting reconstructed to Monte Carlo truth values
  edm::InputTag srcGenVertex_;
  bool fixToGenVertex_;
  bool initializeToGenVertex_;
  mutable AlgebraicVector3 genVertexPos_;
};

#endif


