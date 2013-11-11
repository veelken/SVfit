#ifndef TauAnalysis_SVfit_SVfitLeptonBuilderDummy_h
#define TauAnalysis_SVfit_SVfitLeptonBuilderDummy_h

/** \class SVfitLeptonBuilderDummy
 *
 * Auxiliary class to build SVfitSingleParticleHypothesis
 * representing electrons and muons which do not originate from tau decays;
 * used by SVfit algorithm
 *
 * \author Christian Veelken; LLR
 *
 */

#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/View.h"

#include "TauAnalysis/SVfit/interface/SVfitSingleParticleBuilderBase.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitSingleLeptonHypothesis.h"

template<typename T>
class SVfitLeptonBuilderDummy : public SVfitSingleParticleBuilderBase
{
 public:
  SVfitLeptonBuilderDummy(const edm::ParameterSet&);
  ~SVfitLeptonBuilderDummy() {}

  void beginEvent(const edm::Event&, const edm::EventSetup&);

  SVfitSingleParticleHypothesis* build(const SVfitSingleParticleBuilderBase::inputParticleMap& inputParticles) const;

  void finalize(SVfitSingleParticleHypothesis*) const {}

  bool applyFitParameter(SVfitSingleParticleHypothesis*, const double*) const { return true; } // CV: all solutions are valid

 private:
  void initialize(SVfitSingleLeptonHypothesis*, const reco::Candidate*) const;
    
  int getLeptonPdgId() const;

  /// optional parameters for setting reconstructed to Monte Carlo truth values
  edm::InputTag srcGenLeptons_;
  typedef edm::View<reco::GenParticle> GenParticleView;
  edm::Handle<GenParticleView> genParticles_;
  double dRmatch_;

  bool fixToGenP4_;

  mutable reco::Candidate::LorentzVector genP4_;
};

#endif /* end of include guard: TauAnalysis_SVfit_SVfitLeptonBuilderDummy_h */
