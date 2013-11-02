#ifndef TauAnalysis_SVfit_SVfitTauToLepLikelihoodPhaseSpace_h
#define TauAnalysis_SVfit_SVfitTauToLepLikelihoodPhaseSpace_h

/** \class SVfitTauToLeptonLikelihoodPhaseSpace
 *
 * Plugin to compute likelihood for electron (muon) to be compatible 
 * with tau --> e nu nu (tau --> mu nu nu) three-body decay,
 * assuming constant matrix element, so that energy and angular distribution 
 * of decay products are solely determined by phase-space
 * 
 * \author Evan Friis, Christian Veelken; UC Davis
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/SVfit/interface/SVfitSingleParticleLikelihood.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitSingleParticleHypothesisBase.h"

template <typename T>
class SVfitTauToLepLikelihoodPhaseSpace : public SVfitSingleParticleLikelihood
{
 public:
  SVfitTauToLepLikelihoodPhaseSpace(const edm::ParameterSet&);
  ~SVfitTauToLepLikelihoodPhaseSpace();

  void beginJob(SVfitAlgorithmBase*);

  double operator()(const SVfitSingleParticleHypothesis*, int) const;

 private:
  bool applySinThetaFactor_;

  SVfitAlgorithmBase* algorithm_;
};

#endif
