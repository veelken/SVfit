#ifndef TauAnalysis_SVfit_SVfitTauToLepLikelihoodMatrixElement_h
#define TauAnalysis_SVfit_SVfitTauToLepLikelihoodMatrixElement_h

/** \class SVfitTauToLepLikelihoodMatrixElement
 *
 * Plugin for computing likelihood for tau lepton decay 
 *  tau- --> e- nu nu (tau- --> mu- nu nu)
 * to be compatible with matrix element of V-A electroweak decay
 * (taking tau lepton polarization effects into account)
 *
 * \author Christian Veelken, UC Davis
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/SVfit/interface/SVfitSingleParticleLikelihood.h"

#include "AnalysisDataFormats/TauAnalysis/interface/SVfitSingleParticleHypothesisBase.h"

template <typename T>
class SVfitTauToLepLikelihoodMatrixElement : public SVfitSingleParticleLikelihood
{
 public:
  SVfitTauToLepLikelihoodMatrixElement(const edm::ParameterSet&);
  ~SVfitTauToLepLikelihoodMatrixElement();
  
  void beginJob(SVfitAlgorithmBase*);
  
  double operator()(const SVfitSingleParticleHypothesis*, int) const;

 private:
  bool applySinThetaFactor_;

  SVfitAlgorithmBase* algorithm_;
};

#endif
