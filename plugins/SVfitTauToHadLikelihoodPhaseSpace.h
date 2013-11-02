#ifndef TauAnalysis_SVfit_SVfitTauToHadLikelihoodPhaseSpace_h
#define TauAnalysis_SVfit_SVfitTauToHadLikelihoodPhaseSpace_h

/** \class SVfitTauToLeptonLikelihoodPhaseSpace
 *
 * Plugin to compute likelihood for system of hadrons to be compatible 
 * with tau --> tau-jet + nu two-body decay,
 * assuming constant matrix element, so that energy and angular distribution 
 * of decay products are solely determined by phase-space
 * 
 * \author Evan Friis, Christian Veelken; UC Davis
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/SVfit/interface/SVfitSingleParticleLikelihood.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitSingleParticleHypothesisBase.h"

class SVfitTauToHadLikelihoodPhaseSpace : public SVfitSingleParticleLikelihood
{
 public:
  SVfitTauToHadLikelihoodPhaseSpace(const edm::ParameterSet&);
  ~SVfitTauToHadLikelihoodPhaseSpace();

  void beginJob(SVfitAlgorithmBase*);

  double operator()(const SVfitSingleParticleHypothesis*, int) const;

 private:
  bool applySinThetaFactor_; 

  bool applyVisMassFactor_;
  TH1* histogram_;
  int firstBin_;
  int lastBin_;

  SVfitAlgorithmBase* algorithm_;
};

#endif
