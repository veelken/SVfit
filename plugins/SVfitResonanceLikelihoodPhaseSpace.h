#ifndef TauAnalysis_SVFit_plugins_SVfitResonanceLikelihoodPhaseSpace_h
#define TauAnalysis_SVFit_plugins_SVfitResonanceLikelihoodPhaseSpace_h

/** \class SVfitResonanceLikelihoodPhaseSpace
 *
 * Compute likelihood for resonance to decay into pairs of tau leptons,
 * assuming constant matrix element, so that energy and angular distribution 
 * of decay products are solely determined by phase-space
 *
 * \author Christian Veelken; LLR
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/SVfit/interface/SVfitResonanceLikelihood.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitResonanceHypothesis.h"

class SVfitResonanceLikelihoodPhaseSpace : public SVfitResonanceLikelihood
{
 public:

  SVfitResonanceLikelihoodPhaseSpace(const edm::ParameterSet&);
  ~SVfitResonanceLikelihoodPhaseSpace();

  void beginJob(SVfitAlgorithmBase*) const {}

  void beginCandidate(const SVfitResonanceHypothesis*) const {};

  double operator()(const SVfitResonanceHypothesis*, int) const;

 private:

  bool applySinThetaFactor_; 

  double power_;
};


#endif /* end of include guard: TauAnalysis_SVFit_plugins_SVfitResonanceLikelihoodPhaseSpace_h */
