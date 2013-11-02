#ifndef TauAnalysis_SVfit_plugins_SVfitResonanceLikelihoodBreitWigner_h
#define TauAnalysis_SVfit_plugins_SVfitResonanceLikelihoodBreitWigner_h

/** \class SVfitResonanceLikelihoodBreitWigner
 *
 * Compute likelihood for reconstructed mass of resonance 
 * to be compatible with nominal mass
 *
 * \author Christian Veelken; LLR
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/SVfit/interface/SVfitResonanceLikelihood.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitResonanceHypothesis.h"

class SVfitResonanceLikelihoodBreitWigner : public SVfitResonanceLikelihood
{
 public:

  SVfitResonanceLikelihoodBreitWigner(const edm::ParameterSet&);
  ~SVfitResonanceLikelihoodBreitWigner();

  void beginJob(SVfitAlgorithmBase*) const {}

  void beginCandidate(const SVfitResonanceHypothesis*) const {};

  double operator()(const SVfitResonanceHypothesis*, int) const;

 private:

  double resonance_mass_;
  double resonance_mass2_;
  double resonance_width_;
  double resonance_width2_;
  
  double gamma_;
  double k_;

  double power_;
};


#endif /* end of include guard: TauAnalysis_SVfit_plugins_SVfitResonanceLikelihoodBreitWigner_h */
