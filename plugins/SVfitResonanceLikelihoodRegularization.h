#ifndef TauAnalysis_SVfit_plugins_SVfitResonanceLikelihoodMassPenalty_h
#define TauAnalysis_SVfit_plugins_SVfitResonanceLikelihoodRegularization_h

/** \class SVfitResonanceLikelihoodRegularization
 *
 * Adds a penalty term for high masses.  
 * The return value of operator() is configurable via python.
 *
 * Options used in the past are:
 *  o TMath::Log(mass)
 *  o TMath::Log(TMath::Max([0], [1]*([2] - TMath::Erf((x - [3])*[4]))))
 *    with p0 = 5.00e-3, p1 = 4.21e-2, p2 = 2.52e-2, p3 = 4.40e+1, p4 = -6.90e-3 
 *   (efficiency of gg --> Higgs --> mu + tau_had channel in 2010 analysis)
 *  o TMath::Log(TMath::Max([0], [1]*([2] - TMath::Erf((x - [3])*[4]))))
 *    with p0 = 2.50e-3, p1 = 2.49e-2, p2 = 7.78e-2, p3 = 5.63e+1, p4 = -7.53e-3 
 *   (efficiency of gg --> Higgs --> e + tau_had channel in 2010 analysis)
 #  o TMath::Log(pt)
 *
 * \author Evan Friis, Christian Veelken; UC Davis
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/SVfit/interface/SVfitResonanceLikelihood.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitResonanceHypothesis.h"

#include "TFormula.h"

class SVfitResonanceLikelihoodRegularization : public SVfitResonanceLikelihood
{
 public:

  SVfitResonanceLikelihoodRegularization(const edm::ParameterSet&);
  ~SVfitResonanceLikelihoodRegularization();

  void beginJob(SVfitAlgorithmBase*) const {}

  void beginCandidate(const SVfitResonanceHypothesis*) const {};

  double operator()(const SVfitResonanceHypothesis*, int) const;

 private:

  TFormula* nll_formula_;
  
  double power_;
};


#endif /* end of include guard: TauAnalysis_SVfit_plugins_SVfitResonanceLikelihoodRegularization2_h */
