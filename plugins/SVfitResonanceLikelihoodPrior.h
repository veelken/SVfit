#ifndef TauAnalysis_SVfit_plugins_SVfitResonanceLikelihoodPrior_h
#define TauAnalysis_SVfit_plugins_SVfitResonanceLikelihoodPrior_h

/** \class SVfitResonanceLikelihoodPrior
 *
 * Correct tau decay likelihood for effect of visible Pt cuts
 *
 * \author Christian Veelken; LLR
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/SVfit/interface/SVfitResonanceLikelihood.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitResonanceHypothesis.h"

#include <TFormula.h>

#include <vector>

class SVfitResonanceLikelihoodPrior : public SVfitResonanceLikelihood
{
 public:

  SVfitResonanceLikelihoodPrior(const edm::ParameterSet&);
  ~SVfitResonanceLikelihoodPrior();

  void beginJob(SVfitAlgorithmBase*) const {}

  void beginCandidate(const SVfitResonanceHypothesis*) const {};

  double operator()(const SVfitResonanceHypothesis*, int polHandedness) const;

 private:

  TFormula* function_;
  double xMin_;
  double xMax_;
  int numParameter_;
  std::vector<double> parameter_;

  double power_;
};


#endif /* end of include guard: TauAnalysis_SVfit_plugins_SVfitResonanceLikelihoodPrior_h */
