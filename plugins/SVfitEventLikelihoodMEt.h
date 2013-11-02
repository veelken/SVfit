#ifndef TauAnalysis_SVfit_SVfitEventLikelihoodMEt_h
#define TauAnalysis_SVfit_SVfitEventLikelihoodMEt_h

/** \class SVfitEventLikelihoodMEt
 *
 * Plugin for computing likelihood for neutrinos produced in tau lepton decays
 * to match missing transverse momentum reconstructed in the event
 *
 * \author Evan Friis, Christian Veelken; UC Davis
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/SVfit/interface/SVfitEventLikelihood.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitEventHypothesis.h"

#include <TFormula.h>

class SVfitEventLikelihoodMEt : public SVfitEventLikelihood
{
 public:
  SVfitEventLikelihoodMEt(const edm::ParameterSet&);
  ~SVfitEventLikelihoodMEt();

  void beginJob(SVfitAlgorithmBase*) const;
  void beginCandidate(const SVfitEventHypothesis*) const;

  double operator()(const SVfitEventHypothesis*) const;

 private:
  double power_;

  TFormula* parSigma_;
  TFormula* parBias_;

  TFormula* perpSigma_;
  TFormula* perpBias_;

  mutable double qX_;
  mutable double qY_;
  mutable double qT_;
};

#endif
