#ifndef TauAnalysis_CandidateTools_NSVfitEventLikelihoodMEt_h
#define TauAnalysis_CandidateTools_NSVfitEventLikelihoodMEt_h

/** \class NSVfitEventLikelihoodMEt
 *
 * Plugin for computing likelihood for neutrinos produced in tau lepton decays
 * to match missing transverse momentum reconstructed in the event
 *
 * \author Evan Friis, Christian Veelken; UC Davis
 *
 * \version $Revision: 1.3 $
 *
 * $Id: NSVfitEventLikelihoodMEt.h,v 1.3 2011/04/15 16:54:19 veelken Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/CandidateTools/interface/NSVfitEventLikelihood.h"

#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitEventHypothesis.h"

#include <TFormula.h>

class NSVfitEventLikelihoodMEt : public NSVfitEventLikelihood
{
 public:
  NSVfitEventLikelihoodMEt(const edm::ParameterSet&);
  ~NSVfitEventLikelihoodMEt();

  void beginJob(NSVfitAlgorithmBase*) const;
  void beginCandidate(const NSVfitEventHypothesis*) const;

  double operator()(const NSVfitEventHypothesis*) const;

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
