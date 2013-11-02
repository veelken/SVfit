#ifndef TauAnalysis_SVfit_SVfitTauDecayLikelihoodDummy_h
#define TauAnalysis_SVfit_SVfitTauDecayLikelihoodDummy_h

/** \class SVfitTauDecayLikelihoodDummy
 *
 * Dummy plugin that does actually not compute any likelihood function
 * but requests the fit/integration parameters necessary to build hadronic/leptonic tau decays only.
 *
 * \author Christian Veelken, LLR
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/SVfit/interface/SVfitSingleParticleLikelihood.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitSingleParticleHypothesisBase.h"

#include "RooAbsPdf.h"
#include "RooRealVar.h"

template<typename T>
class SVfitTauDecayLikelihoodDummy : public SVfitSingleParticleLikelihood
{
 public:
  SVfitTauDecayLikelihoodDummy(const edm::ParameterSet&);
  ~SVfitTauDecayLikelihoodDummy();

  void beginJob(SVfitAlgorithmBase*);
  
  double operator()(const SVfitSingleParticleHypothesis*, int) const;
 private:
  SVfitAlgorithmBase* algorithm_;
};

#endif
