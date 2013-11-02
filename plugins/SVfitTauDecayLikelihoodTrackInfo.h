#ifndef TauAnalysis_SVfit_SVfitTauDecayLikelihoodTrackInfo_h
#define TauAnalysis_SVfit_SVfitTauDecayLikelihoodTrackInfo_h

/** \class SVfitTauDecayLikelihoodTrackInfo
 *
 * Plugin to compute likelihood for tracks of tau lepton decay products
 * to be compatible with originating from hypothetic secondary (tau lepton decay) vertex
 *
 * \author Evan Friis, Christian Veelken; UC Davis
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/SVfit/interface/SVfitSingleParticleLikelihood.h"

#include "AnalysisDataFormats/TauAnalysis/interface/SVfitSingleParticleHypothesisBase.h"

#include "TauAnalysis/SVfit/interface/SVfitLegTrackExtractor.h"
#include "TauAnalysis/SVfit/interface/SVfitTrackExtrapolation.h"

class SVfitTauDecayLikelihoodTrackInfo : public SVfitSingleParticleLikelihood
{
 public:
  SVfitTauDecayLikelihoodTrackInfo(const edm::ParameterSet&);
  ~SVfitTauDecayLikelihoodTrackInfo();

  void beginJob(SVfitAlgorithmBase*);
  void beginCandidate(const SVfitSingleParticleHypothesis*);

  double operator()(const SVfitSingleParticleHypothesis*, int) const;

 private:
  bool useLifetimeConstraint_;
  bool useBackwardsPenaltyTerm_;

  bool ignore3Prongs_;
  bool ignore1Prongs_;

  double sfProdVertexCov_;
  double sfDecayVertexCov_;

  SVfitAlgorithmBase* algorithm_;

  double decayDistance_lab_shift_lowerLimit_;
  double decayDistance_lab_shift_upperLimit_;
};

#endif
