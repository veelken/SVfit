#ifndef TauAnalysis_SVfit_SVfitEventLikelihoodTrackInfo_h
#define TauAnalysis_SVfit_SVfitEventLikelihoodTrackInfo_h

/** \class SVfitEventLikelihoodTrackInfo
 *
 * Plugin for computing likelihood for position of primary event vertex
 * refitted after excluding tracks of particles produced in tau lepton decays
 * to be compatible with reconstructed tau lepton decay kinematics
 *
 * \author Evan Friis, Christian Veelken; UC Davis
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/SVfit/interface/SVfitEventLikelihood.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitEventHypothesis.h"

class SVfitEventLikelihoodTrackInfo : public SVfitEventLikelihood
{
 public:
  SVfitEventLikelihoodTrackInfo(const edm::ParameterSet&);
  ~SVfitEventLikelihoodTrackInfo();

  void beginJob(SVfitAlgorithmBase*);

  double operator()(const SVfitEventHypothesis*) const;
};

#endif
