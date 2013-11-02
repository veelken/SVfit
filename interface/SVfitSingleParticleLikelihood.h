#ifndef TauAnalysis_SVfit_SVfitSingleParticleLikelihood_h
#define TauAnalysis_SVfit_SVfitSingleParticleLikelihood_h

/** \class SVfitSingleParticleLikelihoodBase
 *
 * Abstract base-class for plugins computing likelihood of single particle kinematics;
 * used by nSVfit algorithm
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.14 $
 *
 * $Id: SVfitSingleParticleLikelihood.h,v 1.14 2012/03/22 11:27:21 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "TauAnalysis/SVfit/interface/SVfitLikelihoodBase.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitSingleParticleHypothesis.h"
#include "AnalysisDataFormats/SVfit/interface/SVfitResonanceHypothesis.h"

#include <TF1.h>
#include <TFormula.h>
#include <TString.h>

#include <string>
#include <iostream>

class SVfitSingleParticleLikelihood : public SVfitLikelihoodBase
{
 public:
  SVfitSingleParticleLikelihood(const edm::ParameterSet& cfg)
    : SVfitLikelihoodBase(cfg),
      prodParticleLabel_(cfg.getParameter<std::string>("prodParticleLabel")),
      applyVisPtCutCorrection_(false),
      visPtCutThreshold_(0.)
  {
    if ( cfg.exists("applyVisPtCutCorrection") ) {
      applyVisPtCutCorrection_ = cfg.getParameter<bool>("applyVisPtCutCorrection");
      if ( applyVisPtCutCorrection_ ) {
	visPtCutThreshold_ = cfg.getParameter<double>("visPtCutThreshold");
	if ( visPtCutThreshold_ < 0. ) visPtCutThreshold_ = 0.;
      }
    }
  }

  virtual ~SVfitSingleParticleLikelihood() {}

  virtual void beginCandidate(const SVfitSingleParticleHypothesis*) {}

  virtual double operator()(const SVfitSingleParticleHypothesis*, int) const = 0;

 protected:
  std::string prodParticleLabel_;

  bool applyVisPtCutCorrection_;
  double visPtCutThreshold_;
};

#include "FWCore/PluginManager/interface/PluginFactory.h"

typedef edmplugin::PluginFactory<SVfitSingleParticleLikelihood* (const edm::ParameterSet&)> SVfitSingleParticleLikelihoodPluginFactory;

#endif
