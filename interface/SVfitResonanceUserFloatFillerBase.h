#ifndef TauAnalysis_SVfit_SVfitResonanceUserFloatFillerBase_h
#define TauAnalysis_SVfit_SVfitResonanceUserFloatFillerBase_h

/** \class SVfitResonanceUserFloatFillerBase
 *
 * Base-class for computing "user float" variables
 * added to SVfitResonanceHypothesis objects (like in the PAT objects)
 *
 * \author Christian Veelken, LLR/Ecole Polytechnique
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitResonanceHypothesis.h"

#include <TH1.h>

#include <string>

class SVfitAlgorithmBase;

class SVfitResonanceUserFloatFillerBase
{
 public:
  SVfitResonanceUserFloatFillerBase(const edm::ParameterSet&);
  virtual ~SVfitResonanceUserFloatFillerBase();

  //virtual void beginEvent(const edm::Event&, const edm::EventSetup&);
  virtual void beginCandidate(const SVfitResonanceHypothesis*);

  virtual void resetHistograms() = 0;

  virtual void fillHistograms(const SVfitResonanceHypothesis*) = 0;

  virtual void addUserFloats(SVfitResonanceHypothesis*) const = 0;

 protected:
  virtual void addUserFloat(SVfitResonanceHypothesis*, const std::string&, const TH1*) const;

  std::string pluginName_;
  std::string pluginType_;

  int max_or_median_;

  int verbosity_;
};

#include "FWCore/PluginManager/interface/PluginFactory.h"

typedef edmplugin::PluginFactory<SVfitResonanceUserFloatFillerBase* (const edm::ParameterSet&)> SVfitResonanceUserFloatFillerPluginFactory;

#endif
