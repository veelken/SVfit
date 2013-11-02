#ifndef TauAnalysis_SVfit_SVfitBuilderBase_h
#define TauAnalysis_SVfit_SVfitBuilderBase_h

/** \class SVfitBuilderBase
 *
 * Abstract base-class for all plugins building fit hypotheses;
 * used by SVfit algorithm
 *
 * \author Christian Veelken, UC Davis
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <string>
#include <iostream>

class SVfitAlgorithmBase;

class SVfitBuilderBase
{
 public:
  SVfitBuilderBase(const edm::ParameterSet& cfg)
    : pluginName_(cfg.getParameter<std::string>("pluginName")),
      pluginType_(cfg.getParameter<std::string>("pluginType")),
      barcodeCounter_(0)
  {
    verbosity_ = cfg.exists("verbosity") ?
      cfg.getParameter<int>("verbosity") : 0;
  }
  virtual ~SVfitBuilderBase() {}

  virtual void beginJob(SVfitAlgorithmBase*) {}
  virtual void beginEvent(const edm::Event&, const edm::EventSetup&) { barcodeCounter_ = 0; }

  virtual void print(std::ostream&) const {}

 protected:
  int getFitParameterIdx(SVfitAlgorithmBase*, const std::string&, int, bool = false);

  std::string pluginName_;
  std::string pluginType_;

  int verbosity_;

  mutable int barcodeCounter_;
};

#endif
