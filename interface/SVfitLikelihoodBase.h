#ifndef TauAnalysis_SVfit_SVfitLikelihoodBase_h
#define TauAnalysis_SVfit_SVfitLikelihoodBase_h

/** \class SVfitSingleParticleLikelihoodBase
 *
 * Abstract base-class for all likelihood function plugins;
 * used by nSVfit algorithm
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.4 $
 *
 * $Id: SVfitLikelihoodBase.h,v 1.4 2011/03/06 11:31:11 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <string>
#include <iostream>

class SVfitAlgorithmBase;

class SVfitLikelihoodBase
{
 public:
  SVfitLikelihoodBase(const edm::ParameterSet& cfg)
  {
    pluginType_ = cfg.getParameter<std::string>("pluginType");
    pluginName_ = cfg.getParameter<std::string>("pluginName");

    verbosity_ = cfg.exists("verbosity") ?
      cfg.getParameter<int>("verbosity") : 0;
  }
  virtual ~SVfitLikelihoodBase() {}

  const std::string& pluginType() const { return pluginType_; }
  const std::string& pluginName() const { return pluginName_; }

  virtual void beginJob(SVfitAlgorithmBase*) {}
  virtual void beginEvent(const edm::Event&, const edm::EventSetup&) {}
 
  virtual void print(std::ostream& stream) const
  {
    stream << "<SVfitLikelihoodBase::print>:" << std::endl;
    stream << " pluginName = " << pluginName_ << std::endl;
    stream << " pluginType = " << pluginType_ << std::endl;
  }

 protected:
  std::string pluginType_;
  std::string pluginName_;

  int verbosity_;
};

#endif
