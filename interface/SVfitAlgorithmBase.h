#ifndef TauAnalysis_SVfit_SVfitAlgorithmBase_h
#define TauAnalysis_SVfit_SVfitAlgorithmBase_h

/** \class SVfitAlgorithmBase
 *
 * Abstract base-class for plugins finding best (n)SVfit solution,
 * either by integration or by fitting
 *
 * \author Christian Veelken, UC Davis
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/Ptr.h"

#include "TauAnalysis/SVfit/interface/SVfitEventLikelihood.h"
#include "TauAnalysis/SVfit/interface/SVfitResonanceLikelihood.h"
#include "TauAnalysis/SVfit/interface/SVfitSingleParticleLikelihood.h"
#include "TauAnalysis/SVfit/interface/SVfitEventBuilderBase.h"
#include "TauAnalysis/SVfit/interface/SVfitParameter.h"
#include "TauAnalysis/SVfit/interface/SVfitTrackService.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitEventHypothesis.h"
#include "AnalysisDataFormats/SVfit/interface/SVfitResonanceHypothesis.h"
#include "AnalysisDataFormats/SVfit/interface/SVfitSingleParticleHypothesis.h"

#include <vector>
#include <string>

class SVfitAlgorithmBase
{
 public:
  SVfitAlgorithmBase(const edm::ParameterSet&);
  virtual ~SVfitAlgorithmBase();

  virtual void beginJob();
  virtual void beginEvent(const edm::Event&, const edm::EventSetup&);

  virtual void requestFitParameter(const std::string&, int, const std::string&);
  virtual void fixFitParameter(int);

  virtual unsigned getNumFitParameter(const std::string&) const;
  virtual SVfitParameter* getFitParameter(const std::string&, int) const;
  virtual SVfitParameter* getFitParameter(int) const;

  virtual void setFitParameterInitialValue(int, double);
  virtual void setFitParameterLimit(int, double, double);
  virtual void setFitParameterStepSize(int, double);

  virtual bool applyJacobiFactors() const = 0;

  virtual void print(std::ostream&) const {}

  // NOTE: fit method creates a new object of type SVfitEventHypothesisBase (or derrived class);
  //       ownership of this object is held by calling code
  //      --> calling code needs to take-care of deleting this object, in order to avoid memory leak
  typedef edm::Ptr<reco::Candidate> CandidatePtr;
  typedef std::map<std::string, CandidatePtr> inputParticleMap;
  virtual SVfitEventHypothesisBase* fit(const inputParticleMap&, const reco::Vertex*) const;

  virtual bool update(const double* x, const double* param) const;
  virtual double nll(const double* x, const double* param) const;

  static const SVfitAlgorithmBase* gSVfitAlgorithm;

  const SVfitEventHypothesis* currentEventHypothesis() const { return currentEventHypothesis_; }

  friend class SVfitTauLikelihoodTrackInfo;

 protected:
  virtual void fitImp() const = 0;

  void setMassResults(SVfitResonanceHypothesisBase*, double, double, double) const;

  std::string pluginName_;
  std::string pluginType_;

  struct daughterModelType
  {
    daughterModelType(const std::string& daughterName, const edm::ParameterSet& cfg,
		      std::vector<SVfitLikelihoodBase*>& allLikelihoods)
      : daughterName_(daughterName),
	prodParticleLabel_(cfg.getParameter<std::string>("prodParticleLabel"))
    {
      if ( cfg.exists("likelihoodFunctions") ) {
	typedef std::vector<edm::ParameterSet> vParameterSet;
	vParameterSet cfg_likelihoods = cfg.getParameter<vParameterSet>("likelihoodFunctions");
	for ( vParameterSet::iterator cfg_likelihood = cfg_likelihoods.begin();
	      cfg_likelihood != cfg_likelihoods.end(); ++cfg_likelihood ) {
	  cfg_likelihood->addParameter<std::string>("prodParticleLabel", prodParticleLabel_);
	  std::string pluginType = cfg_likelihood->getParameter<std::string>("pluginType");
	  SVfitSingleParticleLikelihood* likelihood =
  	    SVfitSingleParticleLikelihoodPluginFactory::get()->create(pluginType, *cfg_likelihood);
	  likelihoods_.push_back(likelihood);
	  allLikelihoods.push_back(likelihood);
	}
      }
    }
    ~daughterModelType()
    {
      for ( std::vector<SVfitSingleParticleLikelihood*>::iterator it = likelihoods_.begin();
	    it != likelihoods_.end(); ++it ) {
	delete (*it);
      }
    }
    void beginCandidate(const SVfitSingleParticleHypothesis* hypothesis)
    {
      for ( std::vector<SVfitSingleParticleLikelihood*>::const_iterator likelihood = likelihoods_.begin();
	    likelihood != likelihoods_.end(); ++likelihood ) {
	(*likelihood)->beginCandidate(hypothesis);
      }
    }
    double prob(const SVfitSingleParticleHypothesis* hypothesis, int idxPolState) const
    {
      double retVal = 1.;
      for ( std::vector<SVfitSingleParticleLikelihood*>::const_iterator likelihood = likelihoods_.begin();
	    likelihood != likelihoods_.end(); ++likelihood ) {
	retVal *= (**likelihood)(hypothesis, hypothesis->polSign(idxPolState));
      }
      return retVal;
    }
    std::string daughterName_;
    std::string prodParticleLabel_;
    std::vector<SVfitSingleParticleLikelihood*> likelihoods_;
  };

  struct resonanceModelType
  {
    resonanceModelType(const std::string& resonanceName, const edm::ParameterSet& cfg,
		       std::vector<SVfitLikelihoodBase*>& allLikelihoods)
      : resonanceName_(resonanceName)
    {
      if ( cfg.exists("likelihoodFunctions") ) {
	typedef std::vector<edm::ParameterSet> vParameterSet;
	vParameterSet cfg_likelihoods = cfg.getParameter<vParameterSet>("likelihoodFunctions");
        for ( vParameterSet::const_iterator cfg_likelihood = cfg_likelihoods.begin();
	      cfg_likelihood != cfg_likelihoods.end(); ++cfg_likelihood ) {
  	  std::string pluginType = cfg_likelihood->getParameter<std::string>("pluginType");
	  SVfitResonanceLikelihood* likelihood =
	    SVfitResonanceLikelihoodPluginFactory::get()->create(pluginType, *cfg_likelihood);
	    likelihoods_.push_back(likelihood);
	    allLikelihoods.push_back(likelihood);
        }
      }

      edm::ParameterSet cfg_daughters = cfg.getParameter<edm::ParameterSet>("daughters");
      typedef std::vector<std::string> vstring;
      vstring daughterNames = cfg_daughters.getParameterNamesForType<edm::ParameterSet>();
      for ( vstring::const_iterator daughterName = daughterNames.begin();
	    daughterName != daughterNames.end(); ++daughterName ) {
        edm::ParameterSet cfg_daughter = cfg_daughters.getParameter<edm::ParameterSet>(*daughterName);
	cfg_daughter.addParameter<std::string>("prodParticleLabel", *daughterName);
	daughters_.push_back(new daughterModelType(*daughterName, cfg_daughter, allLikelihoods));
      }
      numDaughters_ = daughters_.size();
    }
    ~resonanceModelType()
    {
      for ( std::vector<SVfitResonanceLikelihood*>::iterator it = likelihoods_.begin();
	    it != likelihoods_.end(); ++it ) {
	delete (*it);
      }
      for ( std::vector<daughterModelType*>::iterator it = daughters_.begin();
	    it != daughters_.end(); ++it ) {
	delete (*it);
      }
    }
    void beginCandidate(const SVfitResonanceHypothesis* hypothesis)
    {
      for ( std::vector<SVfitResonanceLikelihood*>::const_iterator likelihood = likelihoods_.begin();
	    likelihood != likelihoods_.end(); ++likelihood ) {
	(*likelihood)->beginCandidate(hypothesis);
      }
      assert(hypothesis->numDaughters() == numDaughters_);
      for ( unsigned iDaughter = 0; iDaughter < numDaughters_; ++iDaughter ) {
	daughters_[iDaughter]->beginCandidate(hypothesis->daughter(iDaughter));
      }
    }
    double prob(const SVfitResonanceHypothesis* hypothesis, int idxPolState) const
    {
      double retVal = 1.;
      for ( std::vector<SVfitResonanceLikelihood*>::const_iterator likelihood = likelihoods_.begin();
	    likelihood != likelihoods_.end(); ++likelihood ) {
	retVal *= (**likelihood)(hypothesis, hypothesis->polHandedness(idxPolState));
      }
      assert(hypothesis->numDaughters() == numDaughters_);
      for ( unsigned iDaughter = 0; iDaughter < numDaughters_; ++iDaughter ) {
	retVal *= daughters_[iDaughter]->prob(hypothesis->daughter(iDaughter), idxPolState);
      }
      return retVal;
    }
    std::string resonanceName_;
    std::vector<SVfitResonanceLikelihood*> likelihoods_;
    std::vector<daughterModelType*> daughters_;
    unsigned numDaughters_;
  };

  struct eventModelType
  {
    eventModelType(const edm::ParameterSet& cfg, std::vector<SVfitLikelihoodBase*>& allLikelihoods)
    {
      edm::ParameterSet cfg_builder = cfg.getParameter<edm::ParameterSet>("builder");
      cfg_builder.addParameter<edm::ParameterSet>("resonances", cfg.getParameter<edm::ParameterSet>("resonances"));
      std::string pluginType_builder = cfg_builder.getParameter<std::string>("pluginType");
      builder_ = SVfitEventBuilderPluginFactory::get()->create(pluginType_builder, cfg_builder);

      if ( cfg.exists("likelihoodFunctions") ) {
	typedef std::vector<edm::ParameterSet> vParameterSet;
	vParameterSet cfg_likelihoods = cfg.getParameter<vParameterSet>("likelihoodFunctions");
	for ( vParameterSet::const_iterator cfg_likelihood = cfg_likelihoods.begin();
	      cfg_likelihood != cfg_likelihoods.end(); ++cfg_likelihood ) {
	  std::string pluginType = cfg_likelihood->getParameter<std::string>("pluginType");
	  SVfitEventLikelihood* likelihood =
  	    SVfitEventLikelihoodPluginFactory::get()->create(pluginType, *cfg_likelihood);
	  likelihoods_.push_back(likelihood);
	  allLikelihoods.push_back(likelihood);
	}
      }

      edm::ParameterSet cfg_resonances = cfg.getParameter<edm::ParameterSet>("resonances");
      typedef std::vector<std::string> vstring;
      vstring resonanceNames = cfg_resonances.getParameterNamesForType<edm::ParameterSet>();
      for ( vstring::const_iterator resonanceName = resonanceNames.begin();
	    resonanceName != resonanceNames.end(); ++resonanceName ) {
        edm::ParameterSet cfg_resonance = cfg_resonances.getParameter<edm::ParameterSet>(*resonanceName);
	resonances_.push_back(new resonanceModelType(*resonanceName, cfg_resonance, allLikelihoods));
      }
      numResonances_ = resonances_.size();
    }
    ~eventModelType()
    {
      delete builder_;
      for ( std::vector<SVfitEventLikelihood*>::iterator it = likelihoods_.begin();
	    it != likelihoods_.end(); ++it ) {
	delete (*it);
      }
      for ( std::vector<resonanceModelType*>::iterator it = resonances_.begin();
	    it != resonances_.end(); ++it ) {
	delete (*it);
      }
    }
    void beginCandidate(const SVfitEventHypothesis* hypothesis)
    {
      for ( std::vector<SVfitEventLikelihood*>::const_iterator likelihood = likelihoods_.begin();
	    likelihood != likelihoods_.end(); ++likelihood ) {
	(*likelihood)->beginCandidate(hypothesis);
      }
      assert(hypothesis->numResonances() == numResonances_);
      for ( unsigned iResonance = 0; iResonance < numResonances_; ++iResonance ) {
	resonances_[iResonance]->beginCandidate(hypothesis->resonance(iResonance));
      }
    }
    double nll(const SVfitEventHypothesis* hypothesis) const
    {
      double prob = 1.;
      for ( std::vector<SVfitEventLikelihood*>::const_iterator likelihood = likelihoods_.begin();
	    likelihood != likelihoods_.end(); ++likelihood ) {
	prob *= (**likelihood)(hypothesis);
      }
      //std::cout << "prob (1) = " << prob << std::endl;
      assert(hypothesis->numResonances() == numResonances_);
      unsigned numPolStates_event = hypothesis->numPolStates();
      if ( numPolStates_event == 1 ) { // Htautau and WH case
	for ( unsigned iResonance = 0; iResonance < numResonances_; ++iResonance ) {
	  const SVfitResonanceHypothesis* resonance_hypothesis = hypothesis->resonance(iResonance);
	  unsigned numPolStates_resonance = resonance_hypothesis->numPolStates();
	  if ( numPolStates_resonance == 1 ) {          
	    prob *= resonances_[iResonance]->prob(resonance_hypothesis, 0);
	  } else {
	    double prob_polSum = 0.;
	    for ( unsigned idxPolState = 0; idxPolState < numPolStates_resonance; ++idxPolState ) {
	      //std::cout << "idxPolState = " << idxPolState << ":" << std::endl;
	      double prob_pol = resonances_[iResonance]->prob(resonance_hypothesis, idxPolState);
	      //std::cout << " prob_pol = " << prob_pol << std::endl;
	      prob_polSum += prob_pol;
	    }
	    //std::cout << "prob_polSum = " << prob_polSum << std::endl;
	    prob *= prob_polSum;
	  }
	}
      } else { // WW case
	double prob_polSum = 0.;
	for ( unsigned idxPolState = 0; idxPolState < numPolStates_event; ++idxPolState ) {
	  //std::cout << "idxPolState = " << idxPolState << ":" << std::endl;
	  double prob_pol = 1.;
	  for ( unsigned iResonance = 0; iResonance < numResonances_; ++iResonance ) {
	    const SVfitResonanceHypothesis* resonance_hypothesis = hypothesis->resonance(iResonance);
	    unsigned numPolStates_resonance = resonance_hypothesis->numPolStates();
	    assert(numPolStates_resonance == numPolStates_event);
	    prob_pol *= resonances_[iResonance]->prob(resonance_hypothesis, idxPolState);
	  }
	  //std::cout << " prob_pol = " << prob_pol << std::endl;
	  prob_polSum += prob_pol;
	}
	//std::cout << "prob_polSum = " << prob_polSum << std::endl;
	prob *= prob_polSum;
      }
      //std::cout << "--> prob (2) = " << prob << std::endl;
      double retVal;
      if ( prob > 0. ) retVal = -TMath::Log(prob);
      else retVal = std::numeric_limits<float>::max();
      return retVal;
    }
    SVfitEventBuilderBase* builder_;
    std::vector<SVfitEventLikelihood*> likelihoods_;
    std::vector<resonanceModelType*> resonances_;
    unsigned numResonances_;
  };

  eventModelType* eventModel_;

  std::vector<SVfitLikelihoodBase*> allLikelihoods_; // list of all event, resonance and single particle likelihood plugins
                                                     // NOTE: plugins in list are **not** owned by SVfitAlgorithmBase

  edm::Service<SVfitTrackService> trackService_;
  const edm::EventSetup* currentEventSetup_;
  mutable SVfitEventHypothesis* currentEventHypothesis_;
  mutable bool currentEventHypothesis_isValidSolution_;
  mutable SVfitEventHypothesisBase* fittedEventHypothesis_;
  mutable double fittedEventHypothesis_nll_;

  mutable std::vector<SVfitParameter> fitParameters_;
  int fitParameterCounter_;

  int verbosity_;
};

#include "FWCore/PluginManager/interface/PluginFactory.h"

typedef edmplugin::PluginFactory<SVfitAlgorithmBase* (const edm::ParameterSet&)> SVfitAlgorithmPluginFactory;

#endif

