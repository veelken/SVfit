#ifndef TauAnalysis_SVfit_SVfitProducer_h
#define TauAnalysis_SVfit_SVfitProducer_h

/** \class SVfitProducer
 *
 * Produce data-formats storing solutions of nSVfit algorithm
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.3 $
 *
 * $Id: SVfitProducerT.h,v 1.3 2012/04/09 16:48:48 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "TauAnalysis/SVfit/interface/SVfitAlgorithmBase.h"

#include <TStopwatch.h>

#include <vector>
#include <string>

template <typename T>
class SVfitProducerT : public edm::EDProducer 
{
 public:
  explicit SVfitProducerT(const edm::ParameterSet&);
  ~SVfitProducerT();

  void beginJob();

  void produce(edm::Event&, const edm::EventSetup&);

  void endJob();

 private:
  std::string moduleLabel_;

  std::string instanceLabel_;

  SVfitAlgorithmBase* algorithm_;
  
  typedef std::vector<T> SVfitEventHypothesisCollection;

  typedef std::vector<std::string> vstring;
  vstring inputParticleNames_;
  typedef std::vector<edm::InputTag> vInputTag;
  vInputTag srcInputParticles_;
  unsigned numInputParticles_;

  double dRmin_; // minimum eta-phi separation between any pair of input particles

  edm::InputTag srcMEt_;
  edm::InputTag srcPrimaryVertex_;

  TStopwatch* timer_;
  long numSVfitCalls_;
  unsigned instanceId_;
  static unsigned instanceCounter_;
};

#endif

