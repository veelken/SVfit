#ifndef TauAnalysis_SVfit_PFMEtSignCovMatrixUnembedder_h
#define TauAnalysis_SVfit_PFMEtSignCovMatrixUnembedder_h

/** \class PFMEtSignCovMatrixUnembedder
 *
 * Extract covariance matrix representing expected PFMET resolution from reco::MET object
 * and add corresponding object of type PFMEtSignCovMatrix to the event
 *
 * \author Christian Veelken, LLR
 *
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class PFMEtSignCovMatrixUnembedder : public edm::EDProducer 
{
 public:
  
  PFMEtSignCovMatrixUnembedder(const edm::ParameterSet&);
  ~PFMEtSignCovMatrixUnembedder();

  void produce(edm::Event&, const edm::EventSetup&);

 private:

  std::string moduleLabel_;

  edm::InputTag src_;
};

#endif


