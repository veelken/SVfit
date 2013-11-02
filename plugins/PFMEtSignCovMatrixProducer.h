#ifndef TauAnalysis_SVfit_PFMEtSignCovMatrixProducer_h  
#define TauAnalysis_SVfit_PFMEtSignCovMatrixProducer_h

/** \class PFMEtSignCovMatrixProducer
 *
 * Produce covariance matrix representing expected PFMET resolution,
 * computed on event-by-event bases by the (PF)MET significance algorithm
 * (see CMS AN-10/400 for description of the (PF)MEt significance computation)
 *
 * \author Christian Veelken, LLR
 *
 * \version $Revision: 1.3 $
 *
 * $Id: PFMEtSignCovMatrixProducer.h,v 1.3 2012/09/01 08:45:45 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "JetMETCorrections/METPUSubtraction/interface/PFMEtSignInterfaceBase.h"

class PFMEtSignCovMatrixProducer : public edm::EDProducer 
{
 public:
  
  PFMEtSignCovMatrixProducer(const edm::ParameterSet&);
  ~PFMEtSignCovMatrixProducer();

  void produce(edm::Event&, const edm::EventSetup&);

 private:

  std::string moduleLabel_;

  typedef std::vector<edm::InputTag> vInputTag;
  vInputTag src_;

  PFMEtSignInterfaceBase pfMEtSignInterface_;
};

#endif


