#ifndef TauAnalysis_SVfit_SVfitEventLikelihoodMEt3_h
#define TauAnalysis_SVfit_SVfitEventLikelihoodMEt3_h

/** \class SVfitEventLikelihoodMEt3
 *
 * Plugin for computing likelihood for neutrinos produced in tau lepton decays
 * to match missing transverse momentum reconstructed in the event
 *
 * New version using covariance matrix of (PF)MET significance calculation
 * (CMS AN-10/400) to compute the likehood
 *
 * \author Christian Veelken, UC Davis
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "RecoMET/METAlgorithms/interface/SignAlgoResolutions.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "TauAnalysis/SVfit/interface/SVfitEventLikelihood.h"
#include "TauAnalysis/SVfit/interface/PFMEtSignInterface.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitEventHypothesis.h"

#include <TH2.h>
#include <TAxis.h>
#include <TRandom3.h>

#include <string>
#include <list>

class SVfitEventLikelihoodMEt3 : public SVfitEventLikelihood
{
 public:
  SVfitEventLikelihoodMEt3(const edm::ParameterSet&);
  ~SVfitEventLikelihoodMEt3();

  void beginJob(SVfitAlgorithmBase*);
  void beginEvent(const edm::Event&, const edm::EventSetup&);
  void beginCandidate(const SVfitEventHypothesis*) const;

  double operator()(const SVfitEventHypothesis*) const;

 private:

  double power_;

  edm::InputTag srcPFJets_;
  edm::InputTag srcPFCandidates_;

  double pfCandPtThreshold_;
  double pfJetPtThreshold_;
  TH2* lut_;
  mutable TAxis* xAxis_;
  mutable int numBinsX_;
  mutable TAxis* yAxis_;
  mutable int numBinsY_;

  unsigned numToys_;

  std::list<const reco::PFJet*> pfJetListForCovMatrix_;
  std::list<const reco::PFCandidate*> pfCandidateListForCovMatrix_;
  std::list<const reco::PFJet*> pfJetListForToys_;
  std::list<const reco::PFCandidate*> pfCandidateListForToys_;
  
  double dRoverlapPFJet_;
  double dRoverlapPFCandidate_;

  PFMEtSignInterfaceBase* pfMEtSign_;

  mutable TRandom3 rnd_;

  bool monitorMEtUncertainty_;
  std::string monitorFilePath_;
  std::string monitorFileName_;
};

#endif
