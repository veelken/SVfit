#include "TauAnalysis/SVfit/plugins/SVfitLeptonBuilderDummy.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Math/interface/deltaR.h"

#include <TMath.h>

template<typename T>
SVfitLeptonBuilderDummy<T>::SVfitLeptonBuilderDummy(const edm::ParameterSet& cfg)
  : SVfitSingleParticleBuilderBase(cfg)
{
  fixToGenP4_ = ( cfg.exists("fixToGenP4") ) ? cfg.getParameter<bool>("fixToGenP4") : false;
  if ( fixToGenP4_ ) {
    srcGenLeptons_ = cfg.getParameter<edm::InputTag>("srcGenLeptons");
    dRmatch_ = cfg.getParameter<double>("dRmatch");
  }
}

template<typename T>
void SVfitLeptonBuilderDummy<T>::beginEvent(const edm::Event& evt, const edm::EventSetup& es)
{
  SVfitBuilderBase::beginEvent(evt, es);
  
  if ( fixToGenP4_ ) {
    evt.getByLabel(srcGenLeptons_, genParticles_);
  }
}

template<typename T>
void SVfitLeptonBuilderDummy<T>::initialize(SVfitSingleLeptonHypothesis* hypothesis, const reco::Candidate* visCandidate) const
{
  if ( this->verbosity_ ) std::cout << "<SVfitLeptonBuilderDummy::initialize>:" << std::endl;

  hypothesis->p4_ = visCandidate->p4();

  if ( fixToGenP4_ ) {
    bool isMatched = false;
    int idx = 0;
    for ( GenParticleView::const_iterator genParticle = genParticles_->begin();
	  genParticle != genParticles_->end() && !isMatched; ++genParticle ) {
      if ( verbosity_ ) {
	std::cout << "genParticle #" << idx << " (pdgId = " << genParticle->pdgId() << "):" 
		  << " Pt = " << genParticle->pt() << ", eta = " << genParticle->eta() << ", phi = " << genParticle->phi() 
		  << " (status = " << genParticle->status() << ")" << std::endl;
      }

      if ( TMath::Abs(genParticle->pdgId()) == getLeptonPdgId() ) {
	const reco::GenParticle* genLepton = &(*genParticle);
	reco::Candidate::LorentzVector genLeptonP4 = genLepton->p4();
	double dR = deltaR(genLepton->p4(), visCandidate->p4());
	if ( verbosity_ ) std::cout << "dR = " << dR << std::endl;

	if ( dR < dRmatch_ ) {
	  genP4_ = genLeptonP4;

	  if ( verbosity_ ) {
	    std::cout << " found match." << std::endl;
	    std::cout << "initializing:" << std::endl;
	    std::cout << " genP4: Pt = " << genP4_.pt() << ", eta = " << genP4_.eta() << ", phi = " << genP4_.phi() << std::endl;
	  }

	  isMatched = true;	    
	}
      }
      ++idx;
    }
    
    if ( !isMatched )
      throw cms::Exception("SVfitLeptonBuilderDummy::initialize")
	<< "Failed to find gen. Lepton matching reconstructed lepton:" 
	<< " Pt = " << visCandidate->pt() << ", eta = " << visCandidate->eta() << ", phi = " << visCandidate->phi() << " !!\n";
    
    hypothesis->p4_ = genP4_;
  }

  hypothesis->p4_fitted_ = hypothesis->p4_;
  hypothesis->dp4_.SetPxPyPzE(0.,0.,0.,0.);
}

template<typename T>
SVfitSingleParticleHypothesis* SVfitLeptonBuilderDummy<T>::build(const SVfitSingleParticleBuilderBase::inputParticleMap& inputParticles) const 
{
  inputParticleMap::const_iterator particlePtr = inputParticles.find(prodParticleLabel_);
  assert(particlePtr != inputParticles.end());

  SVfitSingleLeptonHypothesis* hypothesis = 
    new SVfitSingleLeptonHypothesis(particlePtr->second, prodParticleLabel_, barcodeCounter_);
  if ( verbosity_ ) {
    std::cout << "<SVfitLeptonBuilderDummy::build>:" << std::endl;
    std::cout << " hypothesis #" << barcodeCounter_ << ": " << hypothesis << std::endl;
  }

  ++barcodeCounter_;
  
  initialize(hypothesis, particlePtr->second.get());
  
  return hypothesis;
}

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

template<>
int SVfitLeptonBuilderDummy<pat::Electron>::getLeptonPdgId() const 
{
  return 11;
}

template <>
int SVfitLeptonBuilderDummy<pat::Muon>::getLeptonPdgId() const
{
  return 13;
}

typedef SVfitLeptonBuilderDummy<pat::Electron> SVfitElectronBuilderDummy;
typedef SVfitLeptonBuilderDummy<pat::Muon> SVfitMuonBuilderDummy;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(SVfitSingleParticleBuilderPluginFactory, SVfitElectronBuilderDummy, "SVfitElectronBuilderDummy");
DEFINE_EDM_PLUGIN(SVfitSingleParticleBuilderPluginFactory, SVfitMuonBuilderDummy, "SVfitMuonBuilderDummy");

