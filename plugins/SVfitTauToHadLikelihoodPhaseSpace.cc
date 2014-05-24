#include "TauAnalysis/SVfit/plugins/SVfitTauToHadLikelihoodPhaseSpace.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "TauAnalysis/SVfit/interface/SVfitAlgorithmBase.h"
#include "TauAnalysis/SVfit/interface/svFitAuxFunctions.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitTauToHadHypothesis.h"

using namespace svFit_namespace;
using namespace svFitTauToHadLikelihoodPhaseSpace_namespace;

#define SVFIT_DEBUG 1

//namespace svFitTauToHadLikelihoodPhaseSpace_namespace
//{
//  std::map<std::string, TFile*> lutEntryType::inputFiles_;
//  std::map<std::string, int> lutEntryType::inputFileUsage_;
//}

namespace
{
  std::vector<lutEntryType*> readCfgLUT(const std::string& pluginName, const edm::VParameterSet& cfg) 
  {
    std::vector<lutEntryType*> lutEntries;
    for ( edm::VParameterSet::const_iterator cfg_i = cfg.begin();
	  cfg_i != cfg.end(); ++cfg_i ) {
      lutEntryType* lutEntry = new lutEntryType(pluginName, *cfg_i);
      lutEntries.push_back(lutEntry);
    }	
    return lutEntries;
  }
}

SVfitTauToHadLikelihoodPhaseSpace::SVfitTauToHadLikelihoodPhaseSpace(const edm::ParameterSet& cfg)
  : SVfitSingleParticleLikelihood(cfg),
    lutVisMass_(0),
    lutVisMass_shift_(0),
    lutVisPt_shift_(0),
    algorithm_(0)
{
  applySinThetaFactor_ = cfg.exists("applySinThetaFactor") ?
    cfg.getParameter<bool>("applySinThetaFactor") : false;

  varyVisMass_ = cfg.getParameter<bool>("varyVisMass");
  if ( varyVisMass_ ) {    
    lutVisMass_ = new lutEntryType(pluginName_, cfg);
  }
  
  shiftVisMass_ = cfg.getParameter<bool>("shiftVisMass");
  if ( shiftVisMass_ ) {
    edm::VParameterSet cfgVisMassResolution = cfg.getParameter<edm::VParameterSet>("visMassResolution");
    lutEntriesVisMass_shift_ = readCfgLUT(pluginName_, cfgVisMassResolution);
  }
  shiftVisPt_ = cfg.getParameter<bool>("shiftVisPt");
  if ( shiftVisPt_ ) {
    edm::VParameterSet cfgVisPtResolution = cfg.getParameter<edm::VParameterSet>("visPtResolution");
    lutEntriesVisPt_shift_ = readCfgLUT(pluginName_, cfgVisPtResolution);
  }
}

SVfitTauToHadLikelihoodPhaseSpace::~SVfitTauToHadLikelihoodPhaseSpace()
{
  delete lutVisMass_;
  for ( std::vector<lutEntryType*>::iterator it = lutEntriesVisMass_shift_.begin();
	it != lutEntriesVisMass_shift_.end(); ++it ) {
    delete (*it);
  }
  for ( std::vector<lutEntryType*>::iterator it = lutEntriesVisPt_shift_.begin();
	it != lutEntriesVisPt_shift_.end(); ++it ) {
    delete (*it);
  }
}

void SVfitTauToHadLikelihoodPhaseSpace::beginJob(SVfitAlgorithmBase* algorithm)
{
  algorithm_ = algorithm;

  algorithm->requestFitParameter(prodParticleLabel_,   svFit_namespace::kTau_visEnFracX, pluginName_);
  algorithm->requestFitParameter(prodParticleLabel_,   svFit_namespace::kTau_phi_lab,    pluginName_);
  if ( varyVisMass_ || shiftVisMass_ ) {
    algorithm->requestFitParameter(prodParticleLabel_, svFit_namespace::kTau_visMass,    pluginName_);
  }
  if ( shiftVisPt_ ) {
    algorithm->requestFitParameter(prodParticleLabel_, svFit_namespace::kTau_shiftVisPt, pluginName_);
  }
}

namespace
{
  int getTauDecayMode(const SVfitSingleParticleHypothesis* hypothesis)
  {
    const SVfitTauToHadHypothesis* hypothesis_T = dynamic_cast<const SVfitTauToHadHypothesis*>(hypothesis);
    assert(hypothesis_T != 0);
    
    const pat::Tau* tauJet = dynamic_cast<const pat::Tau*>(hypothesis_T->particle().get());
    assert(tauJet);
    
    int tauDecayMode = tauJet->decayMode();
    
    return tauDecayMode;
  }

  const lutEntryType* findLUT(const SVfitSingleParticleHypothesis* hypothesis, const std::vector<lutEntryType*>& lutEntries) 
  {
    //std::cout << "<findLUT>:" << std::endl;
    
    int tauDecayMode = getTauDecayMode(hypothesis);
    //std::cout << " tauDecayMode = " << tauDecayMode << std::endl;
    
    const lutEntryType* retVal = 0;
    for ( std::vector<lutEntryType*>::const_iterator lutEntry = lutEntries.begin();
	  lutEntry != lutEntries.end(); ++lutEntry ) {
      for ( std::vector<int>::const_iterator lutEntry_tauDecayMode = (*lutEntry)->tauDecayModes_.begin();
	    lutEntry_tauDecayMode != (*lutEntry)->tauDecayModes_.end(); ++lutEntry_tauDecayMode ) {
	//std::cout << "lutEntry->tauDecayMode = " << (*lutEntry_tauDecayMode) << std::endl;
	if ( (*lutEntry_tauDecayMode) == tauDecayMode || (*lutEntry_tauDecayMode) == -1 ) {
	  retVal = (*lutEntry);
	  break;
	}
      }
    }
    
    //std::cout << "--> retVal = " << retVal << std::endl;
    return retVal;
  }
}

void SVfitTauToHadLikelihoodPhaseSpace::beginCandidate(const SVfitSingleParticleHypothesis* hypothesis)
{
#ifdef SVFIT_DEBUG     
  if ( this->verbosity_ ) {
    std::cout << "<SVfitTauToHadLikelihoodPhaseSpace::beginCandidate (pluginName = " << pluginName_ << ")>:" << std::endl;
    std::cout << " shiftVisMass = " << shiftVisMass_ << std::endl;
    std::cout << " shiftVisPt = " << shiftVisPt_ << std::endl;
  }
#endif
  if ( shiftVisMass_ ) {
    lutVisMass_shift_ = findLUT(hypothesis, lutEntriesVisMass_shift_);
    if ( !lutVisMass_shift_ ) throw cms::Exception("SVfitTauToHadLikelihoodPhaseSpace::beginCandidate") 
      << " No tau mass resolution defined for tau decay mode = " << getTauDecayMode(hypothesis) << " !!\n";
    visMass_unshifted_ = hypothesis->p4().mass();
    if ( visMass_unshifted_ < chargedPionMass ) visMass_unshifted_ = chargedPionMass;
    if ( visMass_unshifted_ > tauLeptonMass   ) visMass_unshifted_ = tauLeptonMass;
  }
  if ( shiftVisPt_ ) {
    lutVisPt_shift_ = findLUT(hypothesis, lutEntriesVisPt_shift_);
    if ( !lutVisPt_shift_ ) throw cms::Exception("SVfitTauToHadLikelihoodPhaseSpace::beginCandidate") 
      << " No tau Pt resolution defined for tau decay mode = " << getTauDecayMode(hypothesis) << " !!\n";
  }
#ifdef SVFIT_DEBUG     
  if ( this->verbosity_ ) {
    std::cout << " lutVisMass_shift = " << lutVisMass_shift_ << std::endl;
    std::cout << " lutVisPt_shift = " << lutVisPt_shift_ << std::endl;
  }
#endif
}

double SVfitTauToHadLikelihoodPhaseSpace::operator()(const SVfitSingleParticleHypothesis* hypothesis, int polSign) const
{
//--- compute negative log-likelihood for tau lepton decay "leg"
//    to be compatible with three-body decay,
//    assuming constant matrix element,
//    so that energy and angular distribution of decay products is solely determined by phase-space
//
//    NOTE: the parametrization of the three-body decay phase-space is taken from the PDG:
//          K. Nakamura et al. (Particle Data Group), J. Phys. G 37, 075021 (2010);
//          formula 38.20a
//
  const SVfitTauToHadHypothesis* hypothesis_T = dynamic_cast<const SVfitTauToHadHypothesis*>(hypothesis);
  assert(hypothesis_T != 0);
  
  double decayAngle = hypothesis_T->gjAngle();  
  double visEnFracX = hypothesis_T->visEnFracX();
  double visMass = hypothesis_T->visMass();
  if ( !(varyVisMass_ || shiftVisMass_) ) {
    if ( visMass < chargedPionMass ) visMass = chargedPionMass;
    if ( visMass > tauLeptonMass   ) visMass = tauLeptonMass;
  }
  double visMass2 = square(visMass);
  double Pvis_rf = hypothesis_T->p4vis_rf().P();
#ifdef SVFIT_DEBUG     
  if ( this->verbosity_ ) {
    std::cout << "<SVfitTauToHadLikelihoodPhaseSpace::operator()>:" << std::endl;
    std::cout << " decayAngle = " << decayAngle << std::endl;
    std::cout << " visEnFracX = " << visEnFracX << std::endl;
    std::cout << " visMass = " << visMass << std::endl;
    std::cout << " angle(vis, tau) = " << TMath::ACos(compScalarProduct(normalize(hypothesis_T->p4().Vect()), normalize(hypothesis_T->p4_fitted().Vect()))) << std::endl;
  }
#endif
  // CV: normalize likelihood function such that 
  //               1
  //       integral  prob dX = 1.
  //               0
  double prob = tauLeptonMass/(2.*Pvis_rf);
  if ( visEnFracX < (visMass2/tauLeptonMass2) ) {
    double visEnFracX_limit = visMass2/tauLeptonMass2;
    prob /= (1. + 1.e+6*square(visEnFracX - visEnFracX_limit));
  } else if ( visEnFracX > 1. ) {
    double visEnFracX_limit = 1.;
    prob /= (1. + 1.e+6*square(visEnFracX - visEnFracX_limit));
  }
  if ( applySinThetaFactor_ ) prob *= (0.5*TMath::Sin(decayAngle));

  if ( varyVisMass_ ) {    
    assert(lutVisMass_);
    int bin = lutVisMass_->histogram_->FindBin(visMass);
    if ( bin < lutVisMass_->firstBin_ ) bin = lutVisMass_->firstBin_;
    if ( bin > lutVisMass_->lastBin_  ) bin = lutVisMass_->lastBin_;
    prob *= lutVisMass_->histogram_->GetBinContent(bin);
  }

  if ( shiftVisMass_ ) {    
    assert(lutVisMass_shift_);
    double deltaVisMass = visMass_unshifted_ - visMass;
    int bin = lutVisMass_shift_->histogram_->FindBin(deltaVisMass);
    if ( bin < lutVisMass_shift_->firstBin_ ) bin = lutVisMass_shift_->firstBin_;
    if ( bin > lutVisMass_shift_->lastBin_  ) bin = lutVisMass_shift_->lastBin_;
#ifdef SVFIT_DEBUG 
    if ( verbosity_ ) {
      //std::cout << "lutVisMass = " << lutVisMass_shift_->histogram_->GetName() << " (type = " << lutVisMass_shift_->histogram_->ClassName() << ")" << std::endl;
      std::cout << "deltaVisMass = " << deltaVisMass << ": probCorr = " << lutVisMass_shift_->histogram_->GetBinContent(bin) << std::endl;
      //std::cout << "bin = " << bin << " (firstBin = " << lutVisMass_shift_->firstBin_ << ", lastBin = " << lutVisMass_shift_->lastBin_ << ")" << std::endl;
      //std::cout << "probCorr = " << lutVisMass_shift_->histogram_->GetBinContent(bin) << std::endl;
    }
#endif
    prob *= lutVisMass_shift_->histogram_->GetBinContent(bin);
  }
  if ( shiftVisPt_ ) {    
    assert(lutVisPt_shift_);
    double recTauPtDivGenTauPt = ( hypothesis_T->p4Vis_shifted().pt() > 0. ) ?
      (hypothesis_T->p4().pt()/hypothesis_T->p4Vis_shifted().pt()) : 1.e+3;
    int bin = lutVisPt_shift_->histogram_->FindBin(recTauPtDivGenTauPt);
    if ( bin < lutVisPt_shift_->firstBin_ ) bin = lutVisPt_shift_->firstBin_;
    if ( bin > lutVisPt_shift_->lastBin_  ) bin = lutVisPt_shift_->lastBin_;
#ifdef SVFIT_DEBUG 
    if ( verbosity_ ) {
      //std::cout << "lutVisPt = " << lutVisPt_shift_->histogram_->GetName() << " (type = " << lutVisPt_shift_->histogram_->ClassName() << ")" << std::endl;
      std::cout << "recTauPtDivGenTauPt = " << recTauPtDivGenTauPt << ": probCorr = " << lutVisPt_shift_->histogram_->GetBinContent(bin) << std::endl;
      //std::cout << "bin = " << bin << " (firstBin = " << lutVisPt_shift_->firstBin_ << ", lastBin = " << lutVisPt_shift_->lastBin_ << ")" << std::endl;
      //std::cout << "probCorr = " << lutVisPt_shift_->histogram_->GetBinContent(bin) << std::endl;
    }
#endif
    prob *= lutVisPt_shift_->histogram_->GetBinContent(bin);
    // CV: account for Jacobi factor 
    double genTauPtDivRecTauPt = ( recTauPtDivGenTauPt > 0. ) ? 
      (1./recTauPtDivGenTauPt) : 1.e+1;
    prob *= genTauPtDivRecTauPt;
  }
  
  if ( applyVisPtCutCorrection_ ) {
    double probCorr = 1.;
    const double epsilon_regularization = 1.e-1;
    if ( hypothesis_T->p4_fitted().pt() > visPtCutThreshold_ ) {     
      double xCut = visPtCutThreshold_/hypothesis_T->p4_fitted().pt();      
      probCorr = 1./((1. - xCut) + epsilon_regularization);
    }
#ifdef SVFIT_DEBUG 
    if ( verbosity_ ) std::cout << "probCorr (had) = " << probCorr << std::endl;
#endif
    prob *= probCorr;
  }

  if ( algorithm_->applyJacobiFactors() && visEnFracX > 0. ) {
    double jacobiFactor = 1./(cube(visEnFracX));
#ifdef SVFIT_DEBUG 
    if ( verbosity_ ) std::cout << "jacobiFactor = " << jacobiFactor << std::endl;
#endif
    prob *= jacobiFactor;
  }

#ifdef SVFIT_DEBUG       
  if ( verbosity_ ) std::cout << "--> prob = " << prob << std::endl;
#endif
  return prob;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(SVfitSingleParticleLikelihoodPluginFactory, SVfitTauToHadLikelihoodPhaseSpace, "SVfitTauToHadLikelihoodPhaseSpace");
