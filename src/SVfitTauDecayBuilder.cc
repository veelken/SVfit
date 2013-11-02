#include "TauAnalysis/SVfit/interface/SVfitTauDecayBuilder.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/TauReco/interface/PFTauDecayMode.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/CLHEP/interface/AlgebraicObjects.h"

#include "DataFormats/Math/interface/angle.h"

#include "TauAnalysis/SVfit/interface/SVfitAlgorithmBase.h"
#include "TauAnalysis/SVfit/interface/SVfitParameter.h"
#include "TauAnalysis/SVfit/interface/svFitAuxFunctions.h"
#include "TauAnalysis/SVfit/interface/candidateAuxFunctions.h"

#include <TMath.h>
#include <Math/VectorUtil.h>

using namespace svFit_namespace;

const unsigned minNumTracksFit = 3;

#define SVFIT_DEBUG 1

SVfitTauDecayBuilder::SVfitTauDecayBuilder(const edm::ParameterSet& cfg)
  : SVfitSingleParticleBuilderBase(cfg),
    algorithm_(0),
    trackBuilder_(0),
    decayVertexFitAlgorithm_(0),
    idxFitParameter_nuInvMass_(-1)
{
  edm::ParameterSet cfgTrackQualityCuts = cfg.getParameter<edm::ParameterSet>("trackQualityCuts");
  trackMinNumHits_      = cfgTrackQualityCuts.getParameter<unsigned>("minNumHits");
  trackMinNumPixelHits_ = cfgTrackQualityCuts.getParameter<unsigned>("minNumPixelHits");
  trackMaxChi2DoF_      = cfgTrackQualityCuts.getParameter<double>("maxChi2DoF");
  trackMaxDeltaPoverP_  = cfgTrackQualityCuts.getParameter<double>("maxDeltaPoverP");
  trackMinPt_           = cfgTrackQualityCuts.getParameter<double>("minPt");

  edm::ParameterSet cfgVertexFitAlgorithm(cfg);
  cfgVertexFitAlgorithm.addParameter<unsigned>("minNumTracksFit", minNumTracksFit);
  decayVertexFitAlgorithm_ = new SVfitDecayVertexFitter(cfgVertexFitAlgorithm);
  fitDecayVertex_ = cfg.getParameter<bool>("fitDecayVertex");

  fixToGenVisEnFracX_    = ( cfg.exists("fixToGenVisEnFracX")    ) ? cfg.getParameter<bool>("fixToGenVisEnFracX")    : false;
  fixToGenPhiLab_        = ( cfg.exists("fixToGenPhiLab")        ) ? cfg.getParameter<bool>("fixToGenPhiLab")        : false;
  fixToGenVisMass_       = ( cfg.exists("fixToGenVisMass")       ) ? cfg.getParameter<bool>("fixToGenVisMass")       : false;
  fixToGenNuInvMass_     = ( cfg.exists("fixToGenNuInvMass")     ) ? cfg.getParameter<bool>("fixToGenNuInvMass")     : false;
  fixRecToGenVertex_     = ( cfg.exists("fixRecToGenVertex")     ) ? cfg.getParameter<bool>("fixRecToGenVertex")     : false;
  fixToGenDecayDistance_ = ( cfg.exists("fixToGenDecayDistance") ) ? cfg.getParameter<bool>("fixToGenDecayDistance") : false;
  fixToGenVisP4_         = ( cfg.exists("fixToGenVisP4")         ) ? cfg.getParameter<bool>("fixToGenVisP4")         : false;

  initializeToGen_ = ( cfg.exists("initializeToGen") ) ? cfg.getParameter<bool>("initializeToGen") : false;
  if ( fixToGenVisEnFracX_    || 
       fixToGenPhiLab_        ||
       fixToGenVisMass_       ||
       fixToGenNuInvMass_     ||
       fixRecToGenVertex_     ||
       fixToGenDecayDistance_ ||
       fixToGenVisP4_         ) {
    initializeToGen_ = true;
  }

  if ( initializeToGen_ ) {
    srcGenTaus_ = cfg.getParameter<edm::InputTag>("srcGenTaus");
    dRmatch_ = cfg.getParameter<double>("dRmatch");
  }
}

SVfitTauDecayBuilder::~SVfitTauDecayBuilder()
{
  delete decayVertexFitAlgorithm_;
}

void SVfitTauDecayBuilder::beginJob(SVfitAlgorithmBase* algorithm)
{
  algorithm_ = algorithm;

  // Map the fit parameters to indices.
  idxFitParameter_visEnFracX_          = getFitParameterIdx(algorithm, prodParticleLabel_, nSVfit_namespace::kTau_visEnFracX);
  idxFitParameter_phi_lab_             = getFitParameterIdx(algorithm, prodParticleLabel_, nSVfit_namespace::kTau_phi_lab);
  idxFitParameter_visMass_             = getFitParameterIdx(algorithm, prodParticleLabel_, nSVfit_namespace::kTau_visMass, true);
  idxFitParameter_nuInvMass_           = getFitParameterIdx(algorithm, prodParticleLabel_, nSVfit_namespace::kTau_nuInvMass, true);
  idxFitParameter_decayDistance_shift_ = getFitParameterIdx(algorithm, prodParticleLabel_, nSVfit_namespace::kTau_decayDistance_lab_shift, true);
#ifdef SVFIT_DEBUG 
  if ( verbosity_ ) {
    std::cout << "<SVfitTauDecayBuilder::beginJob>:" << std::endl;
    std::cout << " pluginName = " << pluginName_ << std::endl;
    std::cout << " idxFitParameter_visEnFracX = " << idxFitParameter_visEnFracX_ << std::endl;
    std::cout << " idxFitParameter_phi_lab = " << idxFitParameter_phi_lab_ << std::endl;
    std::cout << " idxFitParameter_visMass = " << idxFitParameter_visMass_ << std::endl;
    std::cout << " idxFitParameter_nuInvMass = " << idxFitParameter_nuInvMass_ << std::endl;
    std::cout << " idxFitParameter_decayDistance_shift = " << idxFitParameter_decayDistance_shift_ << std::endl;
  }
#endif
  if ( fixToGenVisEnFracX_    && idxFitParameter_visEnFracX_          != -1 ) algorithm->fixFitParameter(idxFitParameter_visEnFracX_);
  if ( fixToGenPhiLab_        && idxFitParameter_phi_lab_             != -1 ) algorithm->fixFitParameter(idxFitParameter_phi_lab_);
  if ( fixToGenVisMass_       && idxFitParameter_visMass_             != -1 ) algorithm->fixFitParameter(idxFitParameter_visMass_);
  if ( fixToGenNuInvMass_     && idxFitParameter_nuInvMass_           != -1 ) algorithm->fixFitParameter(idxFitParameter_nuInvMass_);
  if ( fixToGenDecayDistance_ && idxFitParameter_decayDistance_shift_ != -1 ) algorithm->fixFitParameter(idxFitParameter_decayDistance_shift_);
}

void SVfitTauDecayBuilder::beginEvent(const edm::Event& evt, const edm::EventSetup& es)
{
  SVfitBuilderBase::beginEvent(evt, es);

//--- get pointer to TransientTrackBuilder
  edm::ESHandle<TransientTrackBuilder> trackBuilderHandle;
  es.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilderHandle);
  trackBuilder_ = trackBuilderHandle.product();
  if ( !trackBuilder_ ) {
    throw cms::Exception("SVfitTauDecayBuilder::beginEvent")
      << " Failed to access TransientTrackBuilder !!\n";
  }

  decayVertexFitAlgorithm_->beginEvent(evt, es);
  
  if ( initializeToGen_ ) {
    evt.getByLabel(srcGenTaus_, genParticles_);
  }
}

namespace
{
  bool isHigherPt(const reco::Track* track1, const reco::Track* track2)
  {
    return (track1->pt() > track2->pt());
  }
}

void SVfitTauDecayBuilder::initialize(SVfitTauDecayHypothesis* hypothesis, const reco::Candidate* visCandidate) const
{
#ifdef SVFIT_DEBUG 
  if ( this->verbosity_ ) std::cout << "<SVfitTauDecayBuilder::initialize>:" << std::endl;
#endif

  hypothesis->p3Vis_unit_ = visCandidate->p4().Vect().Unit();
  
  hypothesis->visMass_ = visCandidate->mass();

  // Add protection against zero mass:  
  // if lower than the electron mass, set it to the electron mass
  if ( hypothesis->visMass_ < 5.1e-4 ) {
    hypothesis->visMass_ = 5.1e-4;
  }

  // Set tau lepton decay mode
  hypothesis->decayMode_ = getDecayMode(visCandidate);

  // If this is a leptonic tau decay, we need to setup the limits on the
  // neutrino system invariant mass parameter.
  // In case reconstructed mass of visible decay products exceeds 1.6 GeV,
  // assume measurement error and "truncate" @ 1.6 GeV.
  if ( !nuSystemIsMassless() ) {
    SVfitParameter* fitParameter = algorithm_->getFitParameter(idxFitParameter_nuInvMass_);
    assert(fitParameter);
    fitParameter->setUpperLimit(SVfit_namespace::tauLeptonMass - TMath::Min(hypothesis->visMass_, 1.6));
    fitParameter->setInitialValue(0.5*(fitParameter->LowerLimit() + fitParameter->UpperLimit()));
  }

  // Extract the associated tracks, and fit a vertex if possible.
  hypothesis->tracks_ = extractTracks(visCandidate);

  // Sort tracks by **decreasing** Pt and determine "leading" (highest Pt) track
  std::sort(hypothesis->tracks_.begin(), hypothesis->tracks_.end(), isHigherPt);

  unsigned idx = 0;
  for ( std::vector<const reco::Track*>::const_iterator track = hypothesis->tracks_.begin();
	track != hypothesis->tracks_.end(); ++track ) {
#ifdef SVFIT_DEBUG 
    if ( this->verbosity_ ) std::cout << "Track #" << idx << ": Pt = " << (*track)->pt() << ", eta = " << (*track)->eta() << ", phi = " << (*track)->phi() << std::endl;
#endif
    const reco::HitPattern& trackHitPattern = (*track)->hitPattern();
    if ( trackHitPattern.numberOfValidTrackerHits() >= (int)trackMinNumHits_ &&
	 trackHitPattern.numberOfValidPixelHits() >= (int)trackMinNumPixelHits_ &&
	 (*track)->normalizedChi2() < trackMaxChi2DoF_ &&
	 ((*track)->ptError()/(*track)->pt()) < trackMaxDeltaPoverP_ &&
	 (*track)->pt() > trackMinPt_ ) {
#ifdef SVFIT_DEBUG 
      if ( this->verbosity_ ) std::cout << " passes quality cuts." << std::endl;
#endif
      hypothesis->selTracks_.push_back(*track);
    } else {
#ifdef SVFIT_DEBUG 
      if ( this->verbosity_ ) std::cout << " FAILS quality cuts." << std::endl;
#endif
    }
    ++idx;
  }

  if ( hypothesis->selTracks_.size() >= 1 ) hypothesis->leadTrack_ = hypothesis->selTracks_.at(0);
  else hypothesis->leadTrack_ = 0;

  for ( std::vector<const reco::Track*>::const_iterator track = hypothesis->selTracks_.begin();
	track != hypothesis->selTracks_.end(); ++track ) {
    reco::TransientTrack trajectory = trackBuilder_->build(*track);
    hypothesis->selTrackTrajectories_.push_back(trajectory);
  }

  if ( hypothesis->selTrackTrajectories_.size() >= 1 ) hypothesis->leadTrackTrajectory_ = &hypothesis->selTrackTrajectories_.at(0);
  else hypothesis->leadTrackTrajectory_ = 0;

  if ( initializeToGen_ ) {
    bool isMatched = false;
    int idx = 0;
    for ( GenParticleView::const_iterator genParticle = genParticles_->begin();
	  genParticle != genParticles_->end() && !isMatched; ++genParticle ) {
#ifdef SVFIT_DEBUG 
      if ( verbosity_ ) {
	std::cout << "genParticle #" << idx << " (pdgId = " << genParticle->pdgId() << "):" 
		  << " Pt = " << genParticle->pt() << ", eta = " << genParticle->eta() << ", phi = " << genParticle->phi() 
		  << " (status = " << genParticle->status() << ")" << std::endl;
      }
#endif
      if ( TMath::Abs(genParticle->pdgId()) == 15 ) {
	const reco::GenParticle* genTau = &(*genParticle);
	std::vector<const reco::GenParticle*> genTauDecayProducts;
	findDaughters(genTau, genTauDecayProducts, -1);
	reco::Candidate::LorentzVector genVisP4 = getVisMomentum(genTauDecayProducts);
	double dR = deltaR(genParticle->p4(), visCandidate->p4());
#ifdef SVFIT_DEBUG 
	if ( verbosity_ ) std::cout << "dR = " << dR << std::endl;
#endif
	if ( dR < dRmatch_ ) {
	  std::string genTauDecayMode = getGenTauDecayMode(genTau);
#ifdef SVFIT_DEBUG 
	  if ( verbosity_ ) std::cout << "genTauDecayMode = " << genTauDecayMode << std::endl;
#endif
	  std::string recTauDecayMode;
	  if      ( getDecayMode(visCandidate) == reco::PFTauDecayMode::tauDecaysElectron ) recTauDecayMode = "electron";
	  else if ( getDecayMode(visCandidate) == reco::PFTauDecayMode::tauDecayMuon      ) recTauDecayMode = "muon";
	  else                                                                              recTauDecayMode = "hadronic";
#ifdef SVFIT_DEBUG 
	  if ( verbosity_ ) std::cout << "recTauDecayMode = " << recTauDecayMode << std::endl;
#endif
	  // CV: check if reconstructed and generated decay modes match
	  //    (either electron, muon or hadronic)
    	  if ( (genTauDecayMode == "electron" && 
	        getDecayMode(visCandidate) == reco::PFTauDecayMode::tauDecaysElectron) ||
	       (genTauDecayMode == "muon" && 
	        getDecayMode(visCandidate) == reco::PFTauDecayMode::tauDecayMuon) ||
  	       (!(genTauDecayMode == "electron" || genTauDecayMode == "muon") && 
	        !(getDecayMode(visCandidate) == reco::PFTauDecayMode::tauDecaysElectron || getDecayMode(visCandidate) == reco::PFTauDecayMode::tauDecayMuon)) ) {
	    const reco::Candidate::LorentzVector& genTauP4 = genTau->p4();
	    genVisEnFracX_ = genVisP4.energy()/genTauP4.energy();
	    genPhiLab_ = phiLabFromLabMomenta(genTauP4, genVisP4);
	    reco::Candidate::LorentzVector genInvisP4 = getInvisMomentum(genTauDecayProducts);
	    genNuInvMass_ = genInvisP4.mass();
	    // CV: Add protection against rounding errors
	    if ( genNuInvMass_ < 0. ) genNuInvMass_ = 1.e-4;
	    reco::Candidate::Point genProductionVertex = genTau->vertex();
	    genDecayVertex_ = getDecayVertex(genTau);
	    genDecayDistance_ = TMath::Sqrt((genDecayVertex_ - genProductionVertex).mag2());
	    genVisP4_ = genVisP4;
#ifdef SVFIT_DEBUG 
	    if ( verbosity_ ) {
	      std::cout << " found match." << std::endl;
	      std::cout << "initializing:" << std::endl;
	      std::cout << " genVisEnFracX = " << genVisEnFracX_ << std::endl;
	      std::cout << " genPhiLab = " << genPhiLab_ << std::endl;
	      std::cout << " genVisInvMass = " << genVisP4_.mass() << std::endl;
	      std::cout << " genNuInvMass = " << genNuInvMass_ << std::endl;
	      std::cout << " genDecayVertex: x = " << genDecayVertex_.x() << ", y = " << genDecayVertex_.y() << ", z = " << genDecayVertex_.z() << std::endl;
	      std::cout << " genDecayDistance = " << genDecayDistance_ << std::endl;
	      std::cout << " genVisP4: Pt = " << genVisP4_.pt() << ", eta = " << genVisP4_.eta() << ", phi = " << genVisP4_.phi() << std::endl;
	    }
#endif
	    isMatched = true;	    
	  } 
	}
      }
      ++idx;
    }

    if ( !isMatched )
      throw cms::Exception("SVfitTauDecayBuilder::initialize")
	<< "Failed to find gen. Tau matching reconstructed decay products:" 
	<< " Pt = " << visCandidate->pt() << ", eta = " << visCandidate->eta() << ", phi = " << visCandidate->phi() << " !!\n";
    
    if ( fixToGenVisP4_ ) {
      hypothesis->p4_         = genVisP4_;
      hypothesis->p3Vis_unit_ = genVisP4_.Vect().Unit();
      hypothesis->visMass_    = genVisP4_.mass();
    }

    if ( idxFitParameter_visEnFracX_  != -1 ) algorithm_->setFitParameterInitialValue(idxFitParameter_visEnFracX_, genVisEnFracX_);
    if ( idxFitParameter_phi_lab_     != -1 ) algorithm_->setFitParameterInitialValue(idxFitParameter_phi_lab_,    genPhiLab_);
    if ( idxFitParameter_visMass_     != -1 ) algorithm_->setFitParameterInitialValue(idxFitParameter_visMass_,    genVisMass_);
    if ( idxFitParameter_nuInvMass_   != -1 ) algorithm_->setFitParameterInitialValue(idxFitParameter_nuInvMass_,  genNuInvMass_);
  }
}

void
SVfitTauDecayBuilder::finalize(SVfitSingleParticleHypothesis* hypothesis) const
{
  //std::cout << "<SVfitTauDecayBuilder>:" << std::endl;
  //std::cout << " hypothesis " << hypothesis->name() << " #" << hypothesis->barcode() << ": " << hypothesis << std::endl;

  // Cast to the concrete tau decay hypothesis
  SVfitTauDecayHypothesis* hypothesis_T = dynamic_cast<SVfitTauDecayHypothesis*>(hypothesis);
  assert(hypothesis_T);

  const SVfitEventHypothesis* event = hypothesis_T->mother()->eventHypothesis();
  assert(event);

  AlgebraicVector3 eventVertexPos = event->eventVertexPos();
  const AlgebraicMatrix33& eventVertexCov = event->eventVertexCov();
#ifdef SVFIT_DEBUG 
  if ( verbosity_ >= 2 ) {
    printVector("eventVertexPos", eventVertexPos);
    printMatrix("eventVertexCov", eventVertexCov);
  }
#endif
  hypothesis_T->expectedDecayDistance_  = 0.;
  hypothesis_T->expectedDecayVertexPos_ = eventVertexPos;
  hypothesis_T->expectedDecayVertexCov_ = eventVertexCov;

  hypothesis_T->hasDecayVertexFit_ = false;
  if ( idxFitParameter_decayDistance_shift_ != -1 ) {
    if ( hypothesis_T->selTracks_.size() >= minNumTracksFit && fitDecayVertex_ ) {
      TransientVertex decayVertex = decayVertexFitAlgorithm_->fitSecondaryVertex(hypothesis_T->selTracks_);
      if ( decayVertex.isValid() ) {
	if ( fixRecToGenVertex_ ) {
	  hypothesis_T->reconstructedDecayVertexPos_(0) = genDecayVertex_.x();
	  hypothesis_T->reconstructedDecayVertexPos_(1) = genDecayVertex_.y();
	  hypothesis_T->reconstructedDecayVertexPos_(2) = genDecayVertex_.z();
	} else {
	  hypothesis_T->reconstructedDecayVertexPos_(0) = decayVertex.position().x();
	  hypothesis_T->reconstructedDecayVertexPos_(1) = decayVertex.position().y();
	  hypothesis_T->reconstructedDecayVertexPos_(2) = decayVertex.position().z();
	}
	hypothesis_T->reconstructedDecayVertexCov_ = decayVertex.positionError().matrix_new();
	hypothesis_T->hasDecayVertexFit_ = true;
      
	AlgebraicMatrix33 reconstructedDecayVertex_wrt_eventVertexCov;
	for ( unsigned iRow = 0; iRow < 3; ++iRow ) {
	  for ( unsigned iColumn = 0; iColumn < 3; ++iColumn ) {
	    reconstructedDecayVertex_wrt_eventVertexCov(iRow, iColumn) = hypothesis_T->reconstructedDecayVertexCov_(iRow, iColumn) + eventVertexCov(iRow, iColumn);
	  }
	}
	double EigenValue1, EigenValue2, EigenValue3;
	extractEigenValues(reconstructedDecayVertex_wrt_eventVertexCov, EigenValue1, EigenValue2, EigenValue3);
	if ( EigenValue1 < 1.e-4 ) EigenValue1 = 1.e-4; // CV: assume resolution of tau decay vertex reconstruction 
	                                                //     to be never better than 100 microns in tau direction
	double sigma = TMath::Sqrt(EigenValue1);
#ifdef SVFIT_DEBUG 
	if ( verbosity_ >= 2 ) {
	  printVector("recDecayVertexPos", hypothesis_T->reconstructedDecayVertexPos_);
	  printMatrix("recDecayVertexCov", hypothesis_T->reconstructedDecayVertexCov_);
	}
#endif
	double fitParameterLowerLimit = -5.*sigma;
	double fitParameterUpperLimit = +5.*sigma;
	if ( initializeToGen_ ) {
	  double recDecayDistance = TMath::Sqrt(norm2(hypothesis_T->reconstructedDecayVertexPos_ - eventVertexPos));
	  algorithm_->setFitParameterInitialValue(idxFitParameter_decayDistance_shift_, genDecayDistance_ - recDecayDistance);
	  fitParameterLowerLimit = TMath::Min(fitParameterLowerLimit, -3*sigma + (genDecayDistance_ - recDecayDistance));
	  fitParameterUpperLimit = TMath::Max(fitParameterUpperLimit, +3*sigma + (genDecayDistance_ - recDecayDistance));
	} else {
	  algorithm_->setFitParameterInitialValue(idxFitParameter_decayDistance_shift_, 0.);
	}
#ifdef SVFIT_DEBUG 
	if ( verbosity_ >= 2 ) {
	  std::cout << "restricting fitParameter #" << idxFitParameter_decayDistance_shift_ << " to range " << fitParameterLowerLimit << ".." << fitParameterUpperLimit << std::endl;
	}
#endif	
	algorithm_->setFitParameterLimit(idxFitParameter_decayDistance_shift_, fitParameterLowerLimit, fitParameterUpperLimit);
	algorithm_->setFitParameterStepSize(idxFitParameter_decayDistance_shift_, 0.25*sigma);            
      }
    }
    if ( !hypothesis_T->hasDecayVertexFit_ ) {
      hypothesis_T->reconstructedDecayVertexPos_(0) = 0.;
      hypothesis_T->reconstructedDecayVertexPos_(1) = 0.;
      hypothesis_T->reconstructedDecayVertexPos_(2) = 0.;
      for ( int iRow = 0; iRow < 3; ++iRow ) {
	for ( int iColumn = 0; iColumn < 3; ++iColumn ) {
	  hypothesis_T->reconstructedDecayVertexCov_(iRow, iColumn) = 1.e+3;
	}
      }
      
      double fitParameterLowerLimit =  0.;
      double fitParameterUpperLimit = 10.;
      if ( initializeToGen_ ) { 
	double enTau_lab  = genVisP4_.energy()/genVisEnFracX_;
	double pTau_lab = TMath::Sqrt(square(enTau_lab) - tauLeptonMass2);
	double a = (pTau_lab/tauLeptonMass)*cTauLifetime;
#ifdef SVFIT_DEBUG 
	if ( verbosity_ >= 2 ) std::cout << "a(gen) = " << a << " --> setting fitParameter #" << idxFitParameter_decayDistance_shift_ << " to " << (genDecayDistance_/a) << std::endl;
#endif
	algorithm_->setFitParameterInitialValue(idxFitParameter_decayDistance_shift_, genDecayDistance_/a);
	fitParameterUpperLimit = TMath::Max(fitParameterUpperLimit, 5. + (genDecayDistance_/a));
      } else {
	algorithm_->setFitParameterInitialValue(idxFitParameter_decayDistance_shift_, 1.);
      }
      algorithm_->setFitParameterLimit(idxFitParameter_decayDistance_shift_, fitParameterLowerLimit, fitParameterUpperLimit);
      algorithm_->setFitParameterStepSize(idxFitParameter_decayDistance_shift_, 0.25);
    }
  }
}

//------------------------------------------------------------------------------- BEGIN: FOR TESTING ONLY
double getXfromGJangle(double gjAngle_rf, double gjAngle_lab, double pVis_lab, double visMass)
{
  std::cout << "<getXfromGJangle>:" << std::endl;
  std::cout << " gjAngle_rf = " << gjAngle_rf << std::endl;
  std::cout << " gjAngle_lab = " << gjAngle_lab << std::endl;
  std::cout << " pVis_lab = " << pVis_lab << std::endl;
  std::cout << " visMass = " << visMass << std::endl;
  double pVis_lab2 = square(pVis_lab);
  double visMass2 = square(visMass);
  double enVis_lab = TMath::Sqrt(pVis_lab2 + visMass2);
  double enVis_lab2 = square(enVis_lab);
  double sign = ( gjAngle_rf > (0.5*TMath::Pi()) ) ? +1. : -1.; // CV: for fixed pVis, the 'forward' solution is the one corresponding to **lower** pTau !!
  std::cout << "sign = " << sign << std::endl;
  double pTau_lab = ((visMass2 + tauLeptonMass2)*pVis_lab*TMath::Cos(gjAngle_lab) 
		     + sign*TMath::Sqrt(enVis_lab2*(square(visMass2 - tauLeptonMass2) - 4.*tauLeptonMass2*pVis_lab2*square(TMath::Sin(gjAngle_lab)))))/
    (2.*(visMass2 + pVis_lab2*square(TMath::Sin(gjAngle_lab))));
  std::cout << "pTau_lab = " << pTau_lab << std::endl;
  double X = TMath::Sqrt((pVis_lab2 + visMass2)/(square(pTau_lab) + tauLeptonMass2));
  std::cout << "X = " << X << std::endl;
  return X;
}

void runConsistencyCheck(double gjAngle_rf, double gjAngle_lab_input, double pVis_lab, double visMass, double invisMass)
{
  std::cout << "<runConsistencyCheck>:" << std::endl;
  std::cout << " gjAngle_rf = " << gjAngle_rf << std::endl;
  std::cout << " gjAngle_lab (input) = " << gjAngle_lab_input << std::endl;
  std::cout << " pVis_lab = " << pVis_lab << std::endl;
  std::cout << " visMass = " << visMass << std::endl;

  double X = getXfromGJangle(gjAngle_rf, gjAngle_lab_input, pVis_lab, visMass);

  double enVis_lab = TMath::Sqrt(square(pVis_lab) + square(visMass));
  std::cout << " enVis_lab = " << enVis_lab << std::endl;

  bool isValidSolution = true;
  double gjAngle_lab_output = gjAngleFromX_new(X, visMass, invisMass, pVis_lab, enVis_lab, tauLeptonMass, isValidSolution);
  assert(isValidSolution);
  std::cout << " gjAngle_lab (output) = " << gjAngle_lab_output << std::endl;

  const double epsilon = 1.e-6;
       if ( TMath::Abs(gjAngle_lab_output - gjAngle_lab_input) < epsilon ) std::cout << "CHECK passed." << std::endl;
  else std::cout << "CHECK failed !!" << std::endl;
}
//------------------------------------------------------------------------------- END: FOR TESTING ONLY

bool SVfitTauDecayBuilder::applyFitParameter(SVfitSingleParticleHypothesis* hypothesis, const double* param) const
{
#ifdef SVFIT_DEBUG 
  if ( verbosity_ ) {
    std::cout << "<SVfitTauDecayBuilder::applyFitParameter>:" << std::endl;
  }
#endif

  // Cast to the concrete tau decay hypothesis
  SVfitTauDecayHypothesis* hypothesis_T = dynamic_cast<SVfitTauDecayHypothesis*>(hypothesis);
  assert(hypothesis_T);
  
  if ( idxFitParameter_visMass_  != -1 ) {
    if ( fixToGenVisMass_ ) hypothesis_T->visMass_ = genVisMass_;
    else hypothesis_T->visMass_ = param[idxFitParameter_visMass_];
  }

  double visEnFracX = param[idxFitParameter_visEnFracX_];
  double phi_lab    = param[idxFitParameter_phi_lab_];
  double pVis_lab   = hypothesis_T->p4().P();
  double enVis_lab  = hypothesis_T->p4().energy();
  double visMass    = hypothesis_T->visMass();
  double nuInvMass  = ( nuSystemIsMassless() ) ? 0. : param[idxFitParameter_nuInvMass_];

  if ( fixToGenVisEnFracX_ ) visEnFracX = genVisEnFracX_;
  if ( fixToGenPhiLab_     ) phi_lab    = genPhiLab_;
  if ( fixToGenNuInvMass_  ) nuInvMass  = genNuInvMass_;

  const reco::Candidate::Vector& p3Vis_unit = hypothesis_T->p3Vis_unit();
  
  bool isValidSolution = true;

/*
//--- compute energy and momentum of visible decay products in tau lepton rest-frame
  double pVis2_rf = (tauLeptonMass2 - square(visMass + nuInvMass))*(tauLeptonMass2 - square(visMass - nuInvMass))/(4.*tauLeptonMass2);
  if ( pVis2_rf < 0. ) {
    pVis2_rf = 0.;
    isValidSolution = false;
  }
  double pVis_rf = TMath::Sqrt(pVis2_rf);
  double enVis_rf = TMath::Sqrt(square(pVis_rf) + square(visMass));
  //std::cout << "pVis_rf = " << pVis_rf << ", enVis_rf = " << enVis_rf << std::endl;

//--- compute relativistic beta factor
  double beta2 = 1. - square(visEnFracX*tauLeptonMass/enVis_lab);
  if ( !(beta2 >= 0.) ) {
    std::cout << "<SVfitTauDecayBuilder::applyFitParameter>:" << std::endl;
    std::cout << " leg: Pt = " << hypothesis_T->particle()->pt() << "," 
	      << " eta = " << hypothesis_T->particle()->eta() << "," 
	      << " phi = " << hypothesis_T->particle()->phi() << ","
	      << " mass = " << hypothesis_T->particle()->mass() << std::endl;
    std::cout << " visEnFracX = " << visEnFracX << std::endl;
    std::cout << " tauLeptonMass = " << tauLeptonMass << std::endl;
    std::cout << " enVis_lab = " << enVis_lab << std::endl;
    std::cout << "(fixToGenVisEnFracX = " << fixToGenVisEnFracX_ << ")" << std::endl;
    std::cout << "(fixToGenPhiLab = " << fixToGenPhiLab_ << ")" << std::endl;
    std::cout << "(fixToGenVisMass = " << fixToGenVisMass_ << ")" << std::endl;
    std::cout << "(fixToGenNuInvMass = " << fixToGenNuInvMass_ << ")" << std::endl;
    std::cout << "(fixToGenDecayDistance = " << fixToGenDecayDistance_ << ")" << std::endl;
    std::cout << " assert TRIGGERED !!!" << std::endl;
  }
  assert(beta2 >= 0.);
  double beta = TMath::Sqrt(beta2);
  //std::cout << "beta = " << beta << std::endl;

//--- compute Gottfried-Jackson angle
//   (angle of visible decay products wrt. tau lepton flight direction)
  double cosGjAngle_rf = (visEnFracX*tauLeptonMass - enVis_rf)/(beta*pVis_rf);  
  if ( cosGjAngle_rf < -1. ) {
#ifdef SVFIT_DEBUG 
    if ( verbosity_ >= 2 ) std::cout << "cosGjAngle_rf = " << cosGjAngle_rf << " --> setting isValidSolution = false." << std::endl;
#endif
    cosGjAngle_rf = -1.;
    isValidSolution = false;
  } else if ( cosGjAngle_rf > +1. ) {
#ifdef SVFIT_DEBUG 
    if ( verbosity_ >= 2 ) std::cout << "cosGjAngle_rf = " << cosGjAngle_rf << " --> setting isValidSolution = false." << std::endl;
#endif
    cosGjAngle_rf = +1.;
    isValidSolution = false;
  }
  double gjAngle_rf = TMath::ACos(cosGjAngle_rf);
  //std::cout << "gjAngle_rf = " << gjAngle_rf << std::endl;

//--- compute component of visible momentum perpendicular to tau flight direction
//   (NB: the perpendicular component of the visible momentum vector
//        is invariant under Lorentz boost in tau direction, 
//        i.e. is identical in tau rest-frame and laboratory frame)
  double pVis_perp = pVis_rf*TMath::Sin(gjAngle_rf);
  //std::cout << "pVis_perp = " << pVis_perp << std::endl;

  double gjAngle_lab = TMath::ASin(pVis_perp/pVis_lab);
  //std::cout << "gjAngle_lab = " << gjAngle_lab << std::endl;
 */
  
//--- compute Gottfried-Jackson angle
//   (angle of visible decay products wrt. tau lepton flight direction)
  double gjAngle_lab = gjAngleFromX_new(visEnFracX, visMass, nuInvMass, pVis_lab, enVis_lab, tauLeptonMass, isValidSolution);
#ifdef SVFIT_DEBUG 
  if ( verbosity_ >= 2 ) std::cout << "gjAngle_lab = " << gjAngle_lab << std::endl;
#endif

/*
//--- compute energy and momentum of visible decay products in tau lepton rest-frame
  double pVis2_rf = (tauLeptonMass2 - square(visMass + nuInvMass))*(tauLeptonMass2 - square(visMass - nuInvMass))/(4.*tauLeptonMass2);
  if ( pVis2_rf < 0. ) {
    pVis2_rf = 0.;
    isValidSolution = false;
  }
  double pVis_rf = TMath::Sqrt(pVis2_rf);
#ifdef SVFIT_DEBUG 
  if ( verbosity_ >= 2 ) std::cout << "pVis_rf = " << pVis_rf << std::endl;
#endif

//--- compute component of visible momentum perpendicular to tau flight direction
//   (NB: the perpendicular component of the visible momentum vector
//        is invariant under Lorentz boost in tau direction, 
//        i.e. is identical in tau rest-frame and laboratory frame)
  double pVis_perp = pVis_lab*TMath::Sin(gjAngle_lab);
#ifdef SVFIT_DEBUG 
  if ( verbosity_ >= 2 ) std::cout << "pVis_perp = " << pVis_perp << std::endl;
#endif
    
  double gamma_tau = Etau_/Mtau_;
    std::cout << "gamma_tau = " << gamma_tau << std::endl;
    double beta_tau = TMath::Sqrt(1. - 1./(gamma_tau*gamma_tau));
    std::cout << "beta_tau = " << beta_tau << std::endl;
    double pVis_parl_rf = -beta_tau*gamma_tau*Ehad + gamma_tau*TMath::Cos(alpha)*Phad;
    std::cout << "pVis_parl_rf = " << pVis_parl_rf << std::endl;
    double alpha_rf = TMath::ATan2(pVis_perp, pVis_parl_rf);

  // CV: the ratio pVis_perp/pVis_rf determines only sin(gjAngle_rf).
  //     Distinguish 'forward' from 'backward' decays depending on whether the visible tau decay products
  //     carry more or less energy than a tau with maximum Gottfried-Jackson angle 
  //    (cf. Fig. 4.3 of AN-2010/256:
  //       http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2010_256_v2.pdf
  //     for illustration)
    //  double gjAngle_lab_max = TMath::ASin((tauLeptonMass2 - square(visMass))/(2.*tauLeptonMass)

    //as function of X = Etau/Evis
    //     was obtained by solving equation (1) of AN-2010/256:
    //       http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2010_256_v2.pdf
    //     for cosThetaGJ

the visible tau decay products
  //     carry more or or less energy

  //     Distinguish 'forward' from 'backward' decays depending on whether the visible tau decay products
  //     carry more
energy fraction X = Etau/Evis
  //     is greater or smaller than 50%.
  double gjAngle_rf = ( visEnFracX > 0.5 ) ?
    TMath::ASin(pVis_perp/pVis_rf) : (TMath::Pi() - TMath::ASin(pVis_perp/pVis_rf));
#ifdef SVFIT_DEBUG 
  if ( verbosity_ >= 2 ) std::cout << "gjAngle_rf = " << gjAngle_rf << std::endl;
#endif
 */

//--- compute tau lepton energy and momentum in laboratory frame  
  double enTau_lab = enVis_lab/visEnFracX;
  double pTau_lab = TMath::Sqrt(square(enTau_lab) - tauLeptonMass2);
#ifdef SVFIT_DEBUG 
  if ( verbosity_ >= 2 ) std::cout << "pTau_lab = " << pTau_lab << ", enTau_lab = " << enTau_lab << std::endl;
#endif

//--- compute Gottfried-Jackson angle in tau lepton rest frame:
//    the component of the visible momentum vector parallel to the tau flight direction
//    is computed via a Lorentz boost from the laboratory frame to the tau rest-frame;
//    the perpendicular component of the visible momentum vector 
//    is invariant under Lorentz boost in tau flight direction, 
//    i.e. is identical in tau rest-frame and laboratory frame.
  double gamma = enTau_lab/tauLeptonMass;
  //std::cout << "gamma = " << gamma << std::endl;
  double beta = TMath::Sqrt(1. - 1./(gamma*gamma));
  //std::cout << "beta = " << beta << std::endl;
  double pVis_parl_rf = -beta*gamma*enVis_lab + gamma*TMath::Cos(gjAngle_lab)*pVis_lab;
  double pVis_perp = pVis_lab*TMath::Sin(gjAngle_lab);
  double gjAngle_rf = TMath::ATan2(pVis_perp, pVis_parl_rf);
#ifdef SVFIT_DEBUG 
  if ( verbosity_ >= 2 ) {
    std::cout << "pVis_parl_rf = " << pVis_parl_rf << ", pVis_perp = " << pVis_perp << " --> gjAngle_rf = " << gjAngle_rf << std::endl;
  }
#endif

  reco::Candidate::Vector p3Tau_unit = motherDirection(p3Vis_unit, gjAngle_lab, phi_lab);
  reco::Candidate::LorentzVector p4Tau_lab = reco::Candidate::LorentzVector(
    pTau_lab*p3Tau_unit.x(), pTau_lab*p3Tau_unit.y(), pTau_lab*p3Tau_unit.z(), enTau_lab);

  //double phi_labDEBUG = phiLabFromLabMomenta(p4Tau_lab, hypothesis_T->p4_);
  //std::cout << "phi_labDEBUG = " << phi_labDEBUG << std::endl;

  hypothesis_T->p4_fitted_  = p4Tau_lab;
  hypothesis_T->dp4_        = (p4Tau_lab - hypothesis_T->p4_);

  reco::Candidate::Vector boostToTauRF_vector = p4Tau_lab.BoostToCM();
  hypothesis_T->p4invis_rf_ = ROOT::Math::VectorUtil::boost(hypothesis_T->dp4_, boostToTauRF_vector);
  hypothesis_T->p4vis_rf_   = ROOT::Math::VectorUtil::boost(hypothesis_T->p4(), boostToTauRF_vector);

  hypothesis_T->visEnFracX_ = visEnFracX;
  hypothesis_T->gjAngle_    = gjAngle_rf;
  hypothesis_T->phiLab_     = phi_lab;
  
//--- update track based information
  if ( idxFitParameter_decayDistance_shift_ != -1 ) {
    const SVfitEventHypothesis* event = hypothesis_T->mother()->eventHypothesis();
    AlgebraicVector3 eventVertexPos = event->eventVertexPos();
#ifdef SVFIT_DEBUG 
    if ( verbosity_ >= 2 ) printVector(" eventVertexPos", eventVertexPos);
#endif    
    hypothesis_T->expectedFlightPath_unit_ = AlgebraicVector3(p3Tau_unit.x(), p3Tau_unit.y(), p3Tau_unit.z());
    hypothesis_T->expectedDecayDistanceJacobiFactor_ = 0.;
    if ( fixToGenDecayDistance_ ) {
      hypothesis_T->expectedDecayDistance_ = genDecayDistance_;
      hypothesis_T->expectedDecayDistanceJacobiFactor_ = 1.;
    } else {
      if ( hypothesis_T->hasDecayVertexFit_ ) {      
	double decayDistance0 = TMath::Sqrt(norm2(hypothesis_T->reconstructedDecayVertexPos_ - eventVertexPos));
	//std::cout << "decayDistance0 = " << decayDistance0 << std::endl;
	double decayDistance_shift = param[idxFitParameter_decayDistance_shift_];
	//std::cout << "decayDistance_shift = " << decayDistance_shift << std::endl;
	hypothesis_T->expectedDecayDistance_ = decayDistance0 + decayDistance_shift;
	// CV: expectedDecayDistanceJacobiFactor accounts for Jacobi determinant 
	//       | ddecayDistance/dx dphi_lab/dx |
	//       | ddecayDistance/dy dphi_lab/dy |
	//     that arises from the fact that the integration variables are decayDistance, phi_lab
	//     while the likelihood-function implemented in SVfitTauDecayLikelihoodTrackInfo uses x, y coordinates;
	//     c.f. 
	//       https://de.wikipedia.org/wiki/Delta-Distribution (section "Delta-Distribution in krummlinigen Koordinatensystemen")
	//     Sorry, I did not find this section in the English wikipedia.
	//
	//    (using tau lifetime information in the likelihood is numerically stable only if SVfit is run in "integration" mode,
	//     so it is a safe assumption that the "fit" mode, for which the Jacobi factor may not apply, is not used)
	//
	//hypothesis_T->expectedDecayDistanceJacobiFactor_ = hypothesis_T->expectedDecayDistance_;
	hypothesis_T->expectedDecayDistanceJacobiFactor_ = 1.0;
      } else if ( hypothesis_T->leadTrackTrajectory_ ) {
	double a = (pTau_lab/tauLeptonMass)*cTauLifetime;
	// CV: expectedDecayDistanceJacobiFactor accounts for:
	//    1) Jacobi determinant 
	//        | ddecayDistance/dx dphi_lab/dx |
	//        | ddecayDistance/dy dphi_lab/dy |
	//      that arises from the fact that the integration variables are decayDistance, phi_lab
	//      while the likelihood-function implemented in SVfitTauDecayLikelihoodTrackInfo uses x, y coordinates;
	//      c.f. 
	//        https://de.wikipedia.org/wiki/Delta-Distribution (section "Delta-Distribution in krummlinigen Koordinatensystemen")
	//      Sorry, I did not find this section in the English wikipedia.
	//
	//     (using tau lifetime information in the likelihood is numerically stable only if SVfit is run in "integration" mode,
	//      so it is a safe assumption that the "fit" mode, for which the Jacobi factor may not apply, is not used)	
	//
	//    2) parametrization of the decay distance in terms of "expected decay distance" = beta gamma*c tau = (pTau/mTau)*c tau
	//
	//hypothesis_T->expectedDecayDistanceJacobiFactor_ = hypothesis_T->expectedDecayDistance_*a; 
	hypothesis_T->expectedDecayDistanceJacobiFactor_ = a;
	hypothesis_T->expectedDecayDistance_ = param[idxFitParameter_decayDistance_shift_]*a;
      }
    }
    if ( hypothesis_T->expectedDecayDistance_ < 0. || hypothesis_T->expectedDecayDistance_ > 25. ) {
#ifdef SVFIT_DEBUG 
      if ( verbosity_ >= 2 ) std::cout << "expectedDecayDistance = " << hypothesis_T->expectedDecayDistance_ << " --> setting isValidSolution = false." << std::endl;
#endif    
      hypothesis_T->expectedDecayDistance_ = 0.;
      isValidSolution = false;
    }
    
    if ( isValidSolution ) {
      //hypothesis_T->expectedDecayVertexPos_ = compDecayPosition_line(eventVertexPos, hypothesis_T->expectedDecayDistance_, hypothesis_T->expectedFlightPath_unit_);
      hypothesis_T->expectedDecayVertexPos_ = compDecayPosition_helix(eventVertexPos, hypothesis_T->expectedDecayDistance_, p4Tau_lab, hypothesis_T->particle()->charge());

      if ( hypothesis_T->leadTrackTrajectory_ && !hypothesis_T->hasDecayVertexFit_ ) {      	
	SVfitTrackExtrapolation leadTrack_extrapolation(*hypothesis_T->leadTrackTrajectory_, hypothesis_T->expectedDecayVertexPos_);
	hypothesis_T->reconstructedDecayVertexPos_ = leadTrack_extrapolation.point_of_closest_approach();
	hypothesis_T->reconstructedDecayVertexCov_ = leadTrack_extrapolation.covariance();
	hypothesis_T->leadTrackExtrapolationError_ = leadTrack_extrapolation.errorFlag();
#ifdef SVFIT_DEBUG 
	if ( verbosity_ >= 2  ) {
	  std::cout << "SVfitTrackExtrapolation:" << std::endl;
	  printVector("reconstructedDecayVertexPos", hypothesis_T->reconstructedDecayVertexPos_);
	  printMatrix("reconstructedDecayVertexCov", hypothesis_T->reconstructedDecayVertexCov_);
	}
#endif    
      }
    }
  }
#ifdef SVFIT_DEBUG 
  if ( verbosity_ ) {
    std::cout << "<SVfitTauDecayBuilder::applyFitParameter>:" << std::endl;
    std::cout << " hypothesis " << hypothesis->name() << " #" << hypothesis->barcode() << ": " << hypothesis << std::endl;
    std::cout << " visEnFracX = " << visEnFracX << " (getXfromGJangle = " << getXfromGJangle(gjAngle_rf, gjAngle_lab, pVis_lab, visMass) << ")" << std::endl;
/*
std::cout << "getXfromGJangle(0.1, 0.005, 25., 1.2) = " << getXfromGJangle(0.1, 0.005, 25., 1.2) << std::endl;
std::cout << "getXfromGJangle(3.0, 0.005, 25., 1.2) = " << getXfromGJangle(3.0, 0.005, 25., 1.2) << std::endl;
std::cout << "getXfromGJangle(0.1, 0.015, 25., 1.2) = " << getXfromGJangle(0.1, 0.015, 25., 1.2) << std::endl;
std::cout << "getXfromGJangle(3.0, 0.015, 25., 1.2) = " << getXfromGJangle(3.0, 0.015, 25., 1.2) << std::endl;
std::cout << "getXfromGJangle(0.1, 0.005, 40., 1.2) = " << getXfromGJangle(0.1, 0.005, 40., 1.2) << std::endl;
std::cout << "getXfromGJangle(3.0, 0.005, 40., 1.2) = " << getXfromGJangle(3.0, 0.005, 40., 1.2) << std::endl;
std::cout << "getXfromGJangle(0.1, 0.010, 40., 1.2) = " << getXfromGJangle(0.1, 0.010, 40., 1.2) << std::endl;
std::cout << "getXfromGJangle(3.0, 0.010, 40., 1.2) = " << getXfromGJangle(3.0, 0.010, 40., 1.2) << std::endl;
assert(0);
 */
/*
    runConsistencyCheck(0.1, 0.005, 25., 1.2, 0.);
    runConsistencyCheck(3.0, 0.005, 25., 1.2, 0.);
    runConsistencyCheck(0.1, 0.015, 25., 1.2, 0.);
    runConsistencyCheck(3.0, 0.015, 25., 1.2, 0.);
    runConsistencyCheck(0.1, 0.005, 40., 1.2, 0.);
    runConsistencyCheck(3.0, 0.005, 40., 1.2, 0.);
    runConsistencyCheck(0.1, 0.010, 40., 1.2, 0.);
    runConsistencyCheck(3.0, 0.010, 40., 1.2, 0.);    
    assert(0);
 */ 
    std::cout << " phi_lab = " << phi_lab << std::endl;
    std::cout << " visMass = " << visMass << std::endl;
    std::cout << " nuInvMass = " << hypothesis_T->p4invis_rf_.mass() << std::endl;
    std::cout << " gjAngle: rf = " << gjAngle_rf << ", lab = " << gjAngle_lab << std::endl;
    if ( idxFitParameter_decayDistance_shift_ != -1 ) {
      std::cout << "decayVertex: x = " << hypothesis_T->expectedDecayVertexPos_(0) << "," 
		<< " y = " << hypothesis_T->expectedDecayVertexPos_(1) << "," 
		<< " z = " << hypothesis_T->expectedDecayVertexPos_(2) << " (d = " << hypothesis_T->expectedDecayDistance_ << ")" << std::endl;
      AlgebraicVector3 eventVertexPos = hypothesis_T->mother()->eventHypothesis()->eventVertexPos();
      double dxFlight = hypothesis_T->expectedDecayVertexPos_(0) - eventVertexPos(0);
      double dyFlight = hypothesis_T->expectedDecayVertexPos_(1) - eventVertexPos(1);
      double dzFlight = hypothesis_T->expectedDecayVertexPos_(2) - eventVertexPos(2);
      double dFlight  = TMath::Sqrt(square(dxFlight) + square(dyFlight) + square(dzFlight));
      reco::Candidate::LorentzVector p4Flight(dxFlight, dyFlight, dzFlight, dFlight);
      std::cout << " flight-path: dx = " << dxFlight << ", dy = " << dyFlight << ", dz = " << dzFlight 
		<< " (eta = " << p4Flight.eta() << ", phi = " << p4Flight.phi() << ", distance = " << dFlight << ","
		<< " angle wrt. vis = " << angle(p4Flight, hypothesis_T->p4_) << ")" << std::endl;  
    }
    std::cout << "p4Vis (lab): E = " << hypothesis_T->p4_.energy() << ","
	      << " px = " << hypothesis_T->p4_.px() << ", py = " << hypothesis_T->p4_.py() << ","
	      << " pz = " << hypothesis_T->p4_.pz() << " (eta = " << hypothesis_T->p4_.eta() << ", phi = " << hypothesis_T->p4_.phi() << ")" << std::endl;
    std::cout << "chargeTau = " << hypothesis_T->particle()->charge() << std::endl;
    std::cout << "p4Tau (lab): E = " << p4Tau_lab.energy() << ","
	      << " px = " << p4Tau_lab.px() << ", py = " << p4Tau_lab.py() << ","
	      << " pz = " << p4Tau_lab.pz() << " (eta = " << p4Tau_lab.eta() << ", phi = " << p4Tau_lab.phi() << ")" << std::endl;
    std::cout << "p4Vis (rf): E = " << hypothesis_T->p4vis_rf_.energy() << ","
	      << " mass = " << hypothesis_T->p4vis_rf_.mass()
	      << " px = " << hypothesis_T->p4vis_rf_.px() << ", py = " << hypothesis_T->p4vis_rf_.py() << ","
	      << " pz = " << hypothesis_T->p4vis_rf_.pz() << std::endl;
    std::cout << "p4Invis (rf): E = " << hypothesis_T->p4invis_rf_.energy() << ","
              << " mass = " << hypothesis_T->p4invis_rf_.mass()
	      << " px = " << hypothesis_T->p4invis_rf_.px() << ", py = " << hypothesis_T->p4invis_rf_.py() << ","
	      << " pz = " << hypothesis_T->p4invis_rf_.pz() << std::endl;
    std::cout << "isValidSolution = " << isValidSolution << std::endl;
  }
#endif    
  hypothesis_T->isValidSolution_ = isValidSolution;

  return isValidSolution;
}

void SVfitTauDecayBuilder::print(std::ostream& stream) const
{
  stream << "<SVfitTauDecayBuilderBase::print>:" << std::endl;
  stream << " pluginName = " << pluginName_ << std::endl;
  stream << " pluginType = " << pluginType_ << std::endl;
  stream << " prodParticleLabel = " << prodParticleLabel_ << std::endl;
  stream << " idxFitParameter_visEnFracX = " << idxFitParameter_visEnFracX_ << std::endl;
  stream << " idxFitParameter_phi_lab = " << idxFitParameter_phi_lab_ << std::endl;
  stream << " idxFitParameter_visMass = " << idxFitParameter_visMass_ << std::endl;
  stream << " idxFitParameter_nuInvMass = " << idxFitParameter_nuInvMass_ << std::endl;
  stream << " idxFitParameter_decayDistance_shift = " << idxFitParameter_decayDistance_shift_ << std::endl;
}

//
//-------------------------------------------------------------------------------
//

void applyOptionalFitParameter(const double* param, int idxFitParameter, double& value)
{
  if   ( idxFitParameter != -1 ) value = param[idxFitParameter];
  else                           value = 0.;
}

