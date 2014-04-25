#include "TauAnalysis/SVfit/plugins/SVfitTauToHadLikelihoodPhaseSpace.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "TauAnalysis/SVfit/interface/SVfitAlgorithmBase.h"
#include "TauAnalysis/SVfit/interface/svFitAuxFunctions.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitTauToHadHypothesis.h"

#include <TFile.h>
#include <TMath.h>

using namespace svFit_namespace;

#define SVFIT_DEBUG 1

SVfitTauToHadLikelihoodPhaseSpace::SVfitTauToHadLikelihoodPhaseSpace(const edm::ParameterSet& cfg)
  : SVfitSingleParticleLikelihood(cfg),
    histogram_(0),
    algorithm_(0)
{
  applySinThetaFactor_ = cfg.exists("applySinThetaFactor") ?
    cfg.getParameter<bool>("applySinThetaFactor") : false;

  varyVisMass_ = cfg.getParameter<bool>("varyVisMass");
  if ( varyVisMass_ ) {
    edm::FileInPath inputFileName = cfg.getParameter<edm::FileInPath>("inputFileName");
    if ( !inputFileName.isLocal() ) 
      throw cms::Exception("SVfitTauToHadLikelihoodPhaseSpace") 
	<< " Failed to find File = " << inputFileName << " !!\n";
    std::string histogramName = cfg.getParameter<std::string>("histogramName");
    TFile* inputFile = new TFile(inputFileName.fullPath().data());
    TH1* histogramTmp = dynamic_cast<TH1*>(inputFile->Get(histogramName.data()));
    if ( !histogramTmp )
      throw cms::Exception("SVfitTauToHadLikelihoodPhaseSpace") 
	<< " Failed to load visMassHistogram = " << histogramName.data() << " from file = " << inputFileName.fullPath().data() << " !!\n";
    histogram_ = (TH1*)histogramTmp->Clone(std::string(pluginName_).append("_").append(histogramTmp->GetName()).data());
    firstBin_ = 1;
    lastBin_ = histogram_->GetNbinsX();
    delete inputFile;
  }
  varyDeltaVisMass_ = cfg.getParameter<bool>("varyDeltaVisMass");
  if ( varyDeltaVisMass_ ) {

  }
}

SVfitTauToHadLikelihoodPhaseSpace::~SVfitTauToHadLikelihoodPhaseSpace()
{
  delete histogram_;
}

void SVfitTauToHadLikelihoodPhaseSpace::beginJob(SVfitAlgorithmBase* algorithm)
{
  algorithm_ = algorithm;

  algorithm->requestFitParameter(prodParticleLabel_,   svFit_namespace::kTau_visEnFracX, pluginName_);
  algorithm->requestFitParameter(prodParticleLabel_,   svFit_namespace::kTau_phi_lab,    pluginName_);
  if ( varyVisMass_ || varyDeltaVisMass_ ) {
    algorithm->requestFitParameter(prodParticleLabel_, svFit_namespace::kTau_visMass,    pluginName_);
  }
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
  if ( !(varyVisMass_ || varyDeltaVisMass_) ) {
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
    int bin = histogram_->FindBin(visMass);
    if ( bin <= firstBin_ ) bin = firstBin_;
    if ( bin >= lastBin_  ) bin = lastBin_;
    prob *= histogram_->GetBinContent(bin);
  }
  if ( varyDeltaVisMass_ ) {    

  }
  
  if ( applyVisPtCutCorrection_ ) {
    double probCorr = 1.;
    const double epsilon_regularization = 1.e-1;
    if ( hypothesis_T->p4_fitted().pt() > visPtCutThreshold_ ) {     
      double xCut = visPtCutThreshold_/hypothesis_T->p4_fitted().pt();      
      probCorr = 1./((1. - xCut) + epsilon_regularization);
    }
#ifdef SVFIT_DEBUG 
    if ( this->verbosity_ ) std::cout << "probCorr (had) = " << probCorr << std::endl;
#endif
    prob *= probCorr;
  }

  if ( algorithm_->applyJacobiFactors() && visEnFracX > 0. ) {
    double jacobiFactor = 1./(cube(visEnFracX));
#ifdef SVFIT_DEBUG 
    if ( this->verbosity_ ) std::cout << "jacobiFactor = " << jacobiFactor << std::endl;
#endif
    prob *= jacobiFactor;
  }

#ifdef SVFIT_DEBUG       
  if ( this->verbosity_ ) std::cout << "--> prob = " << prob << std::endl;
#endif
  return prob;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(SVfitSingleParticleLikelihoodPluginFactory, SVfitTauToHadLikelihoodPhaseSpace, "SVfitTauToHadLikelihoodPhaseSpace");
