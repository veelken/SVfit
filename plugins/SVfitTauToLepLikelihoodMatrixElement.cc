#include "TauAnalysis/SVfit/plugins/SVfitTauToLepLikelihoodMatrixElement.h"

#include "TauAnalysis/SVfit/interface/SVfitAlgorithmBase.h"
#include "TauAnalysis/SVfit/interface/svFitAuxFunctions.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitTauToDaughtersHypothesisBaseT1T2.h"

#include <TMath.h>

#include <limits>

using namespace svFit_namespace;

template <typename T>
SVfitTauToLepLikelihoodMatrixElement<T>::SVfitTauToLepLikelihoodMatrixElement(const edm::ParameterSet& cfg)
  : SVfitSingleParticleLikelihood(cfg),
    algorithm_(0)
{
  applySinThetaFactor_ = cfg.exists("applySinThetaFactor") ?
    cfg.getParameter<bool>("applySinThetaFactor") : false;
}

template <typename T>
SVfitTauToLepLikelihoodMatrixElement<T>::~SVfitTauToLepLikelihoodMatrixElement()
{
// nothing to be done yet...
}

template <typename T>
void SVfitTauToLepLikelihoodMatrixElement<T>::beginJob(SVfitAlgorithmBase* algorithm)
{
  algorithm_ = algorithm;

  algorithm->requestFitParameter(prodParticleLabel_, svFit_namespace::kTau_visEnFracX, pluginName_);
  algorithm->requestFitParameter(prodParticleLabel_, svFit_namespace::kTau_phi_lab,    pluginName_);
  algorithm->requestFitParameter(prodParticleLabel_, svFit_namespace::kTau_nuInvMass,  pluginName_);
}

template <typename T>
double SVfitTauToLepLikelihoodMatrixElement<T>::operator()(const SVfitSingleParticleHypothesis* hypothesis, int polSign) const
{
//--- compute negative log-likelihood for tau lepton decay 
//    tau- --> e- nu nu (tau- --> mu- nu nu)
//    to be compatible with matrix element of V-A electroweak decay
//   (ignoring tau lepton polarization effects)
//
//    NOTE: The formulas taken from the paper
//           "Tau polarization and its correlations as a probe of new physics",
//           B.K. Bullock, K. Hagiwara and A.D. Martin,
//           Nucl. Phys. B395 (1993) 499.
//
  const SVfitTauToDaughtersHypothesisBaseT1T2<SVfitTauDecayHypothesis, T>* hypothesis_T = 
    dynamic_cast<const SVfitTauToDaughtersHypothesisBaseT1T2<SVfitTauDecayHypothesis, T>*>(hypothesis);
  assert(hypothesis_T != 0);

  double decayAngle = hypothesis_T->gjAngle();
  double nuMass = hypothesis_T->p4invis_rf().mass();
  if ( nuMass < 0. ) nuMass = 0.; // CV: add protection against rounding errors when boosting between laboratory and rest frame
  double nuMass2 = square(nuMass);
  double visEnFracX = hypothesis_T->visEnFracX();
#ifdef SVFIT_DEBUG 
  if ( this->verbosity_ ) {
    std::cout << "<SVfitTauToLepLikelihoodMatrixElement::operator()>:" << std::endl;
    std::cout << " decayAngle = " << decayAngle << std::endl;
    std::cout << " nuMass = " << nuMass << std::endl;
    std::cout << " visEnFracX = " << visEnFracX << std::endl;
    std::cout << " angle(vis, tau) = " << TMath::ACos(compScalarProduct(normalize(hypothesis_T->p4().Vect()), normalize(hypothesis_T->p4_fitted().Vect()))) << std::endl;
  }
#endif
  // LB: normalize likelihood function such that 
  //               1
  //       integral  prob dX dMnunu = 1.
  //               0
  double prob = 1.;
  if ( polSign == +1 || polSign == -1 ) {

    if ( nuMass < TMath::Sqrt((1. - visEnFracX)*tauLeptonMass2) ) {
      prob = (nuMass/(4.*tauLeptonMass4))*((tauLeptonMass2 + 2.*nuMass2)*(tauLeptonMass2 - nuMass2) 
				         + polSign*(tauLeptonMass2*(2.*visEnFracX - 1.) + nuMass2)*(-tauLeptonMass2 + 2.*nuMass2));
    } else {
      double nuMass_limit  = TMath::Sqrt((1. - visEnFracX)*tauLeptonMass2);
      double nuMass2_limit = square(nuMass_limit);
      prob = (nuMass_limit/(4.*tauLeptonMass4))*((tauLeptonMass2 + 2.*nuMass2_limit)*(tauLeptonMass2 - nuMass2_limit) 
	    + polSign*(tauLeptonMass2*(2.*visEnFracX - 1.) + nuMass2_limit)*(-tauLeptonMass2 + 2.*nuMass2_limit));
      prob /= (1. + 1.e+6*square(nuMass - nuMass_limit));
    }

    if ( applyVisPtCutCorrection_ ) {
      double probCorr = 1.;
      if ( hypothesis_T->p4_fitted().pt() > visPtCutThreshold_ ) {
	double xCut = visPtCutThreshold_/hypothesis_T->p4_fitted().pt();
	probCorr = 1./((0.5*(1. + polSign)*(1./3.)*(0.5 - xCut + cube(xCut) - 0.5*fourth(xCut))
		      + 0.5*(1. - polSign)*(1./3.)*(0.75 - xCut + 0.25*fourth(xCut))));
      }
      //if ( this->verbosity_ ) std::cout << "probCorr (lep) = " << probCorr << std::endl;
      prob *= probCorr;
    }
  } else {
    if ( nuMass < TMath::Sqrt((1. - visEnFracX)*tauLeptonMass2) ) { // LB: physical solution
      prob = (13./tauLeptonMass4)*(tauLeptonMass2 - nuMass2)*(tauLeptonMass2 + 2.*nuMass2)*nuMass;
    } else {                                                        // LB: unphysical solution
      double nuMass_limit  = TMath::Sqrt((1. - visEnFracX)*tauLeptonMass2);
      double nuMass2_limit = square(nuMass_limit);
      prob = (13./tauLeptonMass4)*(tauLeptonMass2 - nuMass2_limit)*(tauLeptonMass2 + 2.*nuMass2_limit)*nuMass_limit;
      prob /= (1. + 1.e+6*square(nuMass - nuMass_limit));
    }
    if ( applyVisPtCutCorrection_ ) {
      double probCorr = 1.;
      const double epsilon_regularization = 1.e-3;
      if ( hypothesis_T->p4_fitted().pt() > visPtCutThreshold_ ) {
	double xCut = visPtCutThreshold_/hypothesis_T->p4_fitted().pt();
	probCorr = 1./((3. - 5.*xCut + 3.*cube(xCut) - fourth(xCut)) + epsilon_regularization);
      } 
      //if ( this->verbosity_ ) std::cout << "probCorr (lep) = " << probCorr << std::endl;
      prob *= probCorr;
    }
  }
  if ( applySinThetaFactor_ ) prob *= (0.5*TMath::Sin(decayAngle));

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

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

typedef SVfitTauToLepLikelihoodMatrixElement<pat::Electron> SVfitTauToElecLikelihoodMatrixElement;
typedef SVfitTauToLepLikelihoodMatrixElement<pat::Muon> SVfitTauToMuLikelihoodMatrixElement;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(SVfitSingleParticleLikelihoodPluginFactory, SVfitTauToElecLikelihoodMatrixElement, "SVfitTauToElecLikelihoodMatrixElement");
DEFINE_EDM_PLUGIN(SVfitSingleParticleLikelihoodPluginFactory, SVfitTauToMuLikelihoodMatrixElement, "SVfitTauToMuLikelihoodMatrixElement");

