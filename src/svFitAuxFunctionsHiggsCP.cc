#include "TauAnalysis/SVfit/interface/svFitAuxFunctionsHiggsCP.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitTauDecayHypothesis.h"

#include "TauAnalysis/SVfit/interface/svFitAuxFunctions.h"

#include <TMath.h>
#include <Math/VectorUtil.h>
#include <TRandom3.h>

namespace svFit_namespace
{
  AlgebraicVector3 convertToAlgebraicVector3(const reco::Candidate::LorentzVector& p4)
  {
    return AlgebraicVector3(p4.px(), p4.py(), p4.pz());
  }

  AlgebraicVector3 compNstar(const reco::Candidate::LorentzVector& tauP4, const reco::Candidate::LorentzVector& higgsP4)
  {
    reco::Candidate::LorentzVector tauFlightPath_lab(tauP4.px(), tauP4.py(), tauP4.pz(), 0.); 
    reco::Candidate::LorentzVector tauFlightPath_rf = boostToCOM(higgsP4, tauFlightPath_lab);
    AlgebraicVector3 nStar(tauFlightPath_rf.px(), tauFlightPath_rf.py(), tauFlightPath_rf.pz());
    return nStar;
  }
  AlgebraicVector3 compPhatStar(const reco::Candidate::LorentzVector& visTauP4, const reco::Candidate::LorentzVector& higgsP4)
  {
    reco::Candidate::LorentzVector visP4_rf = boostToCOM(higgsP4, visTauP4);
    AlgebraicVector3 pHatStar = normalize(convertToAlgebraicVector3(visP4_rf));
    return pHatStar;
  }

  AlgebraicVector3 compPerpComponent(const AlgebraicVector3& vec, const AlgebraicVector3& dir)
  {
    AlgebraicVector3 dirHat = normalize(dir);
    double r_parl = vec(0)*dirHat(0) + vec(1)*dirHat(1) + vec(2)*dirHat(2);    
    AlgebraicVector3 vec_perp;
    vec_perp(0) = vec(0) - r_parl*dirHat(0);
    vec_perp(1) = vec(1) - r_parl*dirHat(1);
    vec_perp(2) = vec(2) - r_parl*dirHat(2);
    return vec_perp;
  }

  void compPhiStar_and_PsiStar(const reco::Candidate::LorentzVector& tau1P4, const reco::Candidate::LorentzVector& tau1VisP4, double tau1Charge,
			       const reco::Candidate::LorentzVector& tau2P4, const reco::Candidate::LorentzVector& tau2VisP4, double tau2Charge,
			       const reco::Candidate::LorentzVector& higgsP4, 
			       double& phiStar, double& PsiStar, bool& error)
  {
    //
    // compute Higgs-CP sensitive observables "phiStar" and "PsiStar", 
    // as described in the paper
    //   "Determining the CP parity of Higgs bosons at the LHC in the tau to 1-prong decay channels",
    //   S. Berge and W. Bernreuther,
    //   arXiv:0812.1910v1 
    //
    AlgebraicVector3 nStar1 = compNstar(tau1P4, higgsP4);
    AlgebraicVector3 pHatStar1 = compPhatStar(tau1VisP4, tau1P4);
    AlgebraicVector3 nStar1_perp = compPerpComponent(nStar1, pHatStar1);
    AlgebraicVector3 nHatStar1_perp = normalize(nStar1_perp);

    AlgebraicVector3 nStar2 = compNstar(tau2P4, higgsP4);
    AlgebraicVector3 pHatStar2 = compPhatStar(tau2VisP4, tau2P4);
    AlgebraicVector3 nStar2_perp = compPerpComponent(nStar2, pHatStar2);
    AlgebraicVector3 nHatStar2_perp = normalize(nStar2_perp);

    double cosPhiStar = compScalarProduct(nHatStar1_perp, nHatStar2_perp);
    assert(cosPhiStar >= -1. && cosPhiStar <= +1.);
    phiStar = TMath::ACos(cosPhiStar);

    bool mapPlusTo1;
    if ( tau1Charge > 0. && tau2Charge < 0. ) {
      mapPlusTo1 = true;
    } else if ( tau1Charge < 0. && tau2Charge > 0. ) {
      mapPlusTo1 = false;
    } else {
      static TRandom3 rnd;
      if ( rnd.Rndm() > 0.5 ) mapPlusTo1 = true;
      else mapPlusTo1 = false;
    }
    AlgebraicVector3 nHatStarPlus_perp, pHatStarPlus, nHatStarMinus_perp, pHatStarMinus;
    if ( mapPlusTo1 ) {
      nHatStarPlus_perp = nHatStar1_perp;
      pHatStarPlus = pHatStar1;
      nHatStarMinus_perp = nHatStar2_perp;
      pHatStarMinus = pHatStar2;
    } else {
      nHatStarPlus_perp = nHatStar2_perp;
      pHatStarPlus = pHatStar2;
      nHatStarMinus_perp = nHatStar1_perp;
      pHatStarMinus = pHatStar1;
    } 
    
    AlgebraicVector3 nHatStarPlus_perp_cross_nHatStarMinus_perp = compCrossProduct(nHatStarPlus_perp, nHatStarMinus_perp);
    double cosPsiStar = compScalarProduct(pHatStarMinus, nHatStarPlus_perp_cross_nHatStarMinus_perp);
    assert(cosPsiStar >= -1. && cosPsiStar <= +1.);
    PsiStar = TMath::ACos(cosPsiStar);
  }

  // define spectral functions for spin analyzer power of different tau decay modes,
  // taken from the appendix of arXiv:0812.1910v1 
  double b_pi()
  {
    return 1.;
  }
  double b_rho(const SVfitTauDecayHypothesis* daughter, const reco::Candidate::LorentzVector& chargedPionP4_lab)
  {
    const double p = 4.*chargedPionMass2/tauLeptonMass2;
    const double r = rhoMesonMass2/tauLeptonMass2;
    reco::Candidate::Vector boostToTauRF_vector = daughter->p4_fitted().BoostToCM();
    double enPi_rf = ROOT::Math::VectorUtil::boost(chargedPionP4_lab, boostToTauRF_vector).E();
    double x = 4.*enPi_rf/tauLeptonMass;
    const double enPi_rfMin = 0.25*tauLeptonMass*TMath::Sqrt(1. + r - (1. - r)*TMath::Sqrt(1. - (p/r)));
    const double enPi_rfMax = 0.25*tauLeptonMass*TMath::Sqrt(1. + r + (1. - r)*TMath::Sqrt(1. - (p/r)));
    double b;
    if ( enPi_rf > enPi_rfMin && enPi_rf < enPi_rfMax ) {
      b = (x*square(x - r - 1.) + x*(3. - r)*(r - p) - 4.*(r - p))/TMath::Sqrt(x*x - 4.*p*(square(x - r - 1.) + (1. - r)*(r - p)));
    } else {
      b = 0.;
    }
    return b;
  }
  double b_a1(const SVfitTauDecayHypothesis* daughter, const reco::Candidate::LorentzVector& chargedPionP4_lab)
  {
    reco::Candidate::Vector boostToTauRF_vector = daughter->p4_fitted().BoostToCM();
    double enPi_rf = ROOT::Math::VectorUtil::boost(chargedPionP4_lab, boostToTauRF_vector).E();
    double x = (2.*tauLeptonMass*(enPi_rf - chargedPionMass))/(tauLeptonMass2 - 3.*chargedPionMass2 - 2.*tauLeptonMass*chargedPionMass);
    const double enPi_rfMin = chargedPionMass;
    const double enPi_rfMax = (tauLeptonMass2 - 3.*chargedPionMass2)/(2.*tauLeptonMass);
    double b;
    if ( enPi_rf > enPi_rfMin && enPi_rf < enPi_rfMax ) {
      if ( x < 0. ) x = 0.;
      //b = -5.28726*TMath::Sqrt(x) + x*(9.38612 + x*(-1.26356 + x*(-18.9094 + x*(36.0517 - x*19.4113))));
      b = -5.28726*TMath::Sqrt(x) + 9.38612*x - 1.26356*square(x) - 18.9094*cube(x) + 36.0517*fourth(x) - 19.4113*fifth(x);
    } else {
      b = 0.;
    }
    return b;
  }
  double b_lep(const SVfitTauDecayHypothesis* daughter)
  {
    double enLep_rf = daughter->p4vis_rf().E();
    double x = 2.*enLep_rf/tauLeptonMass;
    const double enLep_rfMin = 0.;
    const double enLep_rfMax = 0.5*tauLeptonMass;
    double b;
    if ( enLep_rf > enLep_rfMin && enLep_rf < enLep_rfMax ) {
      b = (1. - 2*x)/(3. - 2*x);
    } else {
      b = 0.;
    }
    return b;
  }
}
