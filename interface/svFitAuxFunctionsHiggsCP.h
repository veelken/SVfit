#ifndef TauAnalysis_SVfit_svFitAuxFunctionsHiggsCP_h
#define TauAnalysis_SVfit_svFitAuxFunctionsHiggsCP_h

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/AlgebraicROOTObjects.h"

class SVfitTauDecayHypothesis;

namespace svFit_namespace
{
  AlgebraicVector3 convertToAlgebraicVector3(const reco::Candidate::LorentzVector&);

  AlgebraicVector3 compNstar(const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&);
  AlgebraicVector3 compPhatStar(const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&);

  AlgebraicVector3 compPerpComponent(const AlgebraicVector3&, const AlgebraicVector3&);
  
  void compPhiStar_and_PsiStar(const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&, double,
			       const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&, double,
			       const reco::Candidate::LorentzVector&, 
			       double&, double&, bool&);
  
  double b_pi();
  double b_rho(const SVfitTauDecayHypothesis*, const reco::Candidate::LorentzVector&);
  double b_a1(const SVfitTauDecayHypothesis*, const reco::Candidate::LorentzVector&);
  double b_lep(const SVfitTauDecayHypothesis*);
}

#endif
