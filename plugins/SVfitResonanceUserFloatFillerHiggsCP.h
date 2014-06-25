#ifndef TauAnalysis_SVfit_SVfitResonanceUserFloatFillerHiggsCP_h
#define TauAnalysis_SVfit_SVfitResonanceUserFloatFillerHiggsCP_h

/** \class SVfitResonanceBuilder
 *
 * Auxiliary class for computing Higgs CP-sensitive observable, phiStar and PsiStar
 * as described in the paper
 *   "Determining the CP parity of Higgs bosons at the LHC in the tau to 1-prong decay channels",
 *   S. Berge and W. Bernreuther,
 *   arXiv:0812.1910v1 
 *
 * \author Christian Veelken, LLR/Ecole Polytechnique
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "TauAnalysis/SVfit/interface/SVfitResonanceUserFloatFillerBase.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitResonanceHypothesis.h"

#include <TH1.h>

class SVfitResonanceUserFloatFillerHiggsCP : public SVfitResonanceUserFloatFillerBase
{
 public:
  SVfitResonanceUserFloatFillerHiggsCP(const edm::ParameterSet&); 
  ~SVfitResonanceUserFloatFillerHiggsCP();

  void beginCandidate(const SVfitResonanceHypothesis*);

  void resetHistograms();

  void fillHistograms(const SVfitResonanceHypothesis*);

  void addUserFloats(SVfitResonanceHypothesis*) const;

 private:
  TH1* histogramPhiStar_;
  TH1* histogramPhiStarWeighted_;
  TH1* histogramPsiStar_;
  TH1* histogramPsiStarWeighted_;

  int leg1DecayMode_;
  reco::Candidate::LorentzVector leg1DistPionP4_;
  int leg2DecayMode_;
  reco::Candidate::LorentzVector leg2DistPionP4_;
};

#endif


