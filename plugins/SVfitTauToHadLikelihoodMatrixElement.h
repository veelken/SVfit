#ifndef TauAnalysis_SVfit_SVfitTauToHadLikelihoodMatrixElement_h
#define TauAnalysis_SVfit_SVfitTauToHadLikelihoodMatrixElement_h

/** \class SVfitTauToHadLikelihoodMatrixElement
 *
 * Plugin for computing likelihood for tau lepton decay 
 *  tau- --> X nu
 * to be compatible with matrix element of V-A electroweak decay
 * (taking tau lepton polarization effects into account)
 *
 * NOTE: the system of hadrons X can either be pi-, rho- --> pi- pi0,
 *       a1- --> pi- pi0 pi0 or a1- --> pi- pi+ pi-;
 *       tau decays into pi- pi+ pi- pi0 are **not** supported
 *
 * \author Christian Veelken, Lorenzo Bianchini; LLR
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/SVfit/interface/SVfitSingleParticleLikelihood.h"

#include "AnalysisDataFormats/TauAnalysis/interface/SVfitSingleParticleHypothesisBase.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include <TFile.h>
#include <TGraph.h>
#include <TMatrixD.h>
#include <TVectorD.h>

class SVfitTauToHadLikelihoodMatrixElement : public SVfitSingleParticleLikelihood
{
 public:
  SVfitTauToHadLikelihoodMatrixElement(const edm::ParameterSet&);
  ~SVfitTauToHadLikelihoodMatrixElement();

  void beginJob(SVfitAlgorithmBase*);

  void beginCandidate(const SVfitSingleParticleHypothesis*);

  double operator()(const SVfitSingleParticleHypothesis*, int) const;

 private:
  double compProb_pionDecay(double, int, double) const;
  double compProb_rhoDecay(double, double, int, double) const;
  double compProb_a1Decay(double, double, int, double) const;

  bool applySinThetaFactor_; 

  edm::FileInPath inputFileName_;
  TFile* inputFile_;

  TGraph* rhoLPlus_;
  TGraph* rhoNormLPlus_; 
  TGraph* rhoLMinus_;
  TGraph* rhoNormLMinus_; 
  TGraph* rhoTPlus_;
  TGraph* rhoNormTPlus_;
  TGraph* rhoTMinus_;
  TGraph* rhoNormTMinus_;

  TGraph* a1LPlus_;
  TGraph* a1NormLPlus_; 
  TGraph* a1LMinus_;
  TGraph* a1NormLMinus_; 
  TGraph* a1TPlus_;
  TGraph* a1NormTPlus_;
  TGraph* a1TMinus_;
  TGraph* a1NormTMinus_;

  TGraph* a1Lz_;
  TGraph* a1Tz_;

  TGraph* rhoyMinus_;
  TGraph* rhoyNormMinus_;
  TGraph* rhoyPlus_;
  TGraph* rhoyNormPlus_;
  TGraph* a1yMinus_;
  TGraph* a1yNormMinus_;
  TGraph* a1yPlus_;
  TGraph* a1yNormPlus_;

  std::vector<int> supportedTauDecayModes_;
  unsigned numSupportedTauDecayModes_;

  mutable TVectorD vRec_;
  TMatrixD recToGenTauDecayModeMap_;
  mutable TVectorD vGen_;
  mutable TVectorD vProb_;

  mutable long numWarningsUnphysicalDecay_rho_;
  mutable long numWarningsUnphysicalDecay_a1_;

  SVfitAlgorithmBase* algorithm_;
};

#endif
