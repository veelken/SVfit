#include "TauAnalysis/SVfit/plugins/SVfitResonanceUserFloatFillerHiggsCP.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/TauReco/interface/PFTauDecayMode.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitTauDecayHypothesis.h"

#include "TauAnalysis/SVfit/interface/svFitAuxFunctionsHiggsCP.h"
#include "TauAnalysis/SVfit/interface/candidateAuxFunctions.h"

#include <TMath.h>

#include <string>

using namespace svFit_namespace;

SVfitResonanceUserFloatFillerHiggsCP::SVfitResonanceUserFloatFillerHiggsCP(const edm::ParameterSet& cfg)
  : SVfitResonanceUserFloatFillerBase(cfg)    
{
  std::string histogramNamePhiStar = std::string(pluginName_).append("_histogramPhiStar");
  histogramPhiStar_ = new TH1D(histogramNamePhiStar.data(), histogramNamePhiStar.data(), 180, 0., TMath::Pi());
  std::string histogramNamePhiStarWeighted = std::string(pluginName_).append("_histogramPhiStarWeighted");
  histogramPhiStarWeighted_ = new TH1D(histogramNamePhiStarWeighted.data(), histogramNamePhiStarWeighted.data(), 180, 0., TMath::Pi());
  std::string histogramNamePsiStar = std::string(pluginName_).append("_histogramPsiStar");
  histogramPsiStar_ = new TH1D(histogramNamePsiStar.data(), histogramNamePsiStar.data(), 180, 0., TMath::Pi());
  std::string histogramNamePsiStarWeighted = std::string(pluginName_).append("_histogramPsiStarWeighted");
  histogramPsiStarWeighted_ = new TH1D(histogramNamePsiStarWeighted.data(), histogramNamePsiStarWeighted.data(), 180, 0., TMath::Pi());
}

SVfitResonanceUserFloatFillerHiggsCP::~SVfitResonanceUserFloatFillerHiggsCP()
{
  delete histogramPhiStar_;
  delete histogramPsiStar_;
}

namespace
{
  void setTauDecayMode_and_distPionP4(const SVfitSingleParticleHypothesis* daughter, int& decayMode, reco::Candidate::LorentzVector& distPionP4)
  {
    const reco::Candidate* particle = daughter->particle().get();
    assert(particle);
    if ( dynamic_cast<const pat::Tau*>(particle) ) {
      const pat::Tau* tau = dynamic_cast<const pat::Tau*>(particle);
      decayMode = tau->decayMode();
      const reco::Candidate* distPion = getDistPion(*tau);
      if ( distPion ) distPionP4 = distPion->p4();
      else distPionP4.SetXYZT(0.,0.,0.,0.);
    } else if ( dynamic_cast<const pat::Electron*>(particle) ) {
      decayMode = reco::PFTauDecayMode::tauDecaysElectron;
      distPionP4.SetXYZT(0.,0.,0.,0.);
    } else if ( dynamic_cast<const pat::Muon*>(particle) ) {
      decayMode = reco::PFTauDecayMode::tauDecayMuon;
      distPionP4.SetXYZT(0.,0.,0.,0.);
    } else assert(0);
  }

  double getB(const SVfitTauDecayHypothesis* daughter, int decayMode, const reco::Candidate::LorentzVector& distPionP4)
  {
    double b = 1.;
    if ( decayMode == reco::PFTauDecayMode::tauDecaysElectron || decayMode == reco::PFTauDecayMode::tauDecayMuon ) {
      b = b_lep(daughter);
    } else if ( decayMode == reco::PFTau::kOneProng0PiZero ) {
      b = b_pi();
    } else if ( decayMode == reco::PFTau::kOneProng1PiZero || reco::PFTau::kOneProng2PiZero ) {
      b = 0.728*b_rho(daughter, distPionP4) + 0.272*b_a1(daughter, distPionP4); // CV: return average of rho -> pi- pi0 and a1 -> pi- pi0 pi0 decays, weighted by branching fraction
    } else if ( decayMode == reco::PFTau::kThreeProng0PiZero ) {
      b = b_a1(daughter, distPionP4);
    }
    return b;
  }
}

void SVfitResonanceUserFloatFillerHiggsCP::beginCandidate(const SVfitResonanceHypothesis* resonance)
{
  assert(resonance->numDaughters() == 2);
  const SVfitSingleParticleHypothesis* daughter1 = resonance->daughter(0);
  assert(daughter1);
  setTauDecayMode_and_distPionP4(daughter1, leg1DecayMode_, leg1DistPionP4_);
  const SVfitSingleParticleHypothesis* daughter2 = resonance->daughter(1);
  assert(daughter2);
  setTauDecayMode_and_distPionP4(daughter2, leg2DecayMode_, leg2DistPionP4_);
}

void SVfitResonanceUserFloatFillerHiggsCP::resetHistograms()
{
  histogramPhiStar_->Reset();
  histogramPhiStarWeighted_->Reset();
  histogramPsiStar_->Reset();
  histogramPsiStarWeighted_->Reset();
}

void SVfitResonanceUserFloatFillerHiggsCP::fillHistograms(const SVfitResonanceHypothesis* resonance)
{
  assert(resonance->numDaughters() == 2);
  const SVfitSingleParticleHypothesis* daughter1 = resonance->daughter(0);
  const SVfitTauDecayHypothesis* daughter1_nonbase = dynamic_cast<const SVfitTauDecayHypothesis*>(daughter1);
  assert(daughter1_nonbase);
  const SVfitSingleParticleHypothesis* daughter2 = resonance->daughter(1);
  const SVfitTauDecayHypothesis* daughter2_nonbase = dynamic_cast<const SVfitTauDecayHypothesis*>(daughter2);
  assert(daughter2_nonbase);
  double phiStar = -1.;
  double PsiStar = -1.;
  bool error = false;
  compPhiStar_and_PsiStar(daughter1_nonbase->p4_fitted(), daughter1_nonbase->p4(), daughter1_nonbase->particle()->charge(),
			  daughter2_nonbase->p4_fitted(), daughter2_nonbase->p4(), daughter2_nonbase->particle()->charge(),
			  resonance->p4_fitted(),
			  phiStar, PsiStar, error);
  if ( !error ) {
    double leg1_b = getB(daughter1_nonbase, leg1DecayMode_, leg1DistPionP4_);
    double leg2_b = getB(daughter2_nonbase, leg2DecayMode_, leg2DistPionP4_);
    double weight = leg1_b*leg2_b;
    histogramPhiStar_->Fill(phiStar);
    histogramPhiStarWeighted_->Fill(phiStar*weight);
    histogramPsiStar_->Fill(PsiStar);
    histogramPsiStarWeighted_->Fill(PsiStar*weight);
  }
}

void SVfitResonanceUserFloatFillerHiggsCP::addUserFloats(SVfitResonanceHypothesis* resonance) const
{
  addUserFloat(resonance, "phiStar",         histogramPhiStar_);
  addUserFloat(resonance, "phiStarWeighted", histogramPhiStarWeighted_);
  addUserFloat(resonance, "PsiStar",         histogramPsiStar_);
  addUserFloat(resonance, "PsiStarWeighted", histogramPsiStarWeighted_);
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(SVfitResonanceUserFloatFillerPluginFactory, SVfitResonanceUserFloatFillerHiggsCP, "SVfitResonanceUserFloatFillerHiggsCP");
