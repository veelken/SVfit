#ifndef TauAnalysis_SVfit_SVfitTauToHadLikelihoodPhaseSpace_h
#define TauAnalysis_SVfit_SVfitTauToHadLikelihoodPhaseSpace_h

/** \class SVfitTauToLeptonLikelihoodPhaseSpace
 *
 * Plugin to compute likelihood for system of hadrons to be compatible 
 * with tau --> tau-jet + nu two-body decay,
 * assuming constant matrix element, so that energy and angular distribution 
 * of decay products are solely determined by phase-space
 * 
 * \author Evan Friis, Christian Veelken; UC Davis
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/SVfit/interface/SVfitSingleParticleLikelihood.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitSingleParticleHypothesisBase.h"

#include <TFile.h>
#include <TH1.h>
#include <TMath.h>

namespace svFitTauToHadLikelihoodPhaseSpace_namespace
{
  struct lutEntryType
  {
    lutEntryType(const std::string& pluginName, const edm::ParameterSet& cfg)
      : histogram_(0),
	firstBin_(-1),
	lastBin_(-1)
    {
      if ( cfg.exists("tauDecayModes") ) {
	tauDecayModes_ = cfg.getParameter<vint>("tauDecayModes");
      } else {
	tauDecayModes_.push_back(-1); // CV: all tau decay modes
      }
      edm::FileInPath inputFileName = cfg.getParameter<edm::FileInPath>("inputFileName");
      if ( !inputFileName.isLocal() ) 
	throw cms::Exception("SVfitTauToHadLikelihoodPhaseSpace") 
	  << " Failed to find File = " << inputFileName << " !!\n";
      std::string histogramName = cfg.getParameter<std::string>("histogramName");
      inputFileName_full_ = inputFileName.fullPath().data();
      TFile* inputFile = new TFile(inputFileName_full_.data());
      //TFile* inputFile = 0;
      //if ( inputFiles_.find(inputFileName_full_) != inputFiles_.end() ) {
      //  inputFile = inputFiles_[inputFileName_full_];
      //  inputFileUsage_[inputFileName_full_] = inputFileUsage_[inputFileName_full_] + 1;
      //} else {
      //  inputFile = new TFile(inputFileName_full_.data());
      //  inputFiles_[inputFileName_full_] = inputFile;
      //  inputFileUsage_[inputFileName_full_] = 1;
      //}
      TH1* histogramTmp = dynamic_cast<TH1*>(inputFile->Get(histogramName.data()));
      if ( !histogramTmp )
	throw cms::Exception("SVfitTauToHadLikelihoodPhaseSpace") 
	  << " Failed to load histogram = " << histogramName.data() << " from file = " << inputFileName.fullPath().data() << " !!\n";
      histogram_ = (TH1*)histogramTmp->Clone(std::string(pluginName).append("_").append(histogramTmp->GetName()).data());
      //histogram_ = histogramTmp;
      firstBin_ = 1;
      lastBin_ = histogram_->GetNbinsX();
      delete inputFile;
    }
    ~lutEntryType()
    {
      delete histogram_;
      //inputFileUsage_[inputFileName_full_] = inputFileUsage_[inputFileName_full_] - 1;
      //if ( inputFileUsage_[inputFileName_full_] == 0 ) {
      //  delete inputFiles_[inputFileName_full_];
      //  inputFiles_.erase(inputFileName_full_);
      //}
    }
    typedef std::vector<int> vint;
    vint tauDecayModes_;
    TH1* histogram_;
    int firstBin_;
    int lastBin_;
    std::string inputFileName_full_;
    //static std::map<std::string, TFile*> inputFiles_;  // key = inputFileName
    //static std::map<std::string, int> inputFileUsage_; // key = inputFileName
  };
}

class SVfitTauToHadLikelihoodPhaseSpace : public SVfitSingleParticleLikelihood
{
 public:
  SVfitTauToHadLikelihoodPhaseSpace(const edm::ParameterSet&);
  ~SVfitTauToHadLikelihoodPhaseSpace();

  void beginJob(SVfitAlgorithmBase*);

  void beginCandidate(const SVfitSingleParticleHypothesis*);

  double operator()(const SVfitSingleParticleHypothesis*, int) const;

 private:
  bool applySinThetaFactor_; 

  bool varyVisMass_;
  svFitTauToHadLikelihoodPhaseSpace_namespace::lutEntryType* lutVisMass_;

  bool shiftVisMass_;
  std::vector<svFitTauToHadLikelihoodPhaseSpace_namespace::lutEntryType*> lutEntriesVisMass_shift_;
  const svFitTauToHadLikelihoodPhaseSpace_namespace::lutEntryType* lutVisMass_shift_;
  double visMass_unshifted_;
  bool shiftVisPt_;
  std::vector<svFitTauToHadLikelihoodPhaseSpace_namespace::lutEntryType*> lutEntriesVisPt_shift_;
  const svFitTauToHadLikelihoodPhaseSpace_namespace::lutEntryType* lutVisPt_shift_;
  
  SVfitAlgorithmBase* algorithm_;
};

#endif
