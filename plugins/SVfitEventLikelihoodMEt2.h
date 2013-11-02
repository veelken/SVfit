#ifndef TauAnalysis_SVfit_SVfitEventLikelihoodMEt2_h
#define TauAnalysis_SVfit_SVfitEventLikelihoodMEt2_h

/** \class SVfitEventLikelihoodMEt2
 *
 * Plugin for computing likelihood for neutrinos produced in tau lepton decays
 * to match missing transverse momentum reconstructed in the event
 *
 * New version using covariance matrix of (PF)MET significance calculation
 * (CMS AN-10/400) to compute the likehood
 *
 * \author Christian Veelken, UC Davis
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "RecoMET/METAlgorithms/interface/SignAlgoResolutions.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "TauAnalysis/SVfit/interface/SVfitEventLikelihood.h"
#include "TauAnalysis/SVfit/interface/PFMEtSignInterface.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitEventHypothesis.h"

#include <TMatrixD.h>
#include <TVectorD.h>
#include <TFormula.h>

#include <list>

class SVfitEventLikelihoodMEt2 : public SVfitEventLikelihood
{
 public:
  SVfitEventLikelihoodMEt2(const edm::ParameterSet&);
  ~SVfitEventLikelihoodMEt2();

  void beginJob(SVfitAlgorithmBase*);
  void beginEvent(const edm::Event&, const edm::EventSetup&);
  void beginCandidate(const SVfitEventHypothesis*) const;

  double operator()(const SVfitEventHypothesis*) const;

 private:

  double power_;

  edm::InputTag srcMEtCovMatrix_;
  double sfMEtCov_;

  PFMEtSignInterface* pfMEtSign_;

  mutable TMatrixD pfMEtCov_;
  mutable double   pfMEtCovDet_;
  mutable TMatrixD pfMEtCovInverse_;

  mutable double pfMEtCovInverse00_;
  mutable double pfMEtCovInverse01_;
  mutable double pfMEtCovInverse10_;
  mutable double pfMEtCovInverse11_;

  mutable double residual_fitted0_;
  mutable double residual_fitted1_;

  mutable double nllConstTerm_;

  struct tailProbCorrFunctionType
  {
    tailProbCorrFunctionType(const std::string& pluginName, const edm::ParameterSet& cfg)
      : function_(0)
    {
      std::string formula = cfg.getParameter<std::string>("formula");
      std::string functionName = Form("%s_formula", pluginName.data());
      function_ = new TFormula(functionName.data(), formula.data());
      xMin_ = cfg.getParameter<double>("xMin");
      xMax_ = cfg.getParameter<double>("xMax");
      numParameter_ = function_->GetNpar();   
      if ( numParameter_ > 0 ) {
	edm::ParameterSet cfgPars = cfg.getParameter<edm::ParameterSet>("parameter");
	parameter_.resize(numParameter_);
	for ( int iParameter = 0; iParameter < numParameter_; ++iParameter ) {
	  parameter_[iParameter] = cfgPars.getParameter<double>(Form("par%i", iParameter));
	  function_->SetParameter(iParameter, parameter_[iParameter]);
	}
      }
    }

    ~tailProbCorrFunctionType()
    {
      delete function_;
    }

    double eval(double x) const
    {
      double x_limited = x;
      if ( x_limited < xMin_ ) x_limited = xMin_;
      if ( x_limited > xMax_ ) x_limited = xMax_;
      return function_->Eval(x_limited);
    }
    
    TFormula* function_;
    double xMin_;
    double xMax_;
    int numParameter_;
    std::vector<double> parameter_;
  };

  tailProbCorrFunctionType* tailProbCorrFunction_;

  bool monitorMEtUncertainty_;
  std::string monitorFilePath_;
  std::string monitorFileName_;
  unsigned numToys_;
};

#endif
