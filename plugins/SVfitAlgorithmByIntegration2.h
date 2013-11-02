#ifndef TauAnalysis_SVfit_SVfitAlgorithmByIntegration2_h
#define TauAnalysis_SVfit_SVfitAlgorithmByIntegration2_h

/** \class SVfitAlgorithmByIntegration2
 *
 * Concrete implementation of (n)SVfit algorithm
 * by integration of likelihood functions
 * (based on Markov Chain integration algorithm)
 *
 * \author Christian Veelken, UC Davis
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/Ptr.h"

#include "CommonTools/Utils/interface/StringObjectFunction.h"

#include "TauAnalysis/SVfit/interface/SVfitAlgorithmBase.h"
#include "TauAnalysis/SVfit/interface/MarkovChainIntegrator.h"
#include "TauAnalysis/SVfit/interface/svFitAuxFunctions.h"

#include <Math/Functor.h>
#include <TH1.h>

#include <vector>
#include <string>

class SVfitAlgorithmByIntegration2 : public SVfitAlgorithmBase
{
 public:
  SVfitAlgorithmByIntegration2(const edm::ParameterSet&);
  ~SVfitAlgorithmByIntegration2();

  void beginJob();
  void beginEvent(const edm::Event&, const edm::EventSetup&);

  unsigned getNumDimensions() const { return numDimensions_; }

  bool applyJacobiFactors() const 
  {
    return applyJacobiFactors_; 
  }

  void print(std::ostream&) const {}

  bool update(const double* x, const double* param) const;
  double nll(const double* x, const double* param) const;

  void fillProbHistograms(const double*);

 protected:
  TH1* bookPtHistogram(const std::string& name);
  TH1* bookEtaHistogram(const std::string& name);
  TH1* bookPhiHistogram(const std::string& name);
  TH1* bookMassHistogram(const std::string& name);

  void fitImp() const;

  void setMassResults(SVfitResonanceHypothesisBase*, const TH1*) const;
    
  bool isDaughter(const std::string&);
  bool isResonance(const std::string&);

  struct replaceParBase
  {
    virtual void beginJob(SVfitAlgorithmByIntegration2*) {}
    virtual double operator()(const double*) const = 0;
    int iPar_;
  };

  TFormula* makeReplacementFormula(const std::string&, const std::string&, std::vector<replaceParBase*>&, int&);
  
  SVfitParameter* getFitParameter(const std::string&);
  
  struct SVfitParameterMappingType
  {
    SVfitParameterMappingType(const SVfitParameter* base)
      : base_(base)
    {}
    const SVfitParameter* base_;
    int idxByIntegration_;
  };

  std::vector<SVfitParameterMappingType> fitParameterMappings_;

  edm::RunNumber_t currentRunNumber_;
  edm::LuminosityBlockNumber_t currentLumiSectionNumber_;
  edm::EventNumber_t currentEventNumber_;

  struct replaceParByFitParameter : replaceParBase
  {
    void beginJob(SVfitAlgorithmByIntegration2* algorithm)
    {
      idx_ = algorithm->getFitParameter(fitParameterName_)->index();
    }
    double operator()(const double* param) const { return param[idx_]; }
    std::string fitParameterName_;
    int idx_;
  };

  struct replaceParByResonanceHypothesis : replaceParBase
  {
    replaceParByResonanceHypothesis()
      : valueExtractor_(0)
    {}
    ~replaceParByResonanceHypothesis()
    {
      delete valueExtractor_;
    }
    double operator()(const double* param) const { return value_; }
    std::string resonanceName_;
    SVfitResonanceHypothesis* resonanceHypothesis_;
    StringObjectFunction<SVfitResonanceHypothesis>* valueExtractor_;
    mutable double value_;
  };

  struct fitParameterReplacementType
  {    
    fitParameterReplacementType()
      : replaceBy_(0),
	deltaFuncDerrivative_(0)
    {}
    ~fitParameterReplacementType() 
    {
      delete replaceBy_;
      delete deltaFuncDerrivative_;
      for ( std::vector<replaceParBase*>::iterator it = parForReplacements_.begin();
	    it != parForReplacements_.end(); ++it ) {
	delete (*it);
      }
      for ( std::vector<replaceParBase*>::iterator it = parForDeltaFuncDerrivative_.begin();
	    it != parForDeltaFuncDerrivative_.end(); ++it ) {
	delete (*it);
      }
    }
    void beginJob(SVfitAlgorithmByIntegration2* algorithm)
    {
      SVfitParameter* fitParameterToReplace = algorithm->getFitParameter(toReplace_);
      if ( !fitParameterToReplace ) {
	throw cms::Exception("fitParameterReplacementType::beginJob")
	  << " No fitParameter of name = " << toReplace_ << " defined !!";
      }
      idxToReplace_ = fitParameterToReplace->index();

      for ( std::vector<replaceParBase*>::iterator par = parForReplacements_.begin();
	    par != parForReplacements_.end(); ++par ) {
	(*par)->beginJob(algorithm);
      }
      for ( std::vector<replaceParBase*>::iterator par = parForDeltaFuncDerrivative_.begin();
	    par != parForDeltaFuncDerrivative_.end(); ++par ) {
	(*par)->beginJob(algorithm);
      }
    }
    std::string name_;
    std::string toReplace_;
    int idxToReplace_;
    TFormula* replaceBy_;
    TFormula* deltaFuncDerrivative_;
    std::vector<replaceParBase*> parForReplacements_;
    int numParForReplacements_;
    std::vector<replaceParBase*> parForDeltaFuncDerrivative_;
    int numParForDeltaFuncDerrivative_;
  };

  std::vector<fitParameterReplacementType*> fitParameterReplacements_;

  std::vector<SVfitParameterMappingType> constParameterMappings_;
  unsigned numConstParameters_;

  ROOT::Math::Functor* integrand_;
  MarkovChainIntegrator* integrator_;
  ROOT::Math::Functor* auxPhysicalSolutionFinder_;
  mutable std::vector<double> startPosition_;
  mutable std::vector<double> intBoundaryLower_;
  mutable std::vector<double> intBoundaryUpper_;
  unsigned numDimensions_;

  bool monitorMarkovChain_;  
  std::string monitorFilePath_;
  std::vector<ROOT::Math::Functor*> auxResonance_or_DaughterValues_;

  double* fitParameterValues_;

  mutable std::vector<TH1*> probHistFitParameter_;
  struct AuxProbHistogramsDaughter
  {
    size_t idxDaughter_;
    TH1* probHistDaughterPt_;
    TH1* probHistDaughterEta_;
    TH1* probHistDaughterPhi_;
  };
  struct AuxProbHistogramsResonance
  {
    size_t idxResonance_;
    TH1* probHistResonancePt_;
    TH1* probHistResonanceEta_;
    TH1* probHistResonancePhi_;
    TH1* probHistResonanceMass_;
    std::vector<AuxProbHistogramsDaughter*> probHistDaughters_;
  };
  mutable std::vector<AuxProbHistogramsResonance*> probHistResonances_;
  mutable TH1* probHistEventMass_;
  mutable std::vector<double> probListEventMass_;
  ROOT::Math::Functor* auxFillProbHistograms_;
  int max_or_median_;

  bool applyJacobiFactors_;
};

#endif

