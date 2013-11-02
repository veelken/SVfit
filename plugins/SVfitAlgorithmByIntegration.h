#ifndef TauAnalysis_SVfit_SVfitAlgorithmByIntegration_h
#define TauAnalysis_SVfit_SVfitAlgorithmByIntegration_h

/** \class SVfitAlgorithmByIntegration
 *
 * Concrete implementation of (n)SVfit algorithm
 * by integration of likelihood functions
 * (based on VEGAS integration algorithm)
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
#include "TauAnalysis/SVfit/interface/IndepCombinatoricsGeneratorT.h"
#include "TauAnalysis/SVfit/interface/svFitAuxFunctions.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitEventHypothesisByIntegration.h"
#include "AnalysisDataFormats/SVfit/interface/SVfitResonanceHypothesisByIntegration.h"

#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

#include <TArrayF.h>
#include <TFormula.h>
#include <TH1.h>
#include <TMath.h>

#include <vector>
#include <algorithm>
#include <string>

class SVfitAlgorithmByIntegration : public SVfitAlgorithmBase
{
 public:
  SVfitAlgorithmByIntegration(const edm::ParameterSet&);
  ~SVfitAlgorithmByIntegration();

  void beginJob();
  void beginEvent(const edm::Event&, const edm::EventSetup&);

  bool applyJacobiFactors() const 
  { 
    return applyJacobiFactors_; 
  }

  void print(std::ostream&) const {}

  double nll(const double*, const double*) const;

 protected:
  void fitImp() const;

  void setMassResults(SVfitResonanceHypothesisByIntegration*, const TH1*, unsigned) const;

  bool isDaughter(const std::string&);
  bool isResonance(const std::string&);

  struct replaceParBase
  {
    virtual void beginJob(SVfitAlgorithmByIntegration*) {}
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
    void beginJob(SVfitAlgorithmByIntegration* algorithm)
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
      : gridPoints_(0),
	resBinning_(0),
	replaceBy_(0),
	deltaFuncDerrivative_(0)
    {}
    ~fitParameterReplacementType() 
    {
      delete gridPoints_;
      delete resBinning_;
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
    void beginJob(SVfitAlgorithmByIntegration* algorithm)
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
    void beginEvent(double eventLowerLimit)
    {
      std::vector<double> gridPoints_vector;
      double gridPoint = TMath::Max(iterLowerLimit_, eventLowerLimit);
      bool isGridComplete = false;
      while ( !isGridComplete ) {
	if ( gridPoint >= iterUpperLimit_ ) isGridComplete = true;
	gridPoints_vector.push_back(gridPoint);
	double stepSize = TMath::Max((iterStepSizeFactor_ - 1.)*gridPoint, iterMinStepSize_);
	gridPoint += stepSize;	
      }
      std::sort(gridPoints_vector.begin(), gridPoints_vector.end());
      
      numGridPoints_ = gridPoints_vector.size();
      delete gridPoints_;
      gridPoints_ = new TArrayF(numGridPoints_);       
      delete resBinning_;
      resBinning_ = new TArrayF(numGridPoints_ + 1);
      for ( int iGridPoint = 0; iGridPoint < numGridPoints_; ++iGridPoint ) {
	double gridPoint = gridPoints_vector[iGridPoint];
	(*gridPoints_)[iGridPoint] = gridPoint;
	if   ( iGridPoint == 0 ) (*resBinning_)[0]          = gridPoint - 0.5*TMath::Abs(gridPoints_vector[1] - gridPoint);
	else                     (*resBinning_)[iGridPoint] = 0.5*(gridPoints_vector[iGridPoint - 1] + gridPoint);
	if   ( iGridPoint == (numGridPoints_ - 1) ) 
	  (*resBinning_)[numGridPoints_] = gridPoint + 0.5*TMath::Abs(gridPoint - gridPoints_vector[iGridPoint - 1]);
      }
    }
    std::string name_;
    double iterLowerLimit_;
    double iterUpperLimit_;
    double iterStepSizeFactor_;
    double iterMinStepSize_;
    int numGridPoints_;
    TArrayF* gridPoints_;
    TArrayF* resBinning_;
    std::string toReplace_;
    int idxToReplace_;
    TFormula* replaceBy_;
    TFormula* deltaFuncDerrivative_;
    int idxMassParameter_;
    std::vector<replaceParBase*> parForReplacements_;
    int numParForReplacements_;
    std::vector<replaceParBase*> parForDeltaFuncDerrivative_;
    int numParForDeltaFuncDerrivative_;
  };

  std::vector<fitParameterReplacementType*> fitParameterReplacements_;

  double* fitParameterValues_;

  double* xl_;
  double* xu_;

  gsl_monte_function* integrand_;
  gsl_monte_vegas_state* workspace_;
  mutable gsl_rng* rnd_;
  unsigned numCallsGridOpt_;
  unsigned numCallsIntEval_;
  double maxChi2_;
  unsigned maxIntEvalIter_;
  double precision_;
  unsigned numDimensions_;

  unsigned numMassParameters_;
  mutable IndepCombinatoricsGeneratorT<int>* massParForReplacements_;

  int max_or_median_;

  bool applyJacobiFactors_;
};

#endif

