#ifndef TauAnalysis_SVfit_SVfitResonanceBuilder_h
#define TauAnalysis_SVfit_SVfitResonanceBuilder_h

/** \class SVfitResonanceBuilder
 *
 * Auxiliary class for building SVfitResonanceHypothesis objects;
 * used by SVfit algorithm
 *
 * \author Christian Veelken, UC Davis
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/SVfit/interface/SVfitResonanceBuilderBase.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitResonanceHypothesis.h"

class SVfitResonanceBuilder : public SVfitResonanceBuilderBase
{
 public:
  SVfitResonanceBuilder(const edm::ParameterSet&); 
  ~SVfitResonanceBuilder() {}

  virtual SVfitResonanceHypothesis* build(const inputParticleMap&) const;

 private:
  /// different possible polarization states of each tau lepton pair 
  std::vector<int> polHandedness_;
  unsigned numPolStates_;
};

#endif


