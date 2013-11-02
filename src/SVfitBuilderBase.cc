#include "TauAnalysis/SVfit/interface/SVfitBuilderBase.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "TauAnalysis/SVfit/interface/SVfitAlgorithmBase.h"

int SVfitBuilderBase::getFitParameterIdx(SVfitAlgorithmBase* algorithm, const std::string& name, int type, bool isOptional)
{
  SVfitParameter* fitParameter = algorithm->getFitParameter(name, type);
  if ( fitParameter ) {
    return fitParameter->index();
  } else {
    if ( isOptional ) return -1;
    else throw cms::Exception("SVfitBuilderBase::getFitParameterIdx")
      << " No fitParameter = " << get_name_incl_type(name, type) << " defined !!\n"
      << "--> Please check if there is a likelihood which constrains this fitParameter.\n";
  }
}
