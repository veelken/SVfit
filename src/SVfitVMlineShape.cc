#include "TauAnalysis/SVfit/interface/SVfitVMlineShape.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "TauAnalysis/SVfit/interface/svFitAuxFunctions.h"

using namespace svFit_namespace;

SVfitVMlineShape::SVfitVMlineShape(SVfitVMlineShapeIntegrand::VMtype vmType, SVfitVMlineShapeIntegrand::VMpol vmPol)
{
//--- compute lower limit for normalization integral
//   = invariant mass of n-pion system
  if      ( vmType == SVfitVMlineShapeIntegrand::kVMrho       ) minMass2_ = square(chargedPionMass + neutralPionMass);
  else if ( vmType == SVfitVMlineShapeIntegrand::kVMa1Neutral ) minMass2_ = square(chargedPionMass + 2*neutralPionMass);
  else if ( vmType == SVfitVMlineShapeIntegrand::kVMa1Charged ) minMass2_ = square(3*chargedPionMass);
  else throw cms::Exception("SVfitVMlineShape::SVfitVMlineShape")
    << " Invalid vmType = " << vmType << " !!\n";

  integrand_ = new SVfitVMlineShapeIntegrand(minMass2_);
  integrand_->SetVMtype(vmType);
  integrand_->SetVMpol(vmPol);
  integrand_->SetMode(SVfitVMlineShapeIntegrand::kVMnorm);

//--- CV: need to trigger update of ROOT::Math::Integrator by calling integrator->SetFunction
//        after calling any non-const function of SVfitVMlineShapeIntegrand 
  integrator_ = new ROOT::Math::Integrator(*integrand_);
  integrator_->SetFunction(*integrand_);

//--- compute vector meson line-shape normalization factor
  norm_ = integrator_->Integral(minMass2_, tauLeptonMass2); 

//--- set integrand to compute vector meson line-shape integrals in the following...  
  integrand_->SetMode(SVfitVMlineShapeIntegrand::kVMlineShape);
  integrator_->SetFunction(*integrand_);
}

SVfitVMlineShape::SVfitVMlineShape(const SVfitVMlineShape& bluePrint)
{
  integrand_ = new SVfitVMlineShapeIntegrand(*bluePrint.integrand_);
  
  integrator_ = new ROOT::Math::Integrator(*integrand_);
  integrator_->SetFunction(*integrand_);

  minMass2_ = bluePrint.minMass2_;
  norm_ = bluePrint.norm_;
}

SVfitVMlineShape::~SVfitVMlineShape()
{
  delete integrator_;
  delete integrand_;
}

SVfitVMlineShape& SVfitVMlineShape::operator=(const SVfitVMlineShape& bluePrint)
{
  delete integrator_;
  delete integrand_;

  integrand_ = new SVfitVMlineShapeIntegrand(*bluePrint.integrand_);
  
  integrator_ = new ROOT::Math::Integrator(*integrand_);
  integrator_->SetFunction(*integrand_);

  minMass2_ = bluePrint.minMass2_;
  norm_ = bluePrint.norm_;

  return (*this);
}

double SVfitVMlineShape::operator()(double theta, double tauLeptonPol, double z, double vmMass2) const
{
  if ( vmMass2 > minMass2_ && vmMass2 < (square(z)*tauLeptonMass2) ) {
    integrand_->SetParameterZ(z);
    integrand_->SetParameterTauLeptonPol(tauLeptonPol);
    return (*integrand_)(vmMass2)/norm_;
  } else {
    return 0.;
  }
}  
