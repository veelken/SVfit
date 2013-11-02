#include "TauAnalysis/SVfit/interface/SVfitParameter.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "TauAnalysis/SVfit/interface/svFitAuxFunctions.h"

#include <Math/Minimizer.h>

#include <algorithm>
#include <iostream>
#include <iomanip>

using namespace SVfit_namespace;

std::vector<double> SVfitParameter::defaultInitialValues_;
std::vector<std::pair<double, double> > SVfitParameter::defaultLimits_;
std::vector<double> SVfitParameter::defaultStepSizes_;
bool SVfitParameter::defaultValues_initialized_ = false;

SVfitParameter::SVfitParameter(int idx, const std::string& name, int type,
			       double initialValue, double lowerLimit, double upperLimit, double stepSize, bool isFixed)
  : idx_(idx),
    name_(name),
    type_(type),
    initialValue_(initialValue),
    lowerLimit_(lowerLimit),
    upperLimit_(upperLimit),
    stepSize_(stepSize),
    isFixed_(isFixed)
{
  if ( !(type >= 0 && type < svFit_namespace::kNumFitParameter) ) 
    throw cms::Exception("SVfitParameter")
      << "Invalid type = " << type << " !!\n";
  
  reset();
}

SVfitParameter::SVfitParameter(int idx, const std::string& name, int type, bool isFixed)
  : idx_(idx),
    name_(name),
    type_(type),
    isFixed_(isFixed)
{
  if ( !defaultValues_initialized_ ) initializeDefaultValues();

  if ( !(type >= 0 && type < svFit_namespace::kNumFitParameter) ) 
    throw cms::Exception("SVfitParameter")
      << "Invalid type = " << type << " !!\n";
  
  initialValue_ = defaultInitialValues_[type];
  lowerLimit_ = defaultLimits_[type].first;
  upperLimit_ = defaultLimits_[type].second;
  stepSize_ = defaultStepSizes_[type];   

  reset();
}

bool isLimitDisabled(double limit)
{
  if ( limit == TMath::QuietNaN() ) return true;
  else return false;
}

void SVfitParameter::dump(std::ostream& stream) const 
{
  stream << "param #" << idx_ << " (name = " << get_name_incl_type(name_, type_) << "):" 
	 << " value = " << value_ << " + " << errUp_ << " - " << errDown_ << std::endl;
  stream << " initialValue = " << initialValue_ << " +/- " << stepSize_ << ","
	 << " limits = {";
  if ( isLimitDisabled(lowerLimit_) ) stream << "disabled";
  else stream << lowerLimit_;
  stream << ", ";
  if ( isLimitDisabled(upperLimit_) ) stream << "disabled";
  else stream << upperLimit_;
  stream << "}";
  if      ( IsFixed()       ) stream << ", isFixed"; 
  else if ( IsDoubleBound() ) stream << ", isDoubleBound"; 
  else if ( HasLowerLimit() ) stream << ", hasLowerLimit"; 
  else if ( HasUpperLimit() ) stream << ", hasUpperLimit"; 
  stream << std::endl;
}

// Friend helpers for print functions
std::ostream& operator<<(std::ostream& stream, const SVfitParameter& param) 
{
  param.dump(stream);
  return stream;
}

std::string get_name_incl_type(const std::string& name, int type)
{
  std::string retVal = name;
  retVal.append("::");
  if      ( type == svFit_namespace::kPV_shiftX                   ) retVal.append("PV_shiftX");
  else if ( type == svFit_namespace::kPV_shiftY                   ) retVal.append("PV_shiftY");
  else if ( type == svFit_namespace::kPV_shiftZ                   ) retVal.append("PV_shiftZ");
  else if ( type == svFit_namespace::kTau_visEnFracX              ) retVal.append("visEnFracX");
  else if ( type == svFit_namespace::kTau_phi_lab                 ) retVal.append("phi_lab");
  else if ( type == svFit_namespace::kTau_decayDistance_lab_shift ) retVal.append("decayDistance_lab_shift");
  else if ( type == svFit_namespace::kTau_visMass                 ) retVal.append("visMass");
  else if ( type == svFit_namespace::kTau_nuInvMass               ) retVal.append("nuInvMass");
  else if ( type == svFit_namespace::kLep_shiftEn                 ) retVal.append("shiftEn");
  else if ( type == svFit_namespace::kNu_energy_lab               ) retVal.append("energy_lab");
  else if ( type == svFit_namespace::kNu_phi_lab                  ) retVal.append("phi_lab");
  else if ( type == svFit_namespace::kW_theta_lab                 ) retVal.append("W_theta_lab");
  else if ( type == svFit_namespace::kW_phi_lab                   ) retVal.append("W_phi_lab");
  else if ( type == svFit_namespace::kW_mass                      ) retVal.append("W_mass");
  else retVal.append("undefined");
  return retVal;
}

//
//-------------------------------------------------------------------------------
//

void SVfitParameter::initializeDefaultValues()
{
  //std::cout << "<SVfitParameter::initializeDefaultValues>:" << std::endl;
  //std::cout << " numFitParameter = " << svFit_namespace::kNumFitParameter << std::endl;

  defaultInitialValues_.resize(svFit_namespace::kNumFitParameter);
  defaultInitialValues_[svFit_namespace::kPV_shiftX]                   =  0.;
  defaultInitialValues_[svFit_namespace::kPV_shiftY]                   =  0.;
  defaultInitialValues_[svFit_namespace::kPV_shiftZ]                   =  0.;
  defaultInitialValues_[svFit_namespace::kTau_visEnFracX]              =  0.75;
  defaultInitialValues_[svFit_namespace::kTau_phi_lab]                 =  0.;
  defaultInitialValues_[svFit_namespace::kTau_decayDistance_lab_shift] =  0.; 
  defaultInitialValues_[svFit_namespace::kTau_visMass]                 =  0.8; // GeV
  defaultInitialValues_[svFit_namespace::kTau_nuInvMass]               =  0.8; // GeV
  defaultInitialValues_[svFit_namespace::kLep_shiftEn]                 =  0.;
  defaultInitialValues_[svFit_namespace::kNu_energy_lab]               =  0.;
  defaultInitialValues_[svFit_namespace::kNu_phi_lab]                  =  0.;
  defaultInitialValues_[svFit_namespace::kW_theta_lab]                 =  0.50*TMath::Pi();
  defaultInitialValues_[svFit_namespace::kW_phi_lab]                   =  0.;
  defaultInitialValues_[svFit_namespace::kW_mass]                      = 80.399; // GeV

  defaultLimits_.resize(svFit_namespace::kNumFitParameter);
  defaultLimits_[svFit_namespace::kPV_shiftX]                          = pdouble(         -0.1,          +0.1); // cm
  defaultLimits_[svFit_namespace::kPV_shiftY]                          = pdouble(         -0.1,          +0.1); // cm
  defaultLimits_[svFit_namespace::kPV_shiftZ]                          = pdouble(         -2.,           +2.);  // cm
  defaultLimits_[svFit_namespace::kTau_visEnFracX]                     = pdouble(          0.,            1.);  // dimensionless
  defaultLimits_[svFit_namespace::kTau_phi_lab]                        = pdouble(-TMath::Pi(),  +TMath::Pi());  // rad
  defaultLimits_[svFit_namespace::kTau_decayDistance_lab_shift]        = pdouble(        -2.5,           +2.5); // cm
  defaultLimits_[svFit_namespace::kTau_visMass]                        = pdouble(chargedPionMass, tauLeptonMass); // GeV
  defaultLimits_[svFit_namespace::kTau_nuInvMass]                      = pdouble(          0., tauLeptonMass);  // GeV
  defaultLimits_[svFit_namespace::kLep_shiftEn]                        = pdouble(          0.,           10.);  // relative to measured lepton energy
  defaultLimits_[svFit_namespace::kNu_energy_lab]                      = pdouble(          0.,         1.e+3);  // GeV
  defaultLimits_[svFit_namespace::kNu_phi_lab]                         = pdouble(          0.,   TMath::Pi());  // rad
  defaultLimits_[svFit_namespace::kW_theta_lab]                        = pdouble(          0.,   TMath::Pi());  // rad
  defaultLimits_[svFit_namespace::kW_phi_lab]                          = pdouble(-TMath::Pi(),  +TMath::Pi());  // rad
  defaultLimits_[svFit_namespace::kW_mass]                             = pdouble(80.399 - 3.*2.085, 80.399 + 3.*2.085); // GeV

  defaultStepSizes_.resize(svFit_namespace::kNumFitParameter);
  defaultStepSizes_[svFit_namespace::kPV_shiftX]                       =  0.01;
  defaultStepSizes_[svFit_namespace::kPV_shiftY]                       =  0.01;
  defaultStepSizes_[svFit_namespace::kPV_shiftZ]                       =  0.01;
  defaultStepSizes_[svFit_namespace::kTau_visEnFracX]                  =  0.1;
  defaultStepSizes_[svFit_namespace::kTau_phi_lab]                     =  0.25;
  defaultStepSizes_[svFit_namespace::kTau_decayDistance_lab_shift]     =  0.01;
  defaultStepSizes_[svFit_namespace::kTau_visMass]                     =  0.1;
  defaultStepSizes_[svFit_namespace::kTau_nuInvMass]                   =  0.1;
  defaultStepSizes_[svFit_namespace::kLep_shiftEn]                     =  0.01;
  defaultStepSizes_[svFit_namespace::kNu_energy_lab]                   = 10.;
  defaultStepSizes_[svFit_namespace::kNu_phi_lab]                      =  0.25;
  defaultStepSizes_[svFit_namespace::kW_theta_lab]                     =  0.25;
  defaultStepSizes_[svFit_namespace::kW_phi_lab]                       =  0.25;
  defaultStepSizes_[svFit_namespace::kW_mass]                          =  1.;
 
  defaultValues_initialized_ = true;
}
