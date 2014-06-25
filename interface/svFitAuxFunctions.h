#ifndef TauAnalysis_SVfit_svFitAuxFunctions_h
#define TauAnalysis_SVfit_svFitAuxFunctions_h

#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/GlobalError.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include <TH1.h>

namespace svFit_namespace
{
  //-----------------------------------------------------------------------------
  // define masses, widths and lifetimes of particles
  // relevant for computing values of likelihood functions in SVfit algorithm
  //
  // NOTE: the values are taken from
  //        K. Nakamura et al. (Particle Data Group),
  //        J. Phys. G 37, 075021 (2010)
  //
  const double electronMass = 0.51100e-3; // GeV
  const double electronMass2 = electronMass*electronMass;
  const double muonMass = 0.10566; // GeV
  const double muonMass2 = muonMass*muonMass; 
  
  const double chargedPionMass = 0.13957; // GeV
  const double chargedPionMass2 = chargedPionMass*chargedPionMass;
  const double neutralPionMass = 0.13498; // GeV
  const double neutralPionMass2 = neutralPionMass*neutralPionMass;

  const double rhoMesonMass = 0.77549; // GeV
  const double rhoMesonMass2 = rhoMesonMass*rhoMesonMass;
  const double rhoMesonWidth = 0.1491; // GeV

  const double a1MesonMass = 1.230; // GeV
  const double a1MesonMass2 = a1MesonMass*a1MesonMass;
  const double a1MesonWidth = 0.600; // GeV (upper limit of range quoted for width of a1 meson resonance in PDG summary tables)

  const double tauLeptonMass = 1.77685; // GeV
  const double tauLeptonMass2 = tauLeptonMass*tauLeptonMass;
  const double tauLeptonMass3 = tauLeptonMass2*tauLeptonMass;
  const double tauLeptonMass4 = tauLeptonMass3*tauLeptonMass;
  const double cTauLifetime = 8.711e-3; // centimeters

  const double mZ = 91.188; // GeV
  const double gammaZ = 2.495; // GeV

  const double alphaZ = 1./128.9; // fine-structure constant @ Z0 mass

  const double sinTheta_weinberg2 = 0.231;
  const double sinTheta_weinberg = TMath::Sqrt(sinTheta_weinberg2);
  const double cosTheta_weinberg = TMath::Sqrt(1. - sinTheta_weinberg2);
  const double qTau  = -1.;
  const double vTau  = (-1. + 4.*sinTheta_weinberg2)/(4.*sinTheta_weinberg*cosTheta_weinberg);     // -0.044
  const double aTau  = -1./(4.*sinTheta_weinberg*cosTheta_weinberg);                               // -0.593
  const double qUp   = +2./3;
  const double vUp   = (1. - (8./3)*sinTheta_weinberg2)/(4.*sinTheta_weinberg*cosTheta_weinberg);  //  0.227
  const double aUp   = 1./(4.*sinTheta_weinberg*cosTheta_weinberg);                                //  0.593
  const double qDown = -1./3;
  const double vDown = (-1. + (4./3)*sinTheta_weinberg2)/(4.*sinTheta_weinberg*cosTheta_weinberg); // -0.410
  const double aDown = -1./(4.*sinTheta_weinberg*cosTheta_weinberg);                               // -0.593
  //-----------------------------------------------------------------------------

  inline double square(double x)
  {
    return x*x;
  }

  inline double cube(double x)
  {
    return x*x*x;
  }

  inline double fourth(double x)
  {
    return x*x*x*x;
  }

  inline double fifth(double x)
  {
    return x*x*x*x*x;
  }

  /// Boost a lorentz vector given in the laboratory frame into the rest frame of another lorentz vector.
  reco::Candidate::LorentzVector boostToCOM(const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&);

  /// Boost a lorentz vector given in the rest frame of another lorentz vector to the laboratory frame.
  reco::Candidate::LorentzVector boostToLab(const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&);

  inline double energyFromMomentum(double momentum, double mass) {
    return TMath::Sqrt(square(mass) + square(momentum));
  }

  /// Determine Gottfried-Jackson angle from visible energy fraction X
  double gjAngleLabFrameFromX(double, double, double, double, double, double, bool&); 

  double gjAngleLabFrame_max(double, double, double, double);

  /// Determine visible energy fraction X from Gottfried-Jackson angle 
  /// (Note: the function returns 2 solutions, corresponding to "forward"/"backward" decays)
  std::pair<double, double> gjAngleLabFrameToX(double, double, double, double, double, double, bool&);

  /// Determine visible tau rest frame energy given visible mass and neutrino mass
  double pVisRestFrame(double, double, double);

  /// Determine the tau direction given our parameterization
  reco::Candidate::Vector motherDirection(const reco::Candidate::Vector&, double, double);

  /// Compute decay angles in rest frame given momentum of tau lepton and visible decay product in lab frame
  double gjAngleRestFrameFromLabMomenta(const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&);
  reco::Candidate::Vector normalize(const reco::Candidate::Vector&);
  double compScalarProduct(const reco::Candidate::Vector&, const reco::Candidate::Vector&);
  reco::Candidate::Vector compCrossProduct(const reco::Candidate::Vector&, const reco::Candidate::Vector&);
  double phiLabFromLabMomenta(const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&);
  
  void printVector(const std::string&, const AlgebraicVector3&);
  void printVector(const std::string&, const AlgebraicVector2&);
  void printMatrix(const std::string&, const AlgebraicMatrix33&);
  void printMatrix(const std::string&, const AlgebraicMatrix22&);

  /// Compute logarithm of Gaussian probability density function
  double logGaussian(double, double);

  /// Auxiliary functions for track likelihood analysis
  double compScalarProduct(const AlgebraicVector3&, const AlgebraicVector3&);
  AlgebraicVector3 compCrossProduct(const AlgebraicVector3&, const AlgebraicVector3&);
  double norm2(const AlgebraicVector3&);
  double norm2(const AlgebraicVector2&);
  double phi(const AlgebraicVector3&);
  AlgebraicVector3 normalize(const AlgebraicVector3&);
  AlgebraicVector3 compDecayPosition_helix(const AlgebraicVector3&, double, const reco::Candidate::LorentzVector&, double);
  AlgebraicVector3 compDecayPosition_line(const AlgebraicVector3&, double, const AlgebraicVector3&);
  void compLocalCoordinates(const AlgebraicVector3&, AlgebraicVector3&, AlgebraicVector3&, AlgebraicVector3&);
  AlgebraicVector3 transformToLocalCoordinatesPos(const AlgebraicVector3&, const AlgebraicVector3&, const AlgebraicVector3&, const AlgebraicVector3&);
  AlgebraicMatrix33 transformToLocalCoordinatesCov(const AlgebraicMatrix33&, const AlgebraicVector3&, const AlgebraicVector3&, const AlgebraicVector3&);
  void extractEigenValues(const AlgebraicMatrix33&, double&, double&, double&);
  double logGaussian2d(const AlgebraicVector2&, const AlgebraicMatrix22&);
  double logGaussian3d(const AlgebraicVector3&, const AlgebraicMatrix33&);

  /// Compute density of distribution (divide histogram by bin-width)
  TH1* compHistogramDensity(const TH1*);

  /// Extract maximum, mean and { 0.84, 0.50, 0.16 } quantiles of distribution
  void extractHistogramProperties(const TH1*, const TH1*, double&, double&, double&, double&, double&, double&, double&, double&, int = 0);
}

#endif
