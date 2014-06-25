#include "TauAnalysis/SVfit/interface/svFitAuxFunctions.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Vector3D.h"

#include <TMath.h>
#include <Math/VectorUtil.h>
#include <TMatrixD.h>
#include <TVectorD.h>

#include <limits>

namespace svFit_namespace
{
  // Adapted for our vector types from TVector3 class
  reco::Candidate::Vector rotateUz(const math::RThetaPhiVector& toRotate, const reco::Candidate::Vector& newUzVector)
  {
    // NB: newUzVector must be a unit vector !
    Double_t u1 = newUzVector.X();
    Double_t u2 = newUzVector.Y();
    Double_t u3 = newUzVector.Z();
    Double_t up = u1*u1 + u2*u2;

    Double_t fX = toRotate.X();
    Double_t fY = toRotate.Y();
    Double_t fZ = toRotate.Z();

    if ( up > 0. ) {
      up = TMath::Sqrt(up);
      Double_t px = fX;
      Double_t py = fY;
      Double_t pz = fZ;
      fX = (u1*u3*px - u2*py + u1*up*pz)/up;
      fY = (u2*u3*px + u1*py + u2*up*pz)/up;
      fZ = (u3*u3*px -    px + u3*up*pz)/up;
    } else if ( u3 < 0. ) {
      fX = -fX;
      fZ = -fZ;
    } else {}; // phi = 0, theta = pi

    return reco::Candidate::Vector(fX, fY, fZ);
  }

  reco::Candidate::LorentzVector boostToCOM(
      const reco::Candidate::LorentzVector& comSystem,
      const reco::Candidate::LorentzVector& p4ToBoost) {
    reco::Candidate::Vector boost = comSystem.BoostToCM();
    return ROOT::Math::VectorUtil::boost(p4ToBoost, boost);
  }

  reco::Candidate::LorentzVector boostToLab(
      const reco::Candidate::LorentzVector& rfSystem,
      const reco::Candidate::LorentzVector& p4ToBoost) {
    reco::Candidate::Vector boost = rfSystem.BoostToCM();
    return ROOT::Math::VectorUtil::boost(p4ToBoost, -boost);
  }

  double gjAngleLabFrameFromX(double x, double visMass, double invisMass, double pVis_lab, double enVis_lab, double motherMass, bool& isValidSolution) 
  {
    // CV: the expression for the Gottfried-Jackson angle as function of X = Etau/Evis
    //     was obtained by solving equation (1) of AN-2010/256:
    //       http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2010_256_v2.pdf
    //     for cosThetaGJ
    //    (generalized to the case of non-zero mass of the neutrino system in leptonic tau decays, using Mathematica)

    //std::cout << "<gjAngleLabFrameFromX>:" << std::endl;
    //std::cout << " x = " << x << std::endl;
    //std::cout << " visMass = " << visMass << std::endl;
    //std::cout << " invisMass = " << invisMass << std::endl;
    //std::cout << " pVis_lab = " << pVis_lab << std::endl;
    //std::cout << " enVis_lab = " << enVis_lab << std::endl;
    //std::cout << " motherMass = " << motherMass << std::endl;

    double x2 = x*x;
    double visMass2 = visMass*visMass;
    double invisMass2 = invisMass*invisMass;
    double pVis2_lab = pVis_lab*pVis_lab;
    double enVis2_lab = enVis_lab*enVis_lab;
    double motherMass2 = motherMass*motherMass;
    double term1 = enVis2_lab - motherMass2*x2;
    double term2 = 2.*TMath::Sqrt(pVis2_lab*enVis2_lab*enVis2_lab*term1);
    double term3 = ((visMass2 - invisMass2) + motherMass2)*pVis_lab*x*TMath::Sqrt(term1);
    double term4 = 2.*pVis2_lab*term1;
    double cosGjAngle_lab1 =  (term2 - term3)/term4;
    //std::cout << "cosGjAngle_lab1 = " << cosGjAngle_lab1 << std::endl;
    double cosGjAngle_lab2 = -(term2 + term3)/term4;
    //std::cout << "cosGjAngle_lab2 = " << cosGjAngle_lab2 << std::endl;
    double gjAngle = 0.;
    if ( TMath::Abs(cosGjAngle_lab1) <= 1. && TMath::Abs(cosGjAngle_lab2) > 1. ) {
      //std::cout << "taking solution '+-'" << std::endl;
      gjAngle = TMath::ACos(cosGjAngle_lab1);
    } else if ( TMath::Abs(cosGjAngle_lab1) > 1. && TMath::Abs(cosGjAngle_lab2) <= 1. ) {
      //std::cout << "taking solution '--'" << std::endl;
      gjAngle = TMath::ACos(cosGjAngle_lab2);
    } else if ( TMath::Abs(cosGjAngle_lab1) <= 1. && TMath::Abs(cosGjAngle_lab2) <= 1. ) {
      edm::LogWarning ("gjAngleFromX_new")
	<< " Failed to find unique solution !!";
      isValidSolution = false;
    } else {
      //edm::LogWarning ("gjAngleFromX_new")
      //  << " Failed to find valid solution !!";
      isValidSolution = false;
    }
    //std::cout << "--> gjAngle = " << gjAngle << std::endl;

    return gjAngle;
  }

  double gjAngleLabFrame_max(double visMass, double invisMass, double pVis_lab, double motherMass)
  {
    //std::cout << "<gjAngleLabFrame_max>:" << std::endl;
    
    double visMass2 = visMass*visMass;
    double invisMass2 = invisMass*invisMass;
    double motherMass2 = motherMass*motherMass;

    // CV: computation of maximum Gottfried-Jackson angle not fully debugged yet !!    
    double gjAngleMax_lab = TMath::ASin(TMath::Sqrt(square(invisMass2) + square(motherMass2 - visMass2) - 2.*invisMass2*(motherMass2 + visMass2))/(2.*motherMass*pVis_lab));
    //std::cout << "gjAngleMax_lab = " << gjAngleMax_lab << std::endl;

    return gjAngleMax_lab;
  }

  std::pair<double, double> gjAngleLabFrameToX(double gjAngle_lab, double visMass, double invisMass, double pVis_lab, double enVis_lab, double motherMass, bool& isValidSolution)
  {
    // CV: the expression for the Gottfried-Jackson angle as function of X = Etau/Evis
    //     was obtained by solving equation (1) of AN-2010/256:
    //       http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2010_256_v2.pdf
    //     for cosThetaGJ
    //    (generalized to the case of non-zero mass of the neutrino system in leptonic tau decays, using Mathematica)

    //std::cout << "<gjAngleLabFrameToX>:" << std::endl;
    //std::cout << " gjAngle_lab = " << gjAngle_lab << std::endl;
    //std::cout << " visMass = " << visMass << std::endl;
    //std::cout << " invisMass = " << invisMass << std::endl;
    //std::cout << " pVis_lab = " << pVis_lab << std::endl;
    //std::cout << " enVis_lab = " << enVis_lab << std::endl;
    //std::cout << " motherMass = " << motherMass << std::endl;

    double visMass2 = visMass*visMass;
    double invisMass2 = invisMass*invisMass;
    double pVis2_lab = pVis_lab*pVis_lab;
    double enVis2_lab = enVis_lab*enVis_lab;
    double motherMass2 = motherMass*motherMass;

    double gjAngle_max_lab = gjAngleLabFrame_max(visMass, invisMass, pVis_lab, motherMass);
    if ( gjAngle_lab > gjAngle_max_lab ) {
      //edm::LogWarning ("gjAngleToX")
      //  << " Failed to find valid solution !!";
      gjAngle_lab = gjAngle_max_lab;
      isValidSolution = false;
    }
    double cosGjAngle_lab = TMath::Cos(gjAngle_lab);
    double sinGjAngle_lab = TMath::Sin(gjAngle_lab);

    double term1 = (motherMass2 + visMass2) - invisMass2;
    //std::cout << "term1 = " << term1 << std::endl;
    double term2 = 2.*term1*enVis2_lab;
    //std::cout << "term2 = " << term2 << std::endl;
    double term3 = square(cosGjAngle_lab)*pVis2_lab*enVis2_lab*(square(invisMass2) + square(motherMass2 - visMass2) - 2.*invisMass2*(motherMass2 + visMass2) - 4.*square(sinGjAngle_lab)*motherMass2*pVis2_lab);
    //std::cout << "term3 = " << term3 << std::endl;
    assert(term3 >= 0.);
    double term4 = square(term1) + 4.*square(cosGjAngle_lab)*motherMass2*pVis2_lab;
    //std::cout << "term4 = " << term4 << std::endl;
    double x_1 =  (term2 + 2.*TMath::Sqrt(term3))/term4;
    double x_2 =  (term2 - 2.*TMath::Sqrt(term3))/term4;
    double x_3 = -(term2 - 2.*TMath::Sqrt(term3))/term4;
    double x_4 = -(term2 + 2.*TMath::Sqrt(term3))/term4;
    //std::cout << "--> x_1 = " << x_1 << ", x_2 = " << x_2 << ", x_3 = " << x_3 << ", x_4 = " << x_4 << std::endl;

    double x_forward  = -1.;
    double x_backward = -1.;
    int numSolutions = 0;
    if ( x_1 >= 0. && x_1 <= 1. ) {
      x_forward = x_1;
      ++numSolutions;
    }
    if ( x_2 >= 0. && x_2 <= 1. ) {
      if ( numSolutions == 0 ) x_forward = x_2;
      else x_backward = x_2;
      ++numSolutions;
    }
    if ( x_3 >= 0. && x_3 <= 1. ) {
      if ( numSolutions == 0 ) x_forward = x_3;
      else x_backward = x_3;
      ++numSolutions;
    }
    if ( x_4 >= 0. && x_4 <= 1. ) {
      x_backward = x_4;
      ++numSolutions;
    }
    assert(numSolutions == 2);
    //if ( x_forward < x_backward ) {
    //  double temp = x_forward;
    //  x_forward  = x_backward;
    //  x_backward = temp;
    //}

    return std::pair<double, double>(x_forward, x_backward);
  }

  double pVisRestFrame(double visMass, double invisMass, double motherMass)
  {
    double motherMass2 = motherMass*motherMass;
    double pVis = TMath::Sqrt((motherMass2 - square(visMass + invisMass))
                             *(motherMass2 - square(visMass - invisMass)))/(2.*motherMass);
    return pVis;
  }

  reco::Candidate::Vector motherDirection(const reco::Candidate::Vector& pVisLabFrame, double angleVisLabFrame, double phiLab) 
  {
    // The direction is defined using polar coordinates in a system where the visible energy
    // defines the Z axis.
    math::RThetaPhiVector motherDirectionVisibleSystem(1.0, angleVisLabFrame, phiLab);

    // Rotate into the LAB coordinate system
    return rotateUz(motherDirectionVisibleSystem, pVisLabFrame.Unit());
  }

  double gjAngleRestFrameFromLabMomenta(const reco::Candidate::LorentzVector& motherP4, const reco::Candidate::LorentzVector& visP4)
  {
    double gjAngle_rf = 0.;
    reco::Candidate::LorentzVector visP4_rf = boostToCOM(motherP4, visP4);
    if ( (motherP4.pt()*visP4_rf.pt()) > 0. ) {
      double scalarProduct = (motherP4.px()*visP4_rf.px() 
                            + motherP4.py()*visP4_rf.py() 
                            + motherP4.pz()*visP4_rf.pz())/(motherP4.P()*visP4_rf.P());
      gjAngle_rf = TMath::ACos(scalarProduct);
    }
    return gjAngle_rf;
  }
 
  reco::Candidate::Vector normalize(const reco::Candidate::Vector& p)
  {
    double p_x = p.x();
    double p_y = p.y();
    double p_z = p.z();
    double mag2 = square(p_x) + square(p_y) + square(p_z);
    if ( mag2 <= 0. ) return p;
    double mag = TMath::Sqrt(mag2);
    return reco::Candidate::Vector(p_x/mag, p_y/mag, p_z/mag);
  }

  double compScalarProduct(const reco::Candidate::Vector& p1, const reco::Candidate::Vector& p2)
  {
    return (p1.x()*p2.x() + p1.y()*p2.y() + p1.z()*p2.z());
  }
  
  reco::Candidate::Vector compCrossProduct(const reco::Candidate::Vector& p1, const reco::Candidate::Vector& p2)
  {
    double p3_x = p1.y()*p2.z() - p1.z()*p2.y();
    double p3_y = p1.z()*p2.x() - p1.x()*p2.z();
    double p3_z = p1.x()*p2.y() - p1.y()*p2.x();
    return reco::Candidate::Vector(p3_x, p3_y, p3_z);
  }

  double phiLabFromLabMomenta(const reco::Candidate::LorentzVector& motherP4, const reco::Candidate::LorentzVector& visP4)
  {
    reco::Candidate::Vector u_z = normalize(reco::Candidate::Vector(visP4.px(), visP4.py(), visP4.pz()));
    reco::Candidate::Vector u_y = normalize(compCrossProduct(reco::Candidate::Vector(0., 0., 1.), u_z));
    reco::Candidate::Vector u_x = compCrossProduct(u_y, u_z);
    
    reco::Candidate::Vector p3Mother_unit = normalize(reco::Candidate::Vector(motherP4.px(), motherP4.py(), motherP4.pz()));
    
    double phi_lab = TMath::ATan2(compScalarProduct(p3Mother_unit, u_y), compScalarProduct(p3Mother_unit, u_x));
    return phi_lab;
  }

  //
  //-------------------------------------------------------------------------------
  //

  void printVector(const std::string& label, const AlgebraicVector3& p)
  {
    reco::Candidate::Vector v(p(0), p(1), p(2));
    std::cout << label << ": x = " << v.x() << ", y = " << v.y() << ", z = " << v.z() << " (eta = " << v.eta() << ", phi = " << v.phi() << ")" << std::endl;
  }

  void printVector(const std::string& label, const AlgebraicVector2& p)
  {
    reco::Candidate::Vector v(p(0), p(1), 0.);
    std::cout << label << ": x = " << v.x() << ", y = " << v.y() << " (phi = " << v.phi() << ")" << std::endl;
  }

  void printMatrix(const std::string& label, const AlgebraicMatrix33& m)
  {
    std::cout << label << ":" << std::endl;
    for ( unsigned iRow = 0; iRow < 3; ++iRow ) {
      std::cout << " |";
      for ( unsigned iColumn = 0; iColumn < 3; ++iColumn ) {
	std::cout << " " << std::setw(12) << m(iRow, iColumn);
      }
      std::cout << " |" << std::endl;
    }
  }

  void printMatrix(const std::string& label, const AlgebraicMatrix22& m)
  {
    std::cout << label << ":" << std::endl;
    for ( unsigned iRow = 0; iRow < 2; ++iRow ) {
      std::cout << " |";
      for ( unsigned iColumn = 0; iColumn < 2; ++iColumn ) {
	std::cout << " " << std::setw(12) << m(iRow, iColumn);
      }
      std::cout << " |" << std::endl;
    }
  }

  //
  //-------------------------------------------------------------------------------
  //

  double logGaussian(double residual, double sigma)
  {
    if ( sigma > 0. ) {
      return -0.5*(TMath::Log(2*TMath::Pi()*square(sigma)) + square(residual/sigma));
    } else {
      edm::LogError ("logGaussian")
	<< " Parameter sigma must not be zero !!";
      return std::numeric_limits<float>::min();
    }
  }
  
  //
  //-------------------------------------------------------------------------------
  //

  double compScalarProduct(const AlgebraicVector3& p1, const AlgebraicVector3& p2)
  {
    return (p1(0)*p2(0) + p1(1)*p2(1) + p1(2)*p2(2));
  }

  AlgebraicVector3 compCrossProduct(const AlgebraicVector3& p1, const AlgebraicVector3& p2)
  {    
    double p3_x = p1(1)*p2(2) - p1(2)*p2(1);
    double p3_y = p1(2)*p2(0) - p1(0)*p2(2);
    double p3_z = p1(0)*p2(1) - p1(1)*p2(0);
    return AlgebraicVector3(p3_x, p3_y, p3_z);
  }

  double norm2(const AlgebraicVector3& p)
  {
    return compScalarProduct(p, p);
  }

  double norm2(const AlgebraicVector2& p)
  {
    return (p(0)*p(0) + p(1)*p(1));
  }

  double phi(const AlgebraicVector3& p)
  {
    return TMath::ATan2(p(1), p(0));
  }

  AlgebraicVector3 normalize(const AlgebraicVector3& p)
  {
    double p_mag2 = norm2(p);
    if ( p_mag2 > 0. ) {
      double p_mag = TMath::Sqrt(p_mag2);
      return AlgebraicVector3(p(0)/p_mag, p(1)/p_mag, p(2)/p_mag);
    } else {
      return AlgebraicVector3(0., 0., 0.);
    }
  }

  AlgebraicVector3 compDecayPosition_helix(const AlgebraicVector3& origin, double s, const reco::Candidate::LorentzVector& p4, double charge)
  {
    //std::cout << "<compDecayPosition_helix>:" << std::endl;
    //std::cout << " origin: x = " << origin(0) << ", y = " << origin(1) << ", z = " << origin(2) << std::endl;
    //std::cout << " s = " << s << std::endl;
    //std::cout << " p4: E = " << p4.E() << ", Pt = " << p4.pt() << ", eta = " << p4.eta() << ", phi = " << p4.phi() 
    //	        << " (px = " << p4.px() << ", py = " << p4.py() << ", pz = " << p4.pz() << ")" << std::endl;
    //std::cout << " charge = " << charge << std::endl;

    const double B = 3.8;   // Tesla
    const double c = 3.e+8; // m/s
    const double a = 1.e-11*c*B;
    //std::cout << " a = " << a << std::endl;
    double R = p4.pt()/a;
    //std::cout << " R = " << R << std::endl;
    double sign_of_kappa = 0.;
    if      ( charge > +0.5 ) sign_of_kappa = -1.;
    else if ( charge < -0.5 ) sign_of_kappa = +1.;
    double one_div_kappa = sign_of_kappa*R;
    double kappa = ( one_div_kappa != 0. ) ? (1./one_div_kappa) : 0.;
    double phi = p4.phi();
    double x0 = origin(0) - one_div_kappa*TMath::Sin(phi);
    double y0 = origin(1) + one_div_kappa*TMath::Cos(phi);
    double z0 = origin(2);
    double s_xy = s*TMath::Sin(p4.theta());
    double s_z = s*TMath::Cos(p4.theta());

    AlgebraicVector3 retVal = AlgebraicVector3(x0 + one_div_kappa*TMath::Sin(phi + kappa*s_xy), y0 - one_div_kappa*TMath::Cos(phi + kappa*s_xy), z0 + s_z);
    //std::cout << "--> retVal: x = " << retVal(0) << ", y = " << retVal(1) << ", z = " << retVal(2) << " (phi = " << TMath::ATan2(retVal(1), retVal(0)) << ")" << std::endl;
    //AlgebraicVector3 tauFlightPath_unit = normalize(AlgebraicVector3(p4.px(), p4.py(), p4.pz()));
    //compDecayPosition_line(origin, s, tauFlightPath_unit);

    return retVal;
  }

  AlgebraicVector3 compDecayPosition_line(const AlgebraicVector3& origin, double s, const AlgebraicVector3& tauFlightPath_unit)
  {
    //std::cout << "<compDecayPosition_line>:" << std::endl;
    //std::cout << " origin: x = " << origin(0) << ", y = " << origin(1) << ", z = " << origin(2) << std::endl;
    //std::cout << " s = " << s << std::endl;
    //std::cout << " tauFlightPath_unit: x = " << tauFlightPath_unit(0) << ", y = " << tauFlightPath_unit(1) << ", z = " << tauFlightPath_unit(2) << std::endl;

    AlgebraicVector3 retVal = AlgebraicVector3(origin(0) + s*tauFlightPath_unit(0), origin(1) + s*tauFlightPath_unit(1), origin(2) + s*tauFlightPath_unit(2));
    //std::cout << "--> retVal: x = " << retVal(0) << ", y = " << retVal(1) << ", z = " << retVal(2) << " (phi = " << TMath::ATan2(retVal(1), retVal(0)) << ")" << std::endl;

    return retVal;
  }
  
  void compLocalCoordinates(const AlgebraicVector3& direction, AlgebraicVector3& u1, AlgebraicVector3& u2, AlgebraicVector3& u3)
  {
    u3 = normalize(direction);
    u2 = normalize(compCrossProduct(AlgebraicVector3(0., 0., 1.), u3));
    u1 = compCrossProduct(u2, u3);
  }

  AlgebraicVector3 transformToLocalCoordinatesPos(const AlgebraicVector3& p, const AlgebraicVector3& u1, const AlgebraicVector3& u2, const AlgebraicVector3& u3)
  {
    double c1 = compScalarProduct(p, u1);
    double c2 = compScalarProduct(p, u2);
    double c3 = compScalarProduct(p, u3);
    return AlgebraicVector3(c1, c2, c3);    
  }

  AlgebraicMatrix33 transformToLocalCoordinatesCov(const AlgebraicMatrix33& cov, const AlgebraicVector3& u1, const AlgebraicVector3& u2, const AlgebraicVector3& u3)
  {
    AlgebraicMatrix33 R;    
    for ( int i = 0; i < 3; ++i ) {
      R(0, i) = u1(i);
      R(1, i) = u2(i);
      R(2, i) = u3(i);
    }
    AlgebraicMatrix33 R_T;
    for ( int i = 0; i < 3; ++i ) {
      for ( int j = 0; j < 3; ++j ) {
	R_T(i, j) = R(j, i);
      }
    }
    AlgebraicMatrix33 R_times_cov_times_R_T = R*cov*R_T;
    return R_times_cov_times_R_T;
  }

  void extractEigenValues(const AlgebraicMatrix33& cov, double& EigenValue1, double& EigenValue2, double& EigenValue3)
  {
    TMatrixD covMatrix(3, 3);
    for ( int iRow = 0; iRow < 3; ++iRow ) {
      for ( int iColumn = 0; iColumn < 3; ++iColumn ) {
	covMatrix(iRow, iColumn) = cov(iRow, iColumn);
      }
    }
    
    TVectorD EigenValues;
    covMatrix.EigenVectors(EigenValues);

    EigenValue1 = EigenValues(0); // largest EigenValue
    EigenValue2 = EigenValues(1); 
    EigenValue3 = EigenValues(2); // smallest EigenValue
  }

  double logGaussian2d(const AlgebraicVector2& residual, const AlgebraicMatrix22& cov)
  {
    // Warning: This function is specific to the computation of track likelihoods.
    //          Residual and covariance matrix are assumed to be given in a coordinate system
    //          in which the Eigenvector with the largest Eigenvalue of the covariance matrix
    //          corresponds to the 3rd component.
    //          Following the suggestion of John Smith in
    //            https://hypernews.cern.ch/HyperNews/CMS/get/recoTracking/1262/1/1/1.html
    //          the 3rd component is ignored when computing track likelihood values, 
    //          as the uncertainties of the track extrapolation in this direction are ill-defined.

    double det = cov(0, 0)*cov(1, 1) - cov(0, 1)*cov(1, 0);
    if ( det != 0. ) {
      double covInverse00 =  cov(1, 1)/det;
      double covInverse01 = -cov(0, 1)/det;
      double covInverse10 = -cov(1, 0)/det;
      double covInverse11 =  cov(0, 0)/det;
      double pull =  residual(0)*(covInverse00*residual(0) + covInverse01*residual(1))
        	   + residual(1)*(covInverse10*residual(0) + covInverse11*residual(1));
      return -TMath::Log(2.*TMath::Pi()*TMath::Sqrt(det)) - 0.5*pull;
    } else {
       edm::LogError ("logGaussian2d")
	 << " Cannot invert 2x2 covariance matrix, det = " << det << " !!";
       return std::numeric_limits<float>::min();
    }
  }

  double logGaussian3d(const AlgebraicVector3& residual, const AlgebraicMatrix33& cov)
  {
    double det;
    bool statusFlag = cov.Det2(det); // True = calculation successful
    if ( det != 0. && statusFlag ) {
      AlgebraicMatrix33 covInverse = cov;
      bool statusFlag = covInverse.Invert();
      if ( statusFlag ) {
        double pull =  residual(0)*(covInverse(0, 0)*residual(0) + covInverse(0, 1)*residual(1) + covInverse(0, 2)*residual(2))
                     + residual(1)*(covInverse(1, 0)*residual(0) + covInverse(1, 1)*residual(1) + covInverse(1, 2)*residual(2))
                     + residual(2)*(covInverse(2, 0)*residual(0) + covInverse(2, 1)*residual(1) + covInverse(2, 2)*residual(2));      
        return -1.5*TMath::Log(2.*TMath::Pi()) - 0.5*TMath::Log(det) - 0.5*pull;
      }
    }

    edm::LogError ("logGaussian3d")
      << " Cannot invert 3x3 covariance matrix, det = " << det << " !!";
    cov.Print(std::cout);
    return std::numeric_limits<float>::min();
  }

  //
  //-------------------------------------------------------------------------------
  //

  TH1* compHistogramDensity(const TH1* histogram)
  {
    std::string histogramName_density = std::string(histogram->GetName()).append("_density");
    TH1* histogram_density = (TH1*)histogram->Clone(histogramName_density.data());
    for ( int iBin = 1; iBin <= histogram->GetNbinsX(); ++iBin ) {
      double binContent = histogram->GetBinContent(iBin);
      double binError = histogram->GetBinError(iBin);
      double binWidth = histogram->GetBinWidth(iBin);
      assert(binWidth > 0.);
      histogram_density->SetBinContent(iBin, binContent/binWidth);
      histogram_density->SetBinError(iBin, binError/binWidth);
    }
    return histogram_density;
  }

  double getMeanOfBinsAboveThreshold(const TH1* histogram, double threshold, int verbosity)
  {
    //std::cout << "<getMeanOfBinsAboveThreshold>:" << std::endl;
    //std::cout << " threshold = " << threshold << std::endl;

    double mean = 0.;
    double normalization = 0.;
    int numBins = histogram->GetNbinsX();
    for ( int iBin = 1; iBin <= numBins; ++iBin ) {
      double binCenter = histogram->GetBinCenter(iBin);
      double binContent = histogram->GetBinContent(iBin);
      if ( binContent >= threshold ) {
	if ( verbosity >= 3 ) std::cout << " adding binContent = " << binContent << " @ binCenter = " << binCenter << std::endl;
	mean += (binCenter*binContent);
	normalization += binContent;
      }
    }
    if ( normalization > 0. ) mean /= normalization;
    if ( verbosity >= 3 ) std::cout << "--> mean = " << mean << std::endl;
    return mean;
  }

  void extractHistogramProperties(const TH1* histogram, const TH1* histogram_density,
				  double& xMaximum, double& xMaximum_interpol, 
				  double& xMean,
				  double& xQuantile016, double& xQuantile050, double& xQuantile084,
				  double& xMean3sigmaWithinMax, double& xMean5sigmaWithinMax, 
				  int verbosity)
  {
//--- compute median, -1 sigma and +1 sigma limits on reconstructed mass
    if ( verbosity >= 3 ) std::cout << "<extractHistogramProperties>:" << std::endl;

    if ( histogram->Integral() > 0. ) {
      Double_t q[3];
      Double_t probSum[3];
      probSum[0] = 0.16;
      probSum[1] = 0.50;
      probSum[2] = 0.84;
      (const_cast<TH1*>(histogram))->GetQuantiles(3, q, probSum);
      xQuantile016 = q[0];
      xQuantile050 = q[1];
      xQuantile084 = q[2];
    } else {
      xQuantile016 = 0.;
      xQuantile050 = 0.;
      xQuantile084 = 0.;
    }
    
    xMean = histogram->GetMean();
    
    if ( histogram_density->Integral() > 0. ) {
      int binMaximum = histogram_density->GetMaximumBin();
      xMaximum = histogram_density->GetBinCenter(binMaximum);
      double yMaximum = histogram_density->GetBinContent(binMaximum);
      double yMaximumErr = ( histogram->GetBinContent(binMaximum) > 0. ) ?	
	(yMaximum*histogram->GetBinError(binMaximum)/histogram->GetBinContent(binMaximum)) : 0.;
      if ( verbosity >= 3 ) std::cout << "yMaximum = " << yMaximum << " +/- " << yMaximumErr << " @ xMaximum = " << xMaximum << std::endl;
      if ( binMaximum > 1 && binMaximum < histogram_density->GetNbinsX() ) {
	int binLeft       = binMaximum - 1;
	double xLeft      = histogram_density->GetBinCenter(binLeft);
	double yLeft      = histogram_density->GetBinContent(binLeft);    
	
	int binRight      = binMaximum + 1;
	double xRight     = histogram_density->GetBinCenter(binRight);
	double yRight     = histogram_density->GetBinContent(binRight); 
	
	double xMinus     = xLeft - xMaximum;
	double yMinus     = yLeft - yMaximum;
	double xPlus      = xRight - xMaximum;
	double yPlus      = yRight - yMaximum;
	
	xMaximum_interpol = xMaximum + 0.5*(yPlus*square(xMinus) - yMinus*square(xPlus))/(yPlus*xMinus - yMinus*xPlus);
      } else {
	xMaximum_interpol = xMaximum;
      }
      if ( verbosity >= 3 ) std::cout << "computing xMean3sigmaWithinMax:" << std::endl;
      xMean3sigmaWithinMax = getMeanOfBinsAboveThreshold(histogram_density, yMaximum - 3.*yMaximumErr, verbosity);
      if ( verbosity >= 3 ) std::cout << "computing xMean5sigmaWithinMax:" << std::endl;
      xMean5sigmaWithinMax = getMeanOfBinsAboveThreshold(histogram_density, yMaximum - 5.*yMaximumErr, verbosity);
    } else {
      xMaximum = 0.;
      xMaximum_interpol = 0.;
      xMean3sigmaWithinMax = 0.;
      xMean5sigmaWithinMax = 0.;
    }
  }
}
