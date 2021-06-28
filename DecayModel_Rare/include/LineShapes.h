#ifndef LINESHAPES_H
#define LINESHAPES_H 1

#include "HHLineShapeBase.h"

#include <complex>

//Breit-Wigner Lineshape derived class for relativistic resonances
class BWShape : public HHLineShapeBase {
  
 public:
  
  //Constructors
  BWShape();

  BWShape( const double mass,
           const double width,
           const double massmother,
           const double massdaug1,
           const double massdaug2 ); 
  
  std::complex<double> lineshape( const double mhh ) const override;

 private:
  
  double m_mR;
  double m_gR;
  double m_norm{1.0};
  double m_mother;
  double m_daug1;
  double m_daug2;

};

class BWCoupledShape : public HHLineShapeBase {

 public:
  
  //Constructors
  BWCoupledShape();
  
  BWCoupledShape( const double mK, 
		  const double gK, 
		  const double mR, 
		  const double gR, 
		  const std::complex<double> coupling,
		  const double massmother,
		  const double massdaug1, 
		  const double massdaug2) ;
  
  //evaluate
  std::complex<double> lineshape( const double mhh ) const override;
  
 private:
  
  double m_mK;
  double m_gK;
  double m_mR;
  double m_gR;
  
  double m_mother;
  double m_daug1;
  double m_daug2;
  std::complex<double> m_g;
  
  double m_norm{1.0};
};

class BW3CoupledShape : public HHLineShapeBase {

 public:
  
  //Constructors
  BW3CoupledShape();
  
  BW3CoupledShape( const double m1,
				  const double g1,
				  const double m2,
				  const double g2,
				  const double m3,
				  const double g3,
				  const std::complex<double> c12,
				  const std::complex<double> c13,
				  const double massmother,
				  const double massdaug1,
				  const double massdaug2);

  //evaluate
  std::complex<double> lineshape( const double mhh ) const override;
  
 private:
  
 
  double m_m1;
  double m_g1;
  double m_m2;
  double m_g2;
  double m_m3;
  double m_g3;
  
  double m_mother;
  double m_daug1;
  double m_daug2;

  std::complex<double> m_g12;
  std::complex<double> m_g13;
  
  double m_norm{1.0};
};

#endif
