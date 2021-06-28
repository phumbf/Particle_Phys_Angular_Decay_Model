#include "LineShapes.h"
#include "Parameters.h"

#include <cmath>
#include <complex>
#include <iostream>

#include "TRandom.h"

//Breit-Wigner Lineshape derived class for relativistic resonances
BWShape::BWShape( const double mass,
	       	  const double width,
		  const double massmother,
		  const double massdaug1,
		  const double massdaug2 ) : 
  m_mR( mass  ),
  m_gR( width ),
  m_mother ( massmother ),
  m_daug1 ( massdaug1 ),
  m_daug2 ( massdaug2 )
{
  m_norm = normalise(m_mother,m_daug1,m_daug2) ; 
  std::cout << m_norm << std::endl;
}

std::complex<double> BWShape::lineshape( const double mhh ) const{
  
  std::complex<double> numer( m_norm, 0.0 );
  
  std::complex<double> 
    denom( std::pow( m_mR, 2 ) - std::pow( mhh, 2 ), -1.0*m_mR*m_gR );
  
  std::complex<double> bw = numer / denom;

  return bw;
}

BWCoupledShape::BWCoupledShape( const double mK, 
				const double gK, 
				const double mR, 
				const double gR, 
				const std::complex<double> coupling,
	      		        const double massmother,
	      		        const double massdaug1,
	      		        const double massdaug2 ) : 
  m_mK( mK ),
  m_gK( gK ), 
  m_mR( mR ), 
  m_gR( gR ), 
  m_g ( coupling ),
  m_mother ( massmother ),
  m_daug1 ( massdaug1 ),
  m_daug2 ( massdaug2 )
{
  m_norm = normalise(m_mother,m_daug1,m_daug2) ; 
}
  
std::complex<double> BWCoupledShape::lineshape( const double mhh ) const{

  const double msq{mhh*mhh};

  std::complex<double> 
    dK( std::pow( m_mK, 2 ) - msq - 0.25*std::pow( m_gK, 2 ), -1.0*m_mK*m_gK );
  
  std::complex<double> 
    dR( std::pow( m_mR, 2 ) - msq - 0.25*std::pow( m_gR, 2 ), -1.0*m_mR*m_gR ); 
  
  std::complex<double> bwcs =  ( 1.0 / dR ) - ( m_g / dK );
  
  bwcs *= m_norm;
  
  return bwcs;
}

BW3CoupledShape::BW3CoupledShape( const double m1,
				  const double g1,
				  const double m2,
				  const double g2,
				  const double m3,
				  const double g3,
				  const std::complex<double> c12,
				  const std::complex<double> c13,
				  const double massmother,
				  const double massdaug1,
				  const double massdaug2) :

  m_m1( m1 ),
  m_g1( g1 ),
  m_m2( m2 ),
  m_g2( g2 ),
  m_m3( m3 ),
  m_g3( g3 ),
  m_g12( c12 ),
  m_g13( c13 ),
  m_mother( massmother ),
  m_daug1( massdaug1 ),
  m_daug2( massdaug2 )
{
  m_norm = normalise(m_mother,m_daug1,m_daug2);
}

std::complex<double> BW3CoupledShape::lineshape( const double mhh) const{

	const double msq{mhh*mhh};

	std::complex<double> 
		d1( std::pow( m_m1, 2 ) - msq - 0.25*std::pow( m_g1, 2 ), -1.0*m_m1*m_g1 );

	std::complex<double> 
		d2( std::pow( m_m2, 2 ) - msq - 0.25*std::pow( m_g2, 2 ), -1.0*m_m2*m_g2 ); 

	std::complex<double> 
		d3( std::pow( m_m3, 2 ) - msq - 0.25*std::pow( m_g3, 2 ), -1.0*m_m3*m_g3 ); 

	std::complex<double> bwcs =  ( 1.0 / d1 ) - ( m_g12 / d2 ) - ( m_g13 / d3  );

	bwcs *= m_norm;

	return bwcs;
}
