/*
 * Model for the (k pi) lineshape using the LASS model
 * 
 * The model inherits most of its properties 
 * from KPiResonaceBase.
 */

// Local
#include "KPiResonanceLASS.h"
#include "HelperFunctions.h" 
#include "KinematicFunctions.h" 

KPiResonanceLASS::KPiResonanceLASS( const double mean,
                                    const double width,
                                    const double mB,
                                    const double mpsi,
                                    const double mkaon,
                                    const double mpion,
                                    const double radius,
                                    const double a,
                                    const double r,
                                    const double cutoff ) :
  KPiResonanceBase( mB, mpsi, mkaon, mpion ),
  mean_   ( mean ),
  width_  ( width ),
  radius_ ( radius ),
  a_      ( a ),
  r_      ( r ),
  cutoff_ ( cutoff ) {} ;


KPiResonanceLASS::KPiResonanceLASS( const double mB,
                                    const double mpsi,
                                    const double mkaon,
                                    const double mpion ) : 
  KPiResonanceBase( mB, mpsi, mkaon, mpion ), 
  mean_   ( 1.425 ),
  width_  ( 0.270 ),
  radius_ ( 1.600 ),
  a_      ( 3.830 ),
  r_      ( 2.860 ),
  cutoff_ ( mB ){} ; 

KPiResonanceLASS::KPiResonanceLASS( const KPiResonanceLASS& other ) :
  KPiResonanceBase( other ),
  mean_   ( other.mean_ ),
  width_  ( other.width_ ),
  radius_ ( other.radius_ ),
  a_      ( other.a_ ),
  r_      ( other.r_ ),
  cutoff_ ( other.cutoff_ ) {} ;



std::complex<double> KPiResonanceLASS::lineshape( const double m, const double pB, const double pR ) const {
  
  if ( forbidden( m ) ){
    return std::complex<double>(0,0);
  }

  const double pB0  = KinematicFunctions::momentum( mB_  , mpsi_ , mean_  );
  
  // K/pi momentum from resonance decay
  const double pR0  = KinematicFunctions::momentum( mean_, mkaon_, mpion_ ); 
  
  const double bfB  = KinematicFunctions::barrier( pB0, pB, radius_, 1 );
  
  const double gammaf = width_*(mean_/m)*(pR/pR0);
  
  const double norm   = bfB*( pB/pB0 );

    
  const double cotr   = ( mean_*mean_ - m*m )/( gammaf*mean_ );  

  const double cotb   = ( 1.0/( a_*pR ) ) + 0.5*r_*pR;
  
  const double deltab = atan( 1.0/cotb );
  
  std::complex<double> cotbnumer = std::complex<double>( norm,  0. );
  std::complex<double> cotbdenom = std::complex<double>( cotb, -1. );
  
  std::complex<double> cotrnumer = std::polar( norm , 2.0*deltab );
  std::complex<double> cotrdenom = std::complex<double>( cotr , -1.0 );
  
  if ( m > cutoff_ ){
    return (cotrnumer/cotrdenom);
  }

  return ( (cotbnumer/cotbdenom) + (cotrnumer/cotrdenom) );
}
