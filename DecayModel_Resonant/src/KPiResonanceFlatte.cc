/* Model for the lineshape using a Flatte shape */

#include "KPiResonanceFlatte.h"
#include "HelperFunctions.h" 
#include "KinematicFunctions.h" 
#include "Variables.h"
#include "Event.h"

KPiResonanceFlatte::KPiResonanceFlatte( const double mean,
					const double gkk,
					const double gpipi,
					const double alpha,
					const double mB,
					const double mpsi,
					const double mkaon,
					const double mpion,
					const unsigned int LB,
					const double radius ):
  KPiResonanceBase(mB,mpsi,mkaon,mpion),
  mean_   ( mean ),
  gkk_    ( gkk  ),
  gpipi_  ( gpipi ),
  alpha_  ( alpha ),
  radius_ ( radius ),
  LB_     ( LB ){} ;


KPiResonanceFlatte::KPiResonanceFlatte(const KPiResonanceFlatte& other) :
  KPiResonanceBase( other ),
  mean_  ( other.mean_ ),
  gkk_   ( other.gkk_ ),
  gpipi_ ( other.gpipi_ ),
  alpha_ ( other.alpha_ ),
  radius_( other.radius_ ),
  LB_    ( other.LB_ ){};

std::complex<double> KPiResonanceFlatte::lineshape( const double m, 
						    const double pB, 
						    const double /*pR*/ ) const {

  if ( forbidden( m ) ){
    return std::complex<double>(0,0);
  }
  
  //Phase space factors
  const double s = m*m;
  
  const double mKpsq  = std::pow( Variables::mKaon   , 2 ); 
  const double mKzsq  = std::pow( Variables::mKZero  , 2 );
  const double mpipsq = std::pow( Variables::mPion   , 2 );
  const double mpizsq = std::pow( Variables::mPiZero , 2 );
  
  double gamma = 0;
  double delta = 0;

  gamma += (1./3.)*gpipi_*std::sqrt( 1. - 4.*mpizsq/s ); 
  gamma += (2./3.)*gpipi_*std::sqrt( 1. - 4.*mpipsq/s ); 
  
  if ( s > 4.*mKpsq ){ 
    double k = 
      KinematicFunctions::momentum( m, Variables::mKaon, Variables::mKaon );

    double rhokk = std::sqrt( 1. - 4.*mKpsq/s );
    
    gamma += (1./2.)*gkk_*rhokk*std::exp( -2.*alpha_*k*k );
  }
  else { 
    delta += (1./2.)*gkk_*std::sqrt( 4.*mKpsq/s - 1. );
  }
  
  if ( s > 4.*mKzsq ){ 
    double k = 
      KinematicFunctions::momentum( m, Variables::mKZero, Variables::mKZero );
    
    double rhokk = std::sqrt( 1. - 4.*mKpsq/s );
      
    gamma += (1./2.)*gkk_*rhokk*std::exp( -2.*alpha_*k*k );
  }
  else { 
    delta += (1./2.)*gkk_*std::sqrt( 4.*mKzsq/s - 1. );
  }
  
  const double pB0 = KinematicFunctions::momentum( mB_  , mpsi_ , mean_ );
  
  const double bfB = KinematicFunctions::barrier( pB0, pB, radius_, LB_ );
  
  const double orbital = std::pow( pB/mB_ , LB_ );

  std::complex<double> numer( bfB*orbital, 0 );
  std::complex<double> denom( mean_*mean_ - s + delta, -mean_*gamma );
 
  return (numer / denom);
}
