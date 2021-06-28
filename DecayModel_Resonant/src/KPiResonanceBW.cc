/*
 * Model for the lineshape using a relativistic BW. 
 * 
 * The model inherits most of its properties 
 * from KPiResonaceBase.
 */

// Local
#include "KPiResonanceBW.h"
#include "HelperFunctions.h" 
#include "KinematicFunctions.h" 


KPiResonanceBW::KPiResonanceBW( const double mean,
                                const double width,
                                const double mB,
                                const double mpsi,
                                const double mkaon,
                                const double mpion,
                                const unsigned int LR,
                                const unsigned int LB,
                                const double radiusR,
                                const double radiusB  ) :
    KPiResonanceBase( mB, mpsi, mkaon, mpion ),
    mean_    ( mean ),
    width_   ( width ),
    LR_      ( LR ),
    LB_      ( LB ),
    radiusR_ ( radiusR ),
    radiusB_ ( radiusB ){} ;

KPiResonanceBW::KPiResonanceBW( const double mean,
                                const double width,
                                const double mB,
                                const double mpsi,
                                const double mkaon,
                                const double mpion,
                                const unsigned int LR,
                                const unsigned int LB,
                                const double radius) :
    KPiResonanceBase( mB, mpsi, mkaon, mpion ),
    mean_    ( mean ),
    width_   ( width ),
    LR_      ( LR ),
    LB_      ( LB ),
    radiusR_ ( radius ),
    radiusB_ ( radius ){} ;

KPiResonanceBW::KPiResonanceBW( const KPiResonanceBW& other ) :
    KPiResonanceBase( other ),
    mean_    ( other.mean_ ),
    width_   ( other.width_ ),
    LR_      ( other.LR_ ),
    LB_      ( other.LB_ ),
    radiusR_ ( other.radiusR_ ),
    radiusB_ ( other.radiusB_ ){} ;

std::complex<double> KPiResonanceBW::lineshape( const double m, const double pB, const double pR ) const {
  if ( forbidden( m ) ){
    return std::complex<double>(0,0);
  }
  const double meff = KinematicFunctions::meffective( mean_, min_, max_ );
  const double pB0 = KinematicFunctions::momentum( mB_  , mpsi_ , meff  );
  
  // K/pi momentum from resonance decay
  const double pR0 = KinematicFunctions::momentum( mean_, mkaon_, mpion_ ); 
  
  const double bfR = KinematicFunctions::barrier( pR0, pR, radiusR_, LR_ );
  
  const double bfB = KinematicFunctions::barrier( pB0, pB, radiusB_, LB_ );
 
  const double gammaf  = width_*(mean_/m)*bfR*bfR*std::pow( pR/pR0, 2*LR_ + 1 );
  
  const double orbital = std::pow( pB/mB_ , LB_ )*std::pow( pR/m , LR_ );
  std::complex<double> numer( bfB*bfR*orbital, 0 );

  std::complex<double> denom( mean_*mean_ - m*m, -mean_*gammaf );
  return numer/denom;
}
