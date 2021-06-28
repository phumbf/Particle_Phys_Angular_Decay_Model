/*
 * Model for (psi pi) resonances
 * 
 */

#include "PsiPiResonance.h"
#include <cmath>

#include "Wignerd.h"
#include "KinematicFunctions.h" 

/* 
 * 
 * Note: Parity considerations enforce
 * 
 *      A_-lambda = - Pz * (-1)^Jz A_lambda
 * 
 * As a consequence, if Pz * (-1)^Jz is negative A0 = 0. 
 *
 * If A0 is non-zero it can be related to Ap and Am 
 * using LS amplitudes. 
 *
 */
 
PsiPiResonance::PsiPiResonance( const PsiPiResonance& other ) :
    mean_   ( other.mean_ ),
    width_  ( other.width_ ),
    mB_   ( other.mB_ ),
    mpsi_ ( other.mpsi_ ),
    mkaon_( other.mkaon_ ),
    mpion_( other.mpion_ ),
    LR_   ( other.LR_ ) ,
    LB_   ( other.LB_ ) ,
    radiusR_ ( other.radiusR_ ) ,
    radiusB_ ( other.radiusB_ ) ,
    min_  ( other.min_ ),
    max_  ( other.max_ ),
    A0_   ( other.A0_ ),
    Ap_   ( other.Ap_ ),
    Am_   ( other.Am_ ){} ; 

PsiPiResonance::PsiPiResonance( const double mean,
                                const double width,
                                const double mB,
                                const double mpsi,
                                const double mkaon,
                                const double mpion,
                                const unsigned int LR ,
                                const unsigned int LB ,
                                const double radiusR,
                                const double radiusB ) :
    mean_   ( mean ),
    width_  ( width ),
    mB_   ( mB ),
    mpsi_ ( mpsi ),
    mkaon_( mkaon ),
    mpion_( mpion ),
    LR_   ( LR ) ,
    LB_   ( LB ) ,
    radiusR_ ( radiusR ) ,
    radiusB_ ( radiusB ) ,
    min_  ( mpsi + mpion ),
    max_  ( mB - mkaon ),
    A0_   ( 0, 0 ),
    Ap_   ( 0, 0 ),
    Am_   ( 0, 0 ){} ;





void PsiPiResonance::setAmplitudesPolar( const double magA0, const double magAp, const double magAm,
					 const double phaseA0, const double phaseAp, const double phaseAm ){
    
  A0_ = std::polar( magA0, phaseA0 );
  
  if ( LR_ > 0 ){ 
    Ap_ = std::polar( magAp, phaseAp );
    Am_ = std::polar( magAm, phaseAm ); 
  }
  
  return ;
}


void PsiPiResonance::setAmplitudesCartesian( const double reA0, const double imA0, 
					     const double reAp, const double imAp, 
					     const double reAm, const double imAm ){
  
  A0_ = std::complex<double>( reA0, imA0 );
  
  if ( LR_ > 0 ){ 
    Ap_ = std::complex<double>( reAp, imAp );
    Am_ = std::complex<double>( reAm, imAm );
  }

  return ; 
}


std::complex<double> PsiPiResonance::lineshape( const double m, const double pB, const double pR ) const {
    
  const double meff = KinematicFunctions::meffective( mean_, min_, max_ );
    
  const double pB0  = KinematicFunctions::momentum( mB_  , mkaon_ , meff  );
  
  const double pR0  = KinematicFunctions::momentum( mean_, mpsi_  , mpion_ );
  
  const double bfR = KinematicFunctions::barrier( pR0, pR, radiusR_, LR_ );
  
  const double bfB = KinematicFunctions::barrier( pB0, pB, radiusB_, LB_ );
  
  const double gammaf  = width_*(mean_/m)*bfR*bfR*std::pow( pR/pR0, 2*LR_ + 1 );
  
  const double orbital = std::pow( pB/mB_ , LB_ )*std::pow( pR/m , LR_ );
  
  std::complex<double> fn( bfB*bfR*orbital, 0 );
  
  fn /= std::complex<double>( mean_*mean_ - m*m, -mean_*gammaf );
  
  return fn;
}



std::complex<double> PsiPiResonance::evaluate( const Event& event, const int lambda ) const {
    
  static std::complex<double> i(0,1);
  
  unsigned int j = 1;
  
  std::complex<double> result(0,0);
  
  result += ( A0_*
	      Wignerd::djmpm( event.cospsiz_, 1, lambda, 0 )*
	      Wignerd::djmpm( event.cosz_   , j, 0, 0 ) );
  
  result += ( Ap_*
	      Wignerd::djmpm( event.cospsiz_, 1, lambda, 1 )*
	      Wignerd::djmpm( event.cosz_   , j, 0, 1 )*
	      std::exp( -i*event.phiz_ ) );
  
  result += ( Am_*
	      Wignerd::djmpm( event.cospsiz_, 1, lambda, -1 )*
	      Wignerd::djmpm( event.cosz_   , j, 0, -1 )*
	      std::exp( i*event.phiz_ ) );
  
  result *= lineshape( event.mz_, event.pBz_, event.pZ_ ) ;
  
  result *= std::exp( -i*(lambda*event.deltaphi_) ) ;
  
  return result;
}


std::string PsiPiResonance::name() const {
    return name_ ;
}

bool PsiPiResonance::isNamed() const {
    return ( name_.length() > 0 );
}

bool PsiPiResonance::isNamed( const std::string& name ) const {
    return ( name_ == name );
}

bool PsiPiResonance::isNameIn( const std::string& name ) const { 
  return ( name.find( name_ ) != std::string::npos );
}

void PsiPiResonance::setName( const std::string& name ){
    name_ = name;
}


