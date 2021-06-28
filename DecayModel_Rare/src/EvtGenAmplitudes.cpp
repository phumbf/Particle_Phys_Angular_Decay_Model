#include "TMath.h"

#include <iostream>
#include <complex> 
#include <cmath>
#include <limits>

#include "EvtGenAmplitudes.h"
#include "EvtGenFormFactor.h" 
#include "EvtGenWilson.h" 
#include "Parameters.h"
#include "Kinematic.h" 

#include "Math/Functor.h"
#include "Math/Integrator.h"

//Amplitudes class to mirror the treatment found in the EvtGen 
//software package https://evtgen.hepforge.org/

EvtGenAmplitudes::EvtGenAmplitudes( const bool isBs, const bool isbtod ) :
  btod_( isbtod ),
  C9_  ( 0, 0 ) ,
  C10_ ( 0, 0 ) 
{
  if ( isBs ) {
    FF_ = new EvtGenBsFormFactor;
  } 
  else { 
    FF_ = new EvtGenBzFormFactor ; 
  }
} 


EvtGenAmplitudes::EvtGenAmplitudes( const EvtGenAmplitudes& other ) : 
  btod_( other.btod_ ),
  C9_  ( other.C9_ ),
  C10_ ( other.C10_ ),
  FF_  ( other.FF_->clone() ){} 


EvtGenAmplitudes::~EvtGenAmplitudes() {
  if ( FF_ ) delete FF_;
} 

std::complex<double> EvtGenAmplitudes::APERP( const double s, const double LR, const int flavour ) const
{
  if ( Kinematic::isForbidden( s, Parameters::ml, 
			       Parameters::mV, Parameters::mKaon, Parameters::mPion, 
			       Parameters::mB_Bz ) ){
    return std::complex<double>(0,0); 
  }
  
  
  std::complex<double> C9eff = EvtGenWilson::GetC9Eff ( s, true, btod_ ) + C9_; 
  
  std::complex<double> C10   = EvtGenWilson::GetC10Eff( s, true, btod_ ) + C10_;
  
  std::complex<double> C7    = EvtGenWilson::GetC7Eff ( s, true, btod_ );
  
  std::complex<double> result;
  
  result = ( ( C9eff - LR*C10 )*FF_->V(s)/( Parameters::mB_Bz + Parameters::mV ) +
	     ( 2.*Parameters::mb*C7*FF_->T1(s)/s ) );
  
  if ( flavour < 0 ) {
    result = -std::conj( result );
  }

  const double lambda = 
    Kinematic::lambda( s, Parameters::mV, Parameters::mB_Bz ); 

  return std::sqrt( 2.*lambda )*norm( s )*result;
}


std::complex<double> EvtGenAmplitudes::APARA( const double s, const double LR, const int flavour ) const
{
  if ( Kinematic::isForbidden( s, Parameters::ml, 
			       Parameters::mV, Parameters::mKaon, Parameters::mPion, 
			       Parameters::mB_Bz ) ){
    return std::complex<double>(0,0); 
  }    
  
  std::complex<double> C9eff = EvtGenWilson::GetC9Eff ( s, true, btod_ ) + C9_; 
  
  std::complex<double> C10   = EvtGenWilson::GetC10Eff( s, true, btod_ ) + C10_;
  
  std::complex<double> C7    = EvtGenWilson::GetC7Eff ( s, true, btod_ );
  
  std::complex<double> result;
  
  result = ( (C9eff-LR*C10)*FF_->A1(s)/( (Parameters::mB_Bz) - (Parameters::mV) ) +
	     (2*(Parameters::mb)*C7*FF_->T2(s)/s ) );
  
  if ( flavour < 0 ) {
    result = std::conj( result );
  }
  
  const double mBsq = std::pow( Parameters::mB_Bz, 2 );
  const double mVsq = std::pow( Parameters::mV, 2 );
  
  return -std::sqrt(2.)*( mBsq - mVsq )*norm( s )*result;
}

std::complex<double> EvtGenAmplitudes::AZERO( const double s, const double LR, const int flavour ) const
{
  if ( Kinematic::isForbidden( s, Parameters::ml, 
			       Parameters::mV, Parameters::mKaon, Parameters::mPion, 
			       Parameters::mB_Bz ) ){
    return std::complex<double>(0,0); 
  }      
    
   std::complex<double> C9eff = EvtGenWilson::GetC9Eff ( s, true, btod_ ) + C9_ + C10_;  
   
   std::complex<double> C10   = EvtGenWilson::GetC10Eff( s, true, btod_ );
   
   std::complex<double> C7    = EvtGenWilson::GetC7Eff ( s, true, btod_ );
   
   std::complex<double> result;
   
   const double mBsq = std::pow( Parameters::mB_Bz, 2 );
   const double mVsq = std::pow( Parameters::mV, 2 );
   const double mSum = Parameters::mB_Bz + Parameters::mV;

   const double lambda = 
     Kinematic::lambda( s, Parameters::mV, Parameters::mB_Bz ); 
   
   double c9c10  = ( (mBsq - mVsq -s)*mSum*FF_->A1(s) - lambda*FF_->A2(s)/mSum );
   
   double photon = ( ( mBsq +3.*mVsq - s )*FF_->T2(s) - lambda*FF_->T3(s)/(mBsq  - mVsq ) );
   
   result = ( (C9eff-LR*C10)*c9c10 + 2.*(Parameters::mb)*C7*photon );
   
   if ( flavour < 0 ) {
     result = std::conj( result );
   }
   
   return -(0.5*norm(s)/((Parameters::mV)*std::sqrt(s)))*result;
}


std::complex<double> EvtGenAmplitudes::ATIME( const double s, const int flavour ) const { 
  
  if ( Kinematic::isForbidden( s, Parameters::ml, 
			       Parameters::mV, Parameters::mKaon, Parameters::mPion, 
			       Parameters::mB_Bz ) ){
    return std::complex<double>(0,0); 
  }      
  
  const double lambda = 
    Kinematic::lambda( s, Parameters::mV, Parameters::mB_Bz ); 
  
  std::complex<double> C10   = EvtGenWilson::GetC10Eff( s, true, btod_ ) + C10_;
  std::complex<double> val   = 2.*C10*FF_->A0(s);
  
  if ( flavour < 0 ){ 
    val = std::conj( val );
  }
  
  return norm(s)*val*std::sqrt(lambda/s) ;
}


double EvtGenAmplitudes::norm( const double s ) const
{
  const double beta   = Kinematic::beta( s, Parameters::ml ); 
  const double lambda = Kinematic::lambda(s, Parameters::mV, Parameters::mB_Bz);
  
  double Nsq  = s*std::sqrt( lambda )*beta;
    
  return std::sqrt( Nsq );
}


double EvtGenAmplitudes::rate( const double s ) const
{
  double result = 0;
    
  result += std::norm( AZERO( s,  1. ) );
  result += std::norm( AZERO( s, -1. ) );
  result += std::norm( APERP( s,  1. ) );
  result += std::norm( APERP( s, -1. ) );
  result += std::norm( APARA( s,  1. ) );
  result += std::norm( APARA( s, -1. ) );
  
  return result;
}


double EvtGenAmplitudes::I1s( const AmplitudeSet& amps, const double s ) const { 
  
  double result = 0;
  
  const double betasq = 
    Kinematic::betasq( s, Parameters::ml );
  
  result += ( (2. + betasq)/4. )*( std::norm(amps.APERPL) + std::norm(amps.APARAL) + 
					std::norm(amps.APERPR) + std::norm(amps.APARAR) );
  
  result += ( 1. - betasq )*std::real( amps.APERPL*std::conj(amps.APERPR) + 
					    amps.APARAL*std::conj(amps.APARAR) );

  return result; 
}

double EvtGenAmplitudes::I1c( const AmplitudeSet& amps, const double s  ) const {
  
  double result = 0;
  
  const double betasq = 
    Kinematic::betasq( s, Parameters::ml );
  
  result += ( std::norm(amps.AZEROL) + std::norm(amps.AZEROR) );
  
  result += ( 1. - betasq )*2.*std::real(amps.AZEROL*std::conj(amps.AZEROR) );
  
  return result; 
}
	
double EvtGenAmplitudes::I2s( const AmplitudeSet& amps, const double s  ) const {   

  const double betasq = 
     Kinematic::betasq( s, Parameters::ml );
  
  double result = (betasq/4.)*( std::norm(amps.APERPL) + std::norm(amps.APARAL) + 
				std::norm(amps.APERPR) + std::norm(amps.APARAR) );
  
  return result; 
}

double EvtGenAmplitudes::I2c( const AmplitudeSet& amps, const double s  ) const {

  const double betasq = 
     Kinematic::betasq( s, Parameters::ml );
  
  double result = -betasq*( std::norm(amps.AZEROL) + std::norm(amps.AZEROR) );
  
  return result; 
}

double EvtGenAmplitudes::I3( const AmplitudeSet& amps, const double s  ) const {
  
  const double betasq = 
    Kinematic::betasq( s, Parameters::ml );
  
  double result = (betasq/2.)*( std::norm(amps.APERPL) + std::norm(amps.APERPR) -
				std::norm(amps.APARAL) - std::norm(amps.APARAR) );
  
  return result; 
}


double EvtGenAmplitudes::I4( const AmplitudeSet& amps, const double s ) const {    
  
  const double betasq = 
    Kinematic::betasq( s, Parameters::ml );
  
  double result = (betasq/std::sqrt(2.))*( std::real( amps.AZEROL*std::conj(amps.APARAL) + 
						      amps.AZEROR*std::conj(amps.APARAR) ) );
  
  return result; 
}

double EvtGenAmplitudes::I5( const AmplitudeSet& amps, const double s  ) const {  
  
  const double betasq = 
    Kinematic::betasq( s, Parameters::ml );
  
  double result = std::sqrt(2.*betasq)*( std::real( amps.AZEROL*std::conj(amps.APERPL) - 
						    amps.AZEROR*std::conj(amps.APERPR) ) );
  
  return result; 
}
		
double EvtGenAmplitudes::I6s( const AmplitudeSet& amps, const double s  ) const {
  
  const double betasq = 
    Kinematic::betasq( s, Parameters::ml );
  
  double result = 2.*std::sqrt(betasq)*( std::real( amps.APARAL*std::conj(amps.APERPL) - 
						    amps.APARAR*std::conj(amps.APERPR) ) );
  
  return result; 
}


double EvtGenAmplitudes::I7( const AmplitudeSet& amps, const double s  ) const  { 
  
  const double betasq = 
    Kinematic::betasq( s, Parameters::ml );
  
  double result = std::sqrt(2.*betasq)*( std::imag( amps.AZEROL*std::conj(amps.APARAL) - 
						    amps.AZEROR*std::conj(amps.APARAR) ) );
  
  return result; 
}


double EvtGenAmplitudes::I8( const AmplitudeSet& amps, const double s  ) const {
  
  const double betasq = 
    Kinematic::betasq( s, Parameters::ml );
  
  double result = (betasq/std::sqrt(2.))*( std::imag( amps.AZEROL*std::conj(amps.APERPL) + 
						      amps.AZEROR*std::conj(amps.APERPR) ) );
  
  return result; 
}

double EvtGenAmplitudes::I9( const AmplitudeSet& amps, const double s  ) const {
  
  const double betasq = 
    Kinematic::betasq( s, Parameters::ml );

  double result = betasq*( std::imag( std::conj(amps.APARAL)*amps.APERPL + 
				      std::conj(amps.APARAR)*amps.APERPR ) );
  
  return result; 
}



double EvtGenAmplitudes::differential( const AmplitudeSet& amps, 
				       const double s  ,
				       const double ctl, 
				       const double ctk, 
				       const double phi ) const { 
  
  double result = ( I1s( amps, s )*(1. - ctk*ctk ) +
		    I1c( amps, s )*(ctk*ctk) + 
		    I2s( amps, s )*(1. - ctk*ctk )*(2.*ctl*ctl - 1.) +
		    I2c( amps, s )*(ctk*ctk)*(2.*ctl*ctl - 1.) + 
		    I3 ( amps, s )*(1. - ctk*ctk)*(1. - ctl*ctl)*cos(2.*phi) + 
		    I4 ( amps, s )*(2.*ctk*sqrt(1.-ctk*ctk))*(2.*ctl*sqrt(1.-ctl*ctl))*cos(phi) + 
		    I5 ( amps, s )*(2.*ctk*sqrt(1.-ctk*ctk))*sqrt(1. - ctl*ctl)*cos(phi) + 
		    I6s( amps, s )*(1. - ctk*ctk)*ctl + 
		    I7 ( amps, s )*(2.*ctk*sqrt(1.-ctk*ctk))*sqrt(1. - ctl*ctl)*sin(phi) + 
		    I8 ( amps, s )*(2.*ctk*sqrt(1.-ctk*ctk))*(2.*ctl*sqrt(1.-ctl*ctl))*sin(phi) + 
		    I9 ( amps, s )*(1. - ctk*ctk)*(1. - ctl*ctl)*sin(2.*phi) ); 
 
  
  return ( 9.*result )/( 32.*TMath::Pi() ) ; 
}

double EvtGenAmplitudes::operator()( const double* x ) {
  
  const double s = x[0]; 

  AmplitudeSet amps;
  
  amps.AZEROL = AZERO(s, 1, 1); 
  amps.AZEROR = AZERO(s,-1, 1); 
  amps.APARAL = APARA(s, 1, 1); 
  amps.APARAR = APARA(s,-1, 1); 
  amps.APERPL = APERP(s, 1, 1); 
  amps.APERPR = APERP(s,-1, 1); 
  
  double result =  differential( amps, s, x[1], x[2], x[3] );

  return result;
}

double EvtGenAmplitudes::operator()( const std::vector<double> &x ) {
  
  const double s = x[0]; 

  AmplitudeSet amps;
  
  amps.AZEROL = AZERO(s, 1, 1); 
  amps.AZEROR = AZERO(s,-1, 1); 
  amps.APARAL = APARA(s, 1, 1); 
  amps.APARAR = APARA(s,-1, 1); 
  amps.APERPL = APERP(s, 1, 1); 
  amps.APERPR = APERP(s,-1, 1); 
  
  double result =  differential( amps, s, x[1], x[2], x[3] );

  return result;
}

double EvtGenAmplitudes::differential( const double s, const double ctl, const double ctk, const double phi ) const {
  
  AmplitudeSet amps;
  
  amps.AZEROL = AZERO(s, 1, 1); 
  amps.AZEROR = AZERO(s,-1, 1); 
  amps.APARAL = APARA(s, 1, 1); 
  amps.APARAR = APARA(s,-1, 1); 
  amps.APERPL = APERP(s, 1, 1); 
  amps.APERPR = APERP(s,-1, 1); 

  double result =  differential( amps, s, ctl, ctk, phi );

  return result;
}



double EvtGenAmplitudes::integral( const double min, const double max ) const
{
  ROOT::Math::Functor1D  func( this, &EvtGenAmplitudes::rate ); //
  
  ROOT::Math::Integrator ig;
  
  ig.SetFunction(func);

  ig.SetRelTolerance(0.001);

  double result = ig.Integral( min, max );
  
  return result ;
}


