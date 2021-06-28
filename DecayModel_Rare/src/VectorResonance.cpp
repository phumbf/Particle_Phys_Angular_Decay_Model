// Local
#include "VectorResonance.h"
#include "Kinematic.h"
#include "Parameters.h"
#include "EvtGenWilson.h"

// std
#include <cmath>
#include <complex>

#include "TCanvas.h"

//Class for a vector (i.e. spin > 0) resonance 
VectorResonanceBase::VectorResonanceBase(){}

VectorResonanceBase::~VectorResonanceBase(){}

VectorFormFactorBase* VectorResonanceBase::getFFs() const {
  return m_ff.get();
}

void VectorResonanceBase::setFFs( std::unique_ptr<VectorFormFactorBase>& ff ){

  m_ff = std::move(ff);
}

std::complex<double> VectorResonanceBase::APARA( const double s ,
						 const double m ,
						 const double LR ) const {
  
  if ( Kinematic::isForbidden( s, Parameters::ml, m, m_mone, m_mtwo, m_mothermass ) ){
    return std::complex<double>(0,0); 
  }
  
  std::complex<double> C9eff = m_WC(9)  - m_WP(9) + m_WC.Y(s);
  std::complex<double> C10   = m_WC(10) - m_WP(10);
  std::complex<double> C7    = m_WC(7)  - m_WP(7);
  
  if(m_isEvtGenWC){
    C9eff = EvtGenWilson::GetC9Eff(s,true,m_isbtod);
    C10 =   EvtGenWilson::GetC10Eff(s,true,m_isbtod);
    C7  =   EvtGenWilson::GetC7Eff(s,true,m_isbtod);
  }
  
  const double lambdaM    = Kinematic::lambda( s, m, m_mothermass );
  const double lambdaH    = Kinematic::lambda( m*m, m_mone, m_mtwo );
  const double kinematic = std::sqrt(2.)*( m_mothermass*m_mothermass - m*m );
  
  // Eq A4 arXiv:1111.1513 and 11 in arXiv:1207.4004
  std::complex<double> res = 
    ( C9eff - LR*C10 )*m_ff->A1(s)/( m_mothermass - m ) + 
    2.*C7*Parameters::mb*m_ff->T2(s)/s;
    
  res *= kinematic; 
  res *= normalisation( s, m );
  
  // a.m. factor
  res *= m_resfactor*prefactor( 1, m, lambdaM, lambdaH );
  
  // non-factorisable correction
  res *= ( 1. + m_cPARA );

  return res;
}

std::complex<double> VectorResonanceBase::APERP( const double s ,
						 const double m ,
						 const double LR ) const {
  
 
  if ( Kinematic::isForbidden( s, Parameters::ml, m, m_mone, m_mtwo, m_mothermass ) ){
    return std::complex<double>(0,0); 
  }
  
  std::complex<double> C9eff = m_WC(9)  + m_WP(9) + m_WC.Y(s);
  std::complex<double> C10   = m_WC(10) + m_WP(10);
  std::complex<double> C7    = m_WC(7)  + m_WP(7);
  
  if(m_isEvtGenWC){
    C9eff = EvtGenWilson::GetC9Eff(s,true,m_isbtod);
    C10 =   EvtGenWilson::GetC10Eff(s,true,m_isbtod);
    C7  =   EvtGenWilson::GetC7Eff(s,true,m_isbtod);
  }
  
  const double lambdaM    = Kinematic::lambda( s, m, m_mothermass );
  const double lambdaH    = Kinematic::lambda( m*m, m_mone, m_mtwo );
  const double kinematic = std::sqrt( 2.*lambdaM );
 
  // Eq A4 arXiv:1111.1513 and 11 in arXiv:1207.4004
  std::complex<double> res = 
    ( C9eff - (LR*C10) )*m_ff->V(s)/( m_mothermass + m ) + 
    2.0*C7*Parameters::mb*m_ff->T1(s)/s; 

  res *= kinematic;
  res *= normalisation( s , m );
  
  // a.m. factor
  res *= m_resfactor*prefactor( 1, m, lambdaM, lambdaH );

  // non-factorisable correction
  res *= ( 1. + m_cPERP );
  
  return -res;
}

std::complex<double> VectorResonanceBase::AZERO( const double s ,
						 const double m ,
						 const double LR ) const{

  
  if ( Kinematic::isForbidden( s, Parameters::ml, m, m_mone, m_mtwo, m_mothermass ) ){
    return std::complex<double>(0,0); 
  }
  
  std::complex<double> C9eff = m_WC(9)  - m_WP(9) + m_WC.Y(s);
  std::complex<double> C10   = m_WC(10) - m_WP(10);
  std::complex<double> C7    = m_WC(7)  - m_WP(7);

  if(m_isEvtGenWC){
    C9eff = EvtGenWilson::GetC9Eff(s,true,m_isbtod);
    C10 =   EvtGenWilson::GetC10Eff(s,true,m_isbtod);
    C7  =   EvtGenWilson::GetC7Eff(s,true,m_isbtod);
  }
  
  const double lambdaM    = Kinematic::lambda( s, m, m_mothermass );
  const double lambdaH    = Kinematic::lambda( m*m, m_mone, m_mtwo );
  const double kinematic = 2.*m*std::sqrt(s);
  
  const double mVsq = m*m;
  const double mBsq = m_mothermass*m_mothermass;
  

  // Eq 11 in arXiv:1207.4004
  std::complex<double> res = 
    ( C9eff - LR*C10 )*( ( mBsq - mVsq - s )*( m_mothermass + m )*m_ff->A1(s) - 
			 ( lambdaM*m_ff->A2(s)/(m_mothermass + m) ) ) + 
    2.0*Parameters::mb*C7*( ( mBsq + 3.0*mVsq - s )*m_ff->T2(s) - 
			    ( lambdaM*m_ff->T3(s) )/(mBsq - mVsq) );
  
  
  res *= normalisation( s, m );
  res /= kinematic;
  
  // a.m. factor
  res *= m_resfactor*prefactor( 0, m, lambdaM, lambdaH );

  // non-factorisable correction
  res *= ( 1. + m_cZERO );
  
  return res;
}

std::complex<double> VectorResonanceBase::ATIME( const double s,
						 const double m ) const {

  if ( Kinematic::isForbidden( s, Parameters::ml, m, m_mone, m_mtwo, m_mothermass ) ){
    return std::complex<double>(0,0); 
  }
  
  std::complex<double> C10   = m_WC(10) - m_WP(10);
  
  if(m_isEvtGenWC){
    C10 =   EvtGenWilson::GetC10Eff(s,true,m_isbtod);
  }
  const double lambdaM    = Kinematic::lambda( s, m, m_mothermass );
  const double lambdaH    = Kinematic::lambda( m*m, m_mone, m_mtwo );
  const double kinematic = std::sqrt(lambdaM / s);
  
  std::complex<double> res = 2.0 * C10 * m_ff->A0(s);

  res *= kinematic; 
  res *= normalisation( s, m );
  
   // a.m. factor
  res *= m_resfactor*prefactor( 0, m, lambdaM, lambdaH );

  return res;
}

