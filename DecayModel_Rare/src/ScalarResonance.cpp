// Local
#include "ScalarResonance.h"
#include "Parameters.h"
#include "Kinematic.h"
#include "EvtGenWilson.h"

// std
#include <cmath>
#include <complex>

//Class for a scalar (i.e. spin = 0) resonance 
ScalarResonanceBase::ScalarResonanceBase(){}

ScalarResonanceBase::~ScalarResonanceBase(){}

ScalarFormFactorBase* ScalarResonanceBase::getFFs() const {
  
  return m_ff.get();
}

void ScalarResonanceBase::setFFs( std::unique_ptr<ScalarFormFactorBase>& ff ){
  m_ff = std::move(ff);
}

std::complex<double> ScalarResonanceBase::APERP( const double /*s*/ ,
						 const double /*m*/ ,
						 const double /*LR*/ ) const {
  return std::complex<double>(0,0);
}

std::complex<double> ScalarResonanceBase::APARA( const double /*s*/ ,
						 const double /*m*/ ,
						 const double /*LR*/ ) const {
  return std::complex<double>(0,0);
}

std::complex<double> ScalarResonanceBase::AZERO( const double s ,
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
  const double kinematic = std::sqrt(lambdaM) / std::sqrt(s);
  
  // Eq A1 arXiv:1111.1513
  std::complex<double> res =
    (C9eff - LR*C10)*m_ff->fplus(s) + 
    2.0*C7*Parameters::mb*m_ff->fT(s)/(m_mothermass + m);
  
  res *= kinematic;
  res *= normalisation( s, m );

  // non-factorisable correction
  res *= m_resfactor*prefactor(0,m,lambdaM,lambdaH);
  res *= ( 1. + m_cZERO );

  return res;
}

std::complex<double> ScalarResonanceBase::ATIME( const double s,
						 const double m ) const{
  
  if ( Kinematic::isForbidden( s, Parameters::ml, m, m_mone, m_mtwo, m_mothermass ) ){
    return std::complex<double>(0,0); 
  }
  
  const double lambdaM    = Kinematic::lambda( s, m, m_mothermass );
  const double lambdaH    = Kinematic::lambda( m*m, m_mone, m_mtwo );

  std::complex<double> C10   = m_WC(10) - m_WP(10);
  if(m_isEvtGenWC){
    C10 =   EvtGenWilson::GetC10Eff(s,true,m_isbtod);
  }
  
  const double kinematic = ( m_mothermass*m_mothermass - m*m )/std::sqrt(s);
  
  // Eq A1 arXiv:1111.1513
  std::complex<double> res = 2.0 * C10 * m_ff->fzero(s); 
  
  res *= kinematic;
  res *= m_resfactor*prefactor(0,m,lambdaM,lambdaH);
  res *= normalisation( s, m );
  
  return res;
}
