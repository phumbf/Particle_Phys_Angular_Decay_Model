#include "TMath.h"
#include "TRandom3.h"

#include <iostream>
#include <complex> 
#include <cmath>
#include <limits>
#include <functional>

#include "Amplitudes.h"
#include "Parameters.h"
#include "Kinematic.h" 
#include "SphericalHarmonics.h"

#include "Math/Functor.h"
#include "Math/Integrator.h"


#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/Functor.h"

//Amplitude class built using the full amplitude model 

Amplitudes::Amplitudes() { } 

Amplitudes::~Amplitudes() {} 

void Amplitudes::addResonance( std::unique_ptr<ResonanceBase>& r ){ 
  m_res.push_back( std::move(r) );
}

ResonanceBase* Amplitudes::getResonance( const std::string name ) const { 
   for ( const auto & r: m_res ){
     if ( r->hasName( name ) ) {
       return r.get();
     }
   }
   return nullptr; 
}


double Amplitudes::FL( const ResonanceBase* res, const double s, const double mhh ) const { 
  /*
   * Compute the observable FL for 
   * given ResonanceBase at point in q^2 and mhh.
   */
  
  double result = 0;
  double norm   = 0; 
  
  AmplitudeSet amp;
  res->getAmplitudes( s, mhh, amp);
 

  result += std::norm( amp.AZEROL );
  result += std::norm( amp.AZEROR );

  norm += result;
  norm += std::norm( amp.APARAL );
  norm += std::norm( amp.APARAR );
  norm += std::norm( amp.APERPL );
  norm += std::norm( amp.APERPR );

  if ( norm > 0 ){ 
    result /= norm; 
  }
  
  return result;
}

double Amplitudes::AFB( const ResonanceBase* res, const double s, const double mhh) const {

	/* Compute forward-backwards 
	 * asymmetry of dilepton system
	 */

	double result = 0;
	double norm = 0;

  	AmplitudeSet amp;
  
	res->getAmplitudes( s, mhh, amp );

	result += std::real(  amp.APARAL*std::conj(amp.APERPL) 
		            - amp.APARAR*std::conj(amp.APERPR));

	norm += std::norm( amp.AZEROL );
	norm += std::norm( amp.AZEROR );
	norm += std::norm( amp.APARAL );
	norm += std::norm( amp.APARAR );
	norm += std::norm( amp.APERPL );
	norm += std::norm( amp.APERPR );

	if ( norm > 0 ){ 
		result /= norm; 
	}

	return (4./3.)*result;
}

double Amplitudes::AT2( const ResonanceBase* res, const double s, const double mhh ) const { 
  /*
   * Compute the observable AT2 for 
   * given ResonanceBase at point in q^2 and mhh.
   */
  
  double result = 0;
  double norm   = 0; 
  
  AmplitudeSet amp;
  
  res->getAmplitudes( s, mhh, amp );
  
  result += std::norm( amp.APERPL );
  result += std::norm( amp.APERPR );
  result -= std::norm( amp.APARAL );
  result -= std::norm( amp.APARAR );
  
  norm += std::norm( amp.APARAL );
  norm += std::norm( amp.APARAR );
  norm += std::norm( amp.APERPL );
  norm += std::norm( amp.APERPR );

  if(norm > 0){
	  result/=norm;
  }
	  return result;
}

double Amplitudes::rate( const ResonanceBase* res, const double s, const double mhh ) const {
  /*
   * Compute the differential decay rate for 
   * given ResonanceBase at point in q^2 and mhh.
   */
  
  AmplitudeSet amp; 
  AmplitudeSet ampCP; 
  
  res->getAmplitudes( s, mhh, amp );
  
  double result = 0;
 
  result += std::norm( amp.AZEROL );
  result += std::norm( amp.AZEROR );
  result += std::norm( amp.APARAL );
  result += std::norm( amp.APARAR );
  result += std::norm( amp.APERPL );
  result += std::norm( amp.APERPR );

  return result;
}


double Amplitudes::rate( const double s, const double mhh ) const { 
  
  /*
   * Compute the total differential decay rate at 
   * given point in q^2 and mhh.
   */ 
  
  AmplitudeSet amp; 
  
  std::vector< std::complex<double> > a_zero_L; 
  std::vector< std::complex<double> > a_zero_R; 
  std::vector< std::complex<double> > a_perp_L; 
  std::vector< std::complex<double> > a_perp_R;
  std::vector< std::complex<double> > a_para_L;
  std::vector< std::complex<double> > a_para_R;
  
  int lmax = -1;
  
  for ( const auto & r: m_res ){
    if ( r->getL() > lmax ) lmax = r->getL();
  }
  
  for ( int l = 0; l <= lmax; l++ ){ 
    a_zero_L.emplace_back( 0, 0 );
    a_zero_R.emplace_back( 0, 0 );
    a_para_L.emplace_back( 0, 0 );
    a_para_R.emplace_back( 0, 0 );
    a_perp_L.emplace_back( 0, 0 );
    a_perp_R.emplace_back( 0, 0 );
  }

  for ( const auto & r: m_res ){
    r->getAmplitudes( s, mhh, amp );

    a_zero_L[ r->getL() ] += amp.AZEROL;
    a_zero_R[ r->getL() ] += amp.AZEROR;
    a_para_L[ r->getL() ] += amp.APARAL;
    a_para_R[ r->getL() ] += amp.APARAR;
    a_perp_L[ r->getL() ] += amp.APERPL;
    a_perp_R[ r->getL() ] += amp.APERPR;

  }

  double result = 0; 

  for ( int l = 0; l <= lmax; l++ ){ 
    result += std::norm( a_zero_L[l] );
    result += std::norm( a_zero_R[l] );
    result += std::norm( a_para_L[l] );
    result += std::norm( a_para_R[l] );
    result += std::norm( a_perp_L[l] );
    result += std::norm( a_perp_R[l] );
  }
  
  return result; 
}


double Amplitudes::rate_fn( const ResonanceBase* res, const double* x ) const { 
  /*
   * Wrapping function for integration of
   * decay rate for a single resonance. 
   */ 
  
  return rate( res, x[0], x[1] );
}

double Amplitudes::total_rate_fn( const double* x ) const { 
  /*
   * Wrapping function for integration of
   * total decay rate. 
   */ 
  return rate( x[0], x[1] );
}

double Amplitudes::FL( const std::string name, const double s, const double mhh ) const {
  /*
   * Compute the observable FL for the named
   * resonance at point in q^2 and mhh.
   */
  for ( const auto & r: m_res ){
    if ( r->hasName( name ) ) return FL( r.get(), s, mhh ); 
  }
  
  return 0;
}

double Amplitudes::AFB( const std::string name, const double s, const double mhh ) const {
  /*
   * Compute the observable AFB for the named
   * resonance at point in q^2 and mhh.
   */
  for ( const auto & r: m_res ){
    if ( r->hasName( name ) ) return AFB( r.get(), s, mhh ); 
  }
  
  return 0;
}
double Amplitudes::AT2( const std::string name, const double s, const double mhh ) const {
  /*
   * Compute the observable AT2 for the named
   * resonance at point in q^2 and mhh.
   */
  for ( const auto & r: m_res ){
    if ( r->hasName( name ) ) return AT2( r.get(), s, mhh ); 
  }
  
  return 0;
}

double Amplitudes::rate( const std::string name, const double s, const double mhh ) const {
   /*
    * Compute the differential decay rate for
    * the named resonance at point in q^2 and mhh.
    */
  
  for ( const auto & r: m_res ){
    if ( r->hasName( name ) ) return rate( r.get(), s, mhh ); 
  }
    
  return 0;
}



double Amplitudes::I1s(const double ctk,
		       const double mhh,
		       const double s ) const { 
  
  double result = 0;

  std::complex<double> a_perp_L(0,0);
  std::complex<double> a_perp_R(0,0);
  std::complex<double> a_para_L(0,0);
  std::complex<double> a_para_R(0,0);

  AmplitudeSet amp; 
  
  //Loop through the resonances
  for ( const auto & r: m_res ){

    //Retrieve the AmplitudeSet
    r->getAmplitudes( s, mhh, amp );
    
    const int l = r->getL();
    
    a_perp_L += amp.APERPL * SphericalHarmonics::Y(l,-1,ctk);
    a_perp_R += amp.APERPR * SphericalHarmonics::Y(l,-1,ctk);
    a_para_L += amp.APARAL * SphericalHarmonics::Y(l,-1,ctk);
    a_para_R += amp.APARAR * SphericalHarmonics::Y(l,-1,ctk);
  }

  const double betasq = 
    Kinematic::betasq( s, Parameters::ml );
	  
  result += ( (2. + betasq)/4. )*( std::norm(a_perp_L) + 
				   std::norm(a_para_L) + 
				   std::norm(a_perp_R) + 
				   std::norm(a_para_R) );
  
  result += ( 1. - betasq )*std::real( a_perp_L*
				       std::conj(a_perp_R) + 
				       a_para_L*
				       std::conj(a_para_R) );
  
  return result; 
}

double Amplitudes::I1c(const double ctk,
		       const double mhh,
		       const double s ) const {
  
  double result = 0;
  
  std::complex<double> a_zero_L(0,0);
  std::complex<double> a_zero_R(0,0);
  std::complex<double> a_time(0,0);

  AmplitudeSet amp;
  
  for ( const auto & r: m_res ){
    r->getAmplitudes( s, mhh, amp );

    const int l = r->getL();

    a_zero_L += amp.AZEROL * SphericalHarmonics::Y(l,0,ctk);
    a_zero_R += amp.AZEROR * SphericalHarmonics::Y(l,0,ctk);
    a_time += amp.ATIME * SphericalHarmonics::Y(l,0,ctk);
  }
  
  const double betasq = 
    Kinematic::betasq( s, Parameters::ml );

  result += ( std::norm(a_zero_L) + std::norm(a_zero_R) );
  
  result += ( 1. - betasq )*2.*std::real(a_zero_L*
					 std::conj(a_zero_R) );

  result += 4.0*Parameters::ml*Parameters::ml/s * std::norm(a_time);
  
  return result; 
}
	
double Amplitudes::I2s(const double ctk,
		       const double mhh,
		       const double s ) const {

  double result = 0;
  
  std::complex<double> a_perp_L(0,0);
  std::complex<double> a_perp_R(0,0);
  std::complex<double> a_para_L(0,0);
  std::complex<double> a_para_R(0,0);

  AmplitudeSet amp;
  
  for ( const auto & r: m_res ){

    r->getAmplitudes( s, mhh, amp );

    const int l = r->getL();
    
    a_perp_L += amp.APERPL * SphericalHarmonics::Y(l,-1,ctk);
    a_perp_R += amp.APERPR * SphericalHarmonics::Y(l,-1,ctk);
    a_para_L += amp.APARAL * SphericalHarmonics::Y(l,-1,ctk);
    a_para_R += amp.APARAR * SphericalHarmonics::Y(l,-1,ctk);

  }

  const double betasq = 
    Kinematic::betasq( s, Parameters::ml );
  
  result = (betasq/4.)*( std::norm(a_perp_L) +
			 std::norm(a_para_L) + 
			 std::norm(a_perp_R) +
			 std::norm(a_para_R) );
  
  return result; 
}

double Amplitudes::I2c( const double ctk,
			const double mhh,
			const double s ) const {
  
  double result = 0;
  
  std::complex<double> a_zero_L(0,0);
  std::complex<double> a_zero_R(0,0);
  
  AmplitudeSet amp;
  
  for ( const auto & r: m_res ){
    
    r->getAmplitudes( s, mhh, amp );
    
    const int l = r->getL();
    
    a_zero_L += amp.AZEROL * SphericalHarmonics::Y(l,0,ctk);
    a_zero_R += amp.AZEROR * SphericalHarmonics::Y(l,0,ctk);
    
  }
  
  const double betasq = 
    Kinematic::betasq( s, Parameters::ml );
  
  result = -betasq*( std::norm(a_zero_L) +
		     std::norm(a_zero_R) );
  
  return result; 
}

double Amplitudes::I3( const double ctk,
		       const double mhh,
		       const double s ) const {

  double result = 0;
  
  std::complex<double> a_perp_L(0,0);
  std::complex<double> a_perp_R(0,0);
  std::complex<double> a_para_L(0,0);
  std::complex<double> a_para_R(0,0);
  
  AmplitudeSet amp;
  
  for ( const auto & r: m_res ){
    
    r->getAmplitudes( s, mhh, amp );
    
    const int l = r->getL();
    
    a_perp_L += amp.APERPL * SphericalHarmonics::Y(l,-1,ctk);
    a_perp_R += amp.APERPR * SphericalHarmonics::Y(l,-1,ctk);
    a_para_L += amp.APARAL * SphericalHarmonics::Y(l,-1,ctk);
    a_para_R += amp.APARAR * SphericalHarmonics::Y(l,-1,ctk);  
  }
  
  const double betasq = 
    Kinematic::betasq( s, Parameters::ml );

  result += (betasq/2.)*( std::norm(a_perp_L) +
			  std::norm(a_perp_R) -
			  std::norm(a_para_L) -
			  std::norm(a_para_R) );
  
  return result; 
}


double Amplitudes::I4( const double ctk,
		       const double mhh,
		       const double s ) const {
  
  double result = 0;
  
  std::complex<double> a_zero_L(0,0);
  std::complex<double> a_zero_R(0,0);
  std::complex<double> a_para_L(0,0);
  std::complex<double> a_para_R(0,0);

  AmplitudeSet amp; 

  for ( const auto & r: m_res ){
    
    r->getAmplitudes( s, mhh, amp );
    
    const int l = r->getL();
    
    a_zero_L += amp.AZEROL * SphericalHarmonics::Y(l, 0,ctk);
    a_zero_R += amp.AZEROR * SphericalHarmonics::Y(l, 0,ctk);
    a_para_L += amp.APARAL * SphericalHarmonics::Y(l,-1,ctk);
    a_para_R += amp.APARAR * SphericalHarmonics::Y(l,-1,ctk);
    
  }
  
  const double betasq = 
    Kinematic::betasq( s, Parameters::ml );
  
  result += (betasq/std::sqrt(2.))*( std::real( a_zero_L*
						std::conj(a_para_L)    +
						a_zero_R*
						std::conj(a_para_R) ) );
  
  return result; 
}

double Amplitudes::I5( const double ctk,
		       const double mhh,
		       const double s ) const {
  
  double result = 0;
  
  std::complex<double> a_zero_L(0,0);
  std::complex<double> a_zero_R(0,0);
  std::complex<double> a_perp_L(0,0);
  std::complex<double> a_perp_R(0,0);
  
  AmplitudeSet amp;
  
  for ( const auto & r: m_res ){
    
    r->getAmplitudes( s, mhh, amp );
    
    const int l = r->getL();
    
    a_zero_L += amp.AZEROL * SphericalHarmonics::Y(l, 0,ctk);
    a_zero_R += amp.AZEROR * SphericalHarmonics::Y(l, 0,ctk);
    a_perp_L += amp.APERPL * SphericalHarmonics::Y(l,-1,ctk);
    a_perp_R += amp.APERPR * SphericalHarmonics::Y(l,-1,ctk);
    
  }
  
  const double betasq = 
     Kinematic::betasq( s, Parameters::ml );
  
  result += std::sqrt(2.*betasq)*( std::real( a_zero_L*
					      std::conj( a_perp_L ) - 
					      a_zero_R*
					      std::conj(a_perp_R) ) );
 
  return result; 
}
		
double Amplitudes::I6s( const double ctk,
		        const double mhh,
			const double s ) const {

  double result = 0;
  
  std::complex<double> a_perp_L(0,0);
  std::complex<double> a_perp_R(0,0);
  std::complex<double> a_para_L(0,0);
  std::complex<double> a_para_R(0,0);
  
  AmplitudeSet amp;
  
  for ( const auto & r: m_res ){
    
    r->getAmplitudes( s, mhh, amp );
    
    const int l = r->getL();
    
    a_perp_L += amp.APERPL * SphericalHarmonics::Y(l,-1,ctk);
    a_perp_R += amp.APERPR * SphericalHarmonics::Y(l,-1,ctk);
    a_para_L += amp.APARAL * SphericalHarmonics::Y(l,-1,ctk);
    a_para_R += amp.APARAR * SphericalHarmonics::Y(l,-1,ctk);
    
  }
  
  const double betasq = 
    Kinematic::betasq( s, Parameters::ml );
  
  result += 2.*std::sqrt(betasq)*( std::real( a_para_L*
					      std::conj(a_perp_L) - 
					      a_para_R*
					      std::conj(a_perp_R) ) );
  return result; 
}


double Amplitudes::I7( const double ctk,
		       const double mhh,
		       const double s ) const {
  
  double result = 0;
  
  std::complex<double> a_zero_L(0,0);
  std::complex<double> a_zero_R(0,0);
  std::complex<double> a_para_L(0,0);
  std::complex<double> a_para_R(0,0);
  
  AmplitudeSet amp;
  
  for ( const auto & r: m_res ){
    
    r->getAmplitudes( s, mhh, amp );
    
    const int l = r->getL();
    
    a_zero_L += amp.AZEROL * SphericalHarmonics::Y(l, 0,ctk);
    a_zero_R += amp.AZEROR * SphericalHarmonics::Y(l, 0,ctk);
    a_para_L += amp.APARAL * SphericalHarmonics::Y(l,-1,ctk);
    a_para_R += amp.APARAR * SphericalHarmonics::Y(l,-1,ctk);
    
  }
  
  const double betasq = 
    Kinematic::betasq( s, Parameters::ml );

  result += std::sqrt(2.*betasq)*( std::imag( a_zero_L*
					      std::conj(a_para_L) - 
					      a_zero_R*
					      std::conj(a_para_R) ) );
  
  return result; 
}


double Amplitudes::I8( const double ctk,
		       const double mhh,
		       const double s ) const {
  
  double result = 0;
  
  std::complex<double> a_zero_L(0,0);
  std::complex<double> a_zero_R(0,0);
  std::complex<double> a_perp_L(0,0);
  std::complex<double> a_perp_R(0,0);
  
  AmplitudeSet amp;
  
  for ( const auto & r: m_res ){
    
    r->getAmplitudes( s, mhh, amp );
    
    const int l = r->getL();
    
    a_zero_L += amp.AZEROL * SphericalHarmonics::Y(l, 0,ctk);
    a_zero_R += amp.AZEROR * SphericalHarmonics::Y(l, 0,ctk);
    a_perp_L += amp.APERPL * SphericalHarmonics::Y(l,-1,ctk);
    a_perp_R += amp.APERPR * SphericalHarmonics::Y(l,-1,ctk);
    
  }
  
  const double betasq = 
    Kinematic::betasq( s, Parameters::ml );

  result += (betasq/std::sqrt(2.))*( std::imag( a_zero_L*
						std::conj(a_perp_L) + 
						a_zero_R*
						std::conj(a_perp_R) ) );
  return result; 
}

double Amplitudes::I9( const double ctk,
		       const double mhh,
		       const double s ) const {
  
  double result = 0;
  
  std::complex<double> a_perp_L(0,0);
  std::complex<double> a_perp_R(0,0);
  std::complex<double> a_para_L(0,0);
  std::complex<double> a_para_R(0,0);
  
  AmplitudeSet amp;
  
  for ( const auto & r: m_res ){
    
    r->getAmplitudes( s, mhh, amp );
    
    const int l = r->getL();
    
    a_perp_L += amp.APERPL * SphericalHarmonics::Y(l,-1,ctk);
    a_perp_R += amp.APERPR * SphericalHarmonics::Y(l,-1,ctk);
    a_para_L += amp.APARAL * SphericalHarmonics::Y(l,-1,ctk);
    a_para_R += amp.APARAR * SphericalHarmonics::Y(l,-1,ctk);
    
  }

  const double betasq = 
    Kinematic::betasq( s, Parameters::ml );

  result += betasq*( std::imag( std::conj(a_para_L)*
				a_perp_L + 
				std::conj(a_para_R)*
				a_perp_R ) );
  
  return result; 
}

double Amplitudes::differential( const double s,
				 const double mhh,
				 const double ctl,
				 const double ctk,
				 const double phi ) const { 
 
  double result = ( 2.*I1s( ctk, mhh, s ) +
		    I1c( ctk, mhh, s ) + 
		    2.*I2s( ctk, mhh, s )*(2.*ctl*ctl - 1.) +
		    I2c( ctk, mhh, s )*(2.*ctl*ctl - 1.) + 
		    2.*I3 ( ctk, mhh, s )*(1. - ctl*ctl)*cos(2.*phi) + 
		    2.*std::sqrt(2)*I4 ( ctk, mhh, s )*(2.*ctl*std::sqrt(1.-ctl*ctl))*cos(phi) + 
		    2.*std::sqrt(2)*I5 ( ctk, mhh, s )*std::sqrt(1. - ctl*ctl)*cos(phi) + 
		    2.*I6s( ctk, mhh, s )*ctl + 
		    2.*std::sqrt(2)*I7 ( ctk, mhh, s )*std::sqrt(1. - ctl*ctl)*sin(phi) + 
		    2.*std::sqrt(2)*I8 ( ctk, mhh, s )*(2.*ctl*std::sqrt(1.-ctl*ctl))*sin(phi) + 
		    2.*I9 ( ctk, mhh, s )*(1. - ctl*ctl)*sin(2.*phi)) ; 
  
  bool isCP{false};
  for(const auto &r : m_res){
	  if(r->isCP() == true){ isCP = true;}
  }

  if(isCP){
  //Add the Bbar process where Is I5,I6,I8 and I9 flip sign
  //Currently ignoring interference term between B and Bbar
	result +=   ( 2.*I1s( ctk, mhh, s) + 
	            I1c( ctk, mhh, s ) + 
		    2.*I2s( ctk, mhh, s )*(2.*ctl*ctl - 1.) +
		    I2c( ctk, mhh, s )*(2.*ctl*ctl - 1.) + 
		    2.*I3 ( ctk, mhh, s )*(1. - ctl*ctl)*cos(2.*phi) + 
		    2.*std::sqrt(2)*I4 ( ctk, mhh, s )*(2.*ctl*sqrt(1.-ctl*ctl))*cos(phi) -
		    2.*std::sqrt(2)*I5 ( ctk, mhh, s )*sqrt(1. - ctl*ctl)*cos(phi) -
		    2.*I6s( ctk, mhh, s )*ctl + 
		    2.*std::sqrt(2)*I7 ( ctk, mhh, s )*sqrt(1. - ctl*ctl)*sin(phi) -
		    2.*std::sqrt(2)*I8 ( ctk, mhh, s )*(2.*ctl*sqrt(1.-ctl*ctl))*sin(phi) -
		    2.*I9 ( ctk, mhh, s )*(1. - ctl*ctl)*sin(2.*phi)) ; 

    result/=2.;
  }

  return ( 3.*result )/( 8.); 
}

double Amplitudes::operator()( const double* x ) const {
  /*
   * Compute the differential decay rate for given
   * coordinate set x. For use in toy generation. 
   */ 
  const double s = x[0]; 
  const double m = x[1];

  return differential( s, m, x[2], x[3], x[4] );
}

double Amplitudes::operator()( const std::vector<double> x) const {
  /*
   * Compute the differential decay rate for given
   * coordinate set x. For use in toy generation. 
   */ 
  const double s = x[0]; 
  const double m = x[1];
  
  return differential( s, m, x[2], x[3], x[4] );
}


double Amplitudes::integral( const std::string name, 
			     const double qsq_min, 
			     const double qsq_max, 
			     const double mhh_min, 
			     const double mhh_max ) const { 
  
  /*
   * Compute the integrated decay rate for the
   * named resonance.
   */ 
  
  for ( const auto & r: m_res ){
    if ( r->hasName( name ) ) { 
      return integral( r.get(), qsq_min, qsq_max, mhh_min, mhh_max ); 
    }
  }
  
  return 0; 
}

double Amplitudes::integral( const ResonanceBase* res, 
			     const double qsq_min, 
			     const double qsq_max, 
			     const double mhh_min, 
			     const double mhh_max ) const { 
  
  /*
   * Compute the integrated decay rate for resonance
   * ResonanceBase.
   */ 

  auto fn = std::bind( &Amplitudes::rate_fn, this, res, std::placeholders::_1 ); 

  ROOT::Math::Functor func( fn, 2 );
  
  double low[2] = { qsq_min, mhh_min }; 
  double upp[2] = { qsq_max, mhh_max };
  
  
  ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kPLAIN); 
  ig.SetFunction( func );
  
  double result = ig.Integral( low, upp );
  
  return result;
}





double Amplitudes::integral( const double qsq_min, 
			     const double qsq_max, 
			     const double mhh_min, 
			     const double mhh_max ) const { 

  /*
   * Compute the total integrated decay rate across all
   * resonances.
   */
  
  auto fn = std::bind( &Amplitudes::total_rate_fn, this, std::placeholders::_1 );
  
  ROOT::Math::Functor func( fn, 2 );
  
  double low[2] = { qsq_min, mhh_min }; 
  double upp[2] = { qsq_max, mhh_max };
  
  ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kPLAIN); 
  ig.SetFunction( func );
  
  double result = ig.Integral( low, upp );
  
  return result;
}

double Amplitudes::integralMC( const std::string name, 
			       const double qsq_min, 
			       const double qsq_max, 
			       const double mhh_min, 
			       const double mhh_max ) const { 
  
  /*
   * Compute the integrated decay rate for the
   * named resonance.
   */ 
  
  for ( const auto & r: m_res ){
    if ( r->hasName( name ) ) { 
      return integralMC( r.get(), qsq_min, qsq_max, mhh_min, mhh_max ); 
    }
  }
  
  return 0; 
}


double Amplitudes::integralMC( const ResonanceBase* res, 
			       const double qsq_min, 
			       const double qsq_max, 
			       const double mhh_min, 
			       const double mhh_max ) const {
  

  const unsigned int n = 1000000; 

  double result = 0; 

  for ( unsigned int i = 0; i < n; i++ ){
    
    double m = gRandom->Uniform( mhh_min, mhh_max );
    double s = gRandom->Uniform( qsq_min, qsq_max );
    
    result += rate( res, s, m );
  }
  
  result *= ( qsq_max - qsq_min );
  result *= ( mhh_max - mhh_min );
  result /= double( n );

  return result; 
} 
