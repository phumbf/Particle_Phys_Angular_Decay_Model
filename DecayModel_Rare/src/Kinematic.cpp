#include "Kinematic.h"
#include <cmath>

//Useful kinematic terms used to determine the form factor variation with momenta
double Kinematic::z( const double qsq, 
		     const double m, 
		     const double mB ){ 

  const double tp = std::pow( mB + m , 2 );
  const double tm = std::pow( mB - m , 2 );
    
  const double tz = tp - std::sqrt( tp*(tp  - tm) );
    
  return ( (std::sqrt( tp - qsq ) - std::sqrt(tp - tz))/
	   (std::sqrt( tp - qsq ) + std::sqrt(tp - tz)) );
}

double Kinematic::betasq( const double qsq, 
			  const double m ){ 

  return  1.0 - 4.*std::pow( m, 2 )/qsq ; 
}

double Kinematic::beta( const double qsq, 
			const double m ) { 
  return std::sqrt( 1.0 - 4.*std::pow( m, 2 )/qsq );
}

double Kinematic::lambda( const double qsq, 
			  const double m, 
			  const double mB ){ 
  double mBsq = std::pow( mB, 2 );
  double msq  = std::pow( m , 2 );
    
  return ( mBsq*mBsq + msq*msq + qsq*qsq - 2.*(mBsq*msq + msq*qsq + mBsq*qsq ) );
}

double Kinematic::zseries( const double qsq,
			   const double m ,
			   const double mB,
			   const double a0,
			   const double a1,
			   const double a2 ){ 
  double result = a0;
  
  const double zdiff = ( Kinematic::z( qsq, m, mB ) - 
			 Kinematic::z( 0  , m, mB ) );

  result += a1*zdiff ;
  result += a2*std::pow( zdiff, 2 ); 
  
  return result;
}

double Kinematic::pole( const double qsq, const double mpole ){
  return 1./(1. - qsq/(std::pow(mpole,2)));
}

double Kinematic::shat( const double qsq, const double mB ){ 
  return qsq/std::pow( mB, 2 );
}

double Kinematic::ffparam( const double shat, 
			   const double f, 
			   const double a, 
			   const double b ){
  return f/(1. - a*shat + b*shat*shat );
}


bool Kinematic::isForbidden( const double qsq,
			     const double mlepton, 
			     const double mhadron, 
			     const double mhone  ,
			     const double mhtwo  , 
			     const double mB ){
  return ( ( qsq < std::pow( 2.*mlepton  , 2 ) ) || 
	   ( qsq > std::pow( mB - mhadron, 2 ) ) || 
	   ( mhadron < ( mhone + mhtwo ) ) ||
	   ( mhadron > ( mB - std::sqrt(qsq) ) ) );
}
