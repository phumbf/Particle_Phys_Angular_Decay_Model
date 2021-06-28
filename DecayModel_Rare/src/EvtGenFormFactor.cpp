#include "EvtGenFormFactor.h" 
#include "Parameters.h" 

#include <cmath> 

//Class to calculate the Form Factor terms which contribute to the amplitudes
double EvtGenBsFormFactor::T3( const double qsq ) const { 
  
  const double t3tilde = T3tilde( qsq ); 
  const double t2      = T2( qsq );
  
  double t3 = 0.0;
  
  if ( std::fabs(qsq) > 1e-10) {
    t3 = ( std::pow(Parameters::mB_Bs,2) - std::pow(Parameters::mV,2) )*(t3tilde - t2)/qsq;
  }
  
  return t3; 
}

double EvtGenBzFormFactor::T3( const double qsq ) const { 
  
  const double t3tilde = T3tilde( qsq ); 
  const double t2      = T2( qsq );
  
  double t3 = 0.0;
  
  if ( std::fabs(qsq) > 1e-10) {
    t3 = ( std::pow(Parameters::mB_Bz,2) - std::pow(Parameters::mV,2) )*(t3tilde - t2)/qsq;
  }
  
  return t3; 
}

double EvtGenBsFormFactor::A1( const double qsq ) const { 
  return pole( qsq, 0.231, std::sqrt( 32.94 ) );
}

double EvtGenBsFormFactor::A2( const double qsq ) const { 
  return ( pole( qsq  , -0.011, std::sqrt( 40.14 ) ) +
	   sqpole( qsq,  0.192, std::sqrt( 40.14 ) ) ); 
}

double EvtGenBsFormFactor::A0( const double qsq ) const { 
  return  ( pole( qsq,  2.813, 5.37 ) + 
	    pole( qsq, -2.509, std::sqrt( 31.58 ) ) );
}

double EvtGenBsFormFactor::V( const double qsq ) const { 
  return ( pole( qsq,  2.351, 5.42 ) + 
	   pole( qsq, -2.039, std::sqrt(33.10) ) );
}

double EvtGenBsFormFactor::T1( const double qsq ) const { 
  return ( pole( qsq,  2.047, 5.42 ) + 
	   pole( qsq, -1.787, std::sqrt( 32.83 ) ) );
}

double EvtGenBsFormFactor::T2( const double qsq ) const { 
  return pole( qsq, 0.260, std::sqrt( 33.01 ) ); 
}

double EvtGenBsFormFactor::T3tilde( const double qsq ) const { 
  return ( pole( qsq  , 0.043, std::sqrt( 39.38 ) ) + 
	   sqpole( qsq, 0.217, std::sqrt( 39.38 ) ) ); 
}


double EvtGenBsFormFactor::pole( const double qsq, const double r, const double m ) const { 
  return r/( 1. - qsq/std::pow(m,2) );
}


double EvtGenBsFormFactor::sqpole( const double qsq, const double r, const double m ) const { 
  return r/std::pow( 1. - qsq/std::pow( m, 2 ) , 2 );
}

double EvtGenBzFormFactor::pole( const double qsq, const double r, const double m ) const { 
  return r/( 1. - qsq/std::pow(m,2) );
}


double EvtGenBzFormFactor::sqpole( const double qsq, const double r, const double m ) const { 
  return r/std::pow( 1. - qsq/std::pow( m, 2 ) , 2 );
}

double EvtGenBzFormFactor::A1( const double qsq ) const { 
  return pole( qsq, 0.290, std::sqrt(40.38) ); 
}


double EvtGenBzFormFactor::A2( const double qsq ) const { 
  return ( pole( qsq  , -0.084, std::sqrt( 52.00 ) ) + 
	   sqpole( qsq,  0.342, std::sqrt( 52.00 ) ) ); 
}

double EvtGenBzFormFactor::A0( const double qsq ) const { 
  return ( pole( qsq,  1.364, 5.28 ) + 
	   pole( qsq, -0.990, std::sqrt( 36.78 ) ) ); 
}

double EvtGenBzFormFactor::V( const double qsq ) const { 
  return ( pole( qsq, 0.923, 5.32 ) + 
	   pole( qsq, -0.511, std::sqrt( 49.40 ) ) );
}

double EvtGenBzFormFactor::T1( const double qsq ) const { 
  return ( pole( qsq,  0.823, 5.32 ) + 
	   pole( qsq, -0.491, std::sqrt( 46.31 ) ) ); 
}

double EvtGenBzFormFactor::T2( const double qsq ) const { 
  return pole( qsq, 0.332, std::sqrt( 41.41 ) );
}

double EvtGenBzFormFactor::T3tilde( const double qsq ) const { 
  return ( pole( qsq  , -0.036, std::sqrt( 48.10 ) ) + 
	   sqpole( qsq,  0.368, std::sqrt( 48.10 ) ) );
}

