#include "KinematicFunctions.h" 
#include <cmath> 
#include <iostream>

//useful kinematic terms relevant for the maths
//
double KinematicFunctions::meffective( const double mean, const double min, const double max )  {
  return 0.5*(min + ( max - min )*(1. + std::tanh( (mean - 0.5*( min + max ))/(max -min ))));
}

double KinematicFunctions::momentum( const double m, const double m1, const double m2 ) {
    double msq = m*m;
    
    double psq = 0.25*(msq - (m1 + m2)*(m1 + m2))*(msq - (m1 - m2)*(m1 - m2))/msq;
    if(psq < 0){
	std::cout << "ERROR: psq variable is negative in Helper::momentum. Probably using the wrong kinematical variables in Event. Please check." << std::endl;
	std::cout << "DEBUG KINEMATIC FUNCTIONS: msq, m1, m2" << msq << " " << m1 << " " << m2 << std::endl;
	return 0;
    }
    else{
	    return std::sqrt( psq );
    }
}

double KinematicFunctions::parity( const int pz, const unsigned int jz ){
  return -pz*std::pow( -1, jz );
} 

double KinematicFunctions::barrier( const double p0, const double p , const double radius, const unsigned int LX ){
    
  const double z0 = p0*p0*radius*radius;
  const double z  = p*p*radius*radius;
  double result = 0;
  
  switch ( LX ) {
  case 0:
    result = 1.0;
    break;
  case 1:
    result = std::sqrt( (1. + z0 )/
			(1. + z ) );
    break;
  case 2:
    result = std::sqrt( (9. + 3.*z0 + z0*z0 )/
			(9. + 3.*z  + z*z ) );
    break;
  case 3:
    result = std::sqrt( (225. + 45.*z0 + 6.*z0*z0 + z0*z0*z0)/
			(225. + 45.*z  + 6.*z*z   + z*z*z) );
    break;
  case 4:
    result = std::sqrt( (11025. + 1575.*z0 + 135.*z0*z0 + 10.*z0*z0*z0 + z0*z0*z0*z0)/
			(11025. + 1575.*z  + 135.*z*z   + 10.*z*z*z    + z*z*z*z) );
    break;
  case 5:
    result = std::sqrt( (893025. + 99225.*z0 + 6300.*z0*z0 + 315.*z0*z0*z0 + 15*z0*z0*z0*z0 + z0*z0*z0*z0*z0 )/
			(893025. + 99225.*z  + 6300.*z*z   + 315.*z*z*z    + 15*z*z*z*z     + z*z*z*z*z    ) );
    break;
  default:
    result = 0;
    break;
  }
  return result;
}

double KinematicFunctions::beta( const double msq, const double m ){ 
   return std::sqrt( 1.0 - 4.*m*m/msq );
}
