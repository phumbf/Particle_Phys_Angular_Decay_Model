// Local
#include "SphericalHarmonics.h" 

// std
#include <cmath> 

// ROOT
#include "TMath.h" 


//Class to store spherical harmonic terms used in the calculations
double SphericalHarmonics::sintheta( const double costheta ){ 
  return std::sqrt( 1. - std::pow(costheta,2) );
}

double SphericalHarmonics::Y00( const double /*costheta*/ ){
  return 1./std::sqrt( 4. * TMath::Pi() );
}

double SphericalHarmonics::Y10( const double costheta ){ 
  return std::sqrt( 3. / 4. / TMath::Pi() )*costheta; 
}

double SphericalHarmonics::Y11( const double costheta ){ 
  return -std::sqrt( 3. / 8. / TMath::Pi() ) * sintheta( costheta );
}

double SphericalHarmonics::Y20( const double costheta ){ 
  return std::sqrt( 5. / 16. / TMath::Pi() ) * ( 3. * std::pow( costheta, 2 ) - 1. );
}

double SphericalHarmonics::Y21( const double costheta ){ 
  return -std::sqrt( 15. / 32. / TMath::Pi() ) * 2. * sintheta( costheta ) * costheta; 
}


double SphericalHarmonics::Y( const int l, const int m, const double costheta ){ 
  
  if ( m < 0 ){ 
    return -SphericalHarmonics::Y( l, -m, costheta );
  }
  
  if ( 0 == l ) {
    if ( 0 == m ) return SphericalHarmonics::Y00( costheta );
  }
  if ( 1 == l ) { 
    if ( 0 == m ) return SphericalHarmonics::Y10( costheta );
    if ( 1 == m ) return SphericalHarmonics::Y11( costheta );
  }
  if ( 2 == l ) {
    if ( 0 == m ) return SphericalHarmonics::Y20( costheta );
    if ( 1 == m ) return SphericalHarmonics::Y21( costheta );
  }
  
  return 0; 
}


