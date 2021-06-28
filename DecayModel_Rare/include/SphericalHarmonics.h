#ifndef SPHERICALHARMONICS_H 
#define SPHERICALHARMONICS_H 1

//Class to store spherical harmonic terms used in the calculations
namespace SphericalHarmonics { 
  extern double Y( const int l, const int m, const double costheta );

  extern double sintheta( const double costheta );

  extern double Y00( const double costheta );

  extern double Y10( const double costheta );

  extern double Y11( const double costheta );

  extern double Y20( const double costheta );

  extern double Y21( const double costheta );

}

#endif
