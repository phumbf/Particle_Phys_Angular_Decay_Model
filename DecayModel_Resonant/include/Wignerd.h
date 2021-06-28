#include <cmath>

//Wigner d function values used in the calculations
namespace Wignerd {
    
    double d000 ( const double costheta );
    
    double d100 ( const double costheta );
    double d110 ( const double costheta );
    double d11p1( const double costheta );
    double d11m1( const double costheta );
    
    double d22p2( const double costheta );
    double d22p1( const double costheta );
    double d220 ( const double costheta );
    double d22m2( const double costheta );
    double d22m1( const double costheta );
    double d21p1( const double costheta );
    double d21m1( const double costheta );
    double d210 ( const double costheta );
    double d200 ( const double costheta );
    
    double d33p3( const double costheta );
    double d33p2( const double costheta );
    double d33p1( const double costheta );
    double d330 ( const double costheta );
    double d33m1( const double costheta );
    double d33m2( const double costheta );
    double d33m3( const double costheta );
    double d32p2( const double costheta );
    double d32p1( const double costheta );
    double d320 ( const double costheta );
    double d32m2( const double costheta );
    double d32m1( const double costheta );
    double d31p1( const double costheta );
    double d31m1( const double costheta );
    double d310 ( const double costheta );
    double d300 ( const double costheta );
    
    
    double d44p4( const double costheta );
    double d44p3( const double costheta );
    double d44p2( const double costheta );
    double d44p1( const double costheta );
    double d440 ( const double costheta );
    double d44m1( const double costheta );
    double d44m2( const double costheta );
    double d44m3( const double costheta );
    double d44m4( const double costheta );
    double d43p3( const double costheta );
    double d43p2( const double costheta );
    double d43p1( const double costheta );
    double d430 ( const double costheta );
    double d43m1( const double costheta );
    double d43m2( const double costheta );
    double d43m3( const double costheta );
    double d42p2( const double costheta );
    double d42p1( const double costheta );
    double d420 ( const double costheta );
    double d42m2( const double costheta );
    double d42m1( const double costheta );
    double d41p1( const double costheta );
    double d41m1( const double costheta );
    double d410 ( const double costheta );
    double d400 ( const double costheta );
    
    
    double djmpm( const double costheta, const unsigned int j,  int mp,  int m );
    
    
    double generic( const double beta, const int j, int mp, int m );
    

    template<int j, int mp, int n> struct dfunction {
        double operator()( const double /*costheta*/ ) { return 0; }
    };
    
    template<> struct dfunction<0,0,0> {
        double operator()( const double /*costheta*/ ) { return 1.0; }
    };
    
    template<> struct dfunction<1,0,0> {
        double operator()( const double costheta ) {
            return costheta;
        }
    };
    
    template<> struct dfunction<1,1,0> {
        double operator()( const double costheta ) {
            return -std::sqrt( (1. - costheta*costheta )/2. );
        }
    };

    
    template<> struct dfunction<1,1,1> {
        double operator()( const double costheta ) {
            return ( 1. + costheta )/2. ;
        }
    };
                         
    template<> struct dfunction<1,1,-1> {
        double operator()( const double costheta ) {
            return ( 1. - costheta )/2. ;
        }
    };
    
    template<> struct dfunction<2,2,2> {
        double operator()( const double costheta ) {
            return (1. + costheta)*(1. + costheta)/4.;
        }
    };
    
    template<> struct dfunction<2,2,1> {
        double operator()( const double costheta ) {
            double sintheta = std::sqrt( 1. - costheta*costheta );
            return -(1. + costheta)*sintheta/2.;
        }
    };

    template<> struct dfunction<2,2,0> {
        double operator()( const double costheta ) {
            return costheta*costheta*std::sqrt(6.)/4.;
        }
    };
    
    template<> struct dfunction<2,2,-1> {
        double operator()( const double costheta ) {
            double sintheta = std::sqrt( 1. - costheta*costheta );
            return -(1. - costheta)*sintheta/2.;
        }
    };
    
    template<> struct dfunction<2,2,-2> {
        double operator()( const double costheta ) {
            return (1. - costheta)*(1. - costheta)/4.;
        }
    };
    
    template<> struct dfunction<2,1,1> {
        double operator()( const double costheta ) {
            return ( 1. + costheta )*(2.*costheta - 1. )/2.;
        }
    };
    
    template<> struct dfunction<2,1,-1> {
        double operator()( const double costheta ) {
            return ( 1. - costheta )*(2.*costheta + 1. )/2.;
        }
    };
    
    template<> struct dfunction<2,1,0> {
        double operator()( const double costheta ) {
            double sintheta = std::sqrt( 1. - costheta*costheta );
            return -sintheta*costheta*std::sqrt(3./2.);
        }
    };
    
    template<> struct dfunction<2,0,0> {
        double operator()( const double costheta ) {
            return (3./2)*costheta*costheta - (1./2.);
        }
    };



}
