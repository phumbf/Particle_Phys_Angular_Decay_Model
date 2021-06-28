
#include "KPiResonanceNR.h"
#include "boost/multi_array.hpp"

KPiResonanceNR::KPiResonanceNR( const double mB,
                                const double mpsi,
                                const double mkaon,
                                const double mpion ) :
KPiResonanceBase( mB, mpsi, mkaon, mpion ) {} ;

KPiResonanceNR::KPiResonanceNR( const KPiResonanceNR& other ) :
KPiResonanceBase( other ) {} ;


std::complex<double> KPiResonanceNR::lineshape( const double /*m*/, const double pB, const double /*pR*/ ) const {
    return std::complex<double>( pB/mB_, 0 );
}



bool KPiResonanceNR::setParameter( const double /*val*/, const std::string& par ) {
    
    std::cout << " ==> INFO (KPiResonanceNR): Unknown parameter " << par << std::endl;
    return false;
}


bool KPiResonanceNR::hasChanged() const {
    return false;
}
