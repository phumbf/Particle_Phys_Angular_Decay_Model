#ifndef KPIRESONANCENR_H
#define KPIRESONANCENR_H

#include "KPiResonanceBase.h"

#include <complex>
#include <string>

//An experimental lineshape derived class
class KPiResonanceNR : public KPiResonanceBase {
public:
    KPiResonanceNR( const double mB,
                   const double mpsi,
                   const double mkaon,
                   const double mpion ) ;
    
    KPiResonanceNR( const KPiResonanceNR& other ) ;
    
    virtual ~KPiResonanceNR() {} ;
    
    std::complex<double> lineshape( const double m, const double pB, const double pR ) const ;
    
    unsigned int spin() const ;
    
    bool setParameter( const double val, const std::string& par );
    
    bool hasChanged() const ;
    
    KPiResonanceBase* clone() const {
        return new KPiResonanceNR( *this );
    }
    std::string getType() const { return "NR"; }
};

inline unsigned int KPiResonanceNR::spin() const {
    return 0;
}

#endif
