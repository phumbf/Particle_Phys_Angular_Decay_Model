#ifndef KPIRESONANCELASS_H 
#define KPIRESONANCELASS_H

#include "KPiResonanceBase.h"

#include <complex>
#include <string>

//Derived class to implement a LASS lineshape
class KPiResonanceLASS : public KPiResonanceBase {
public:
    KPiResonanceLASS( const double mean,
                     const double width,
                     const double mB,
                     const double mpsi,
                     const double mkaon,
                     const double mpion,
                     const double radius,
                     const double a,
                     const double r,
                     const double cutoff ) ;

    KPiResonanceLASS( const double mB,
		      const double mpsi,
		      const double mkaon,
		      const double mpion ) ;
    
    KPiResonanceLASS( const KPiResonanceLASS& other ) ;
    
    virtual ~KPiResonanceLASS() {} ;
    
    std::complex<double> lineshape( const double m, const double pB, const double pR ) const ;
    
    unsigned int spin() const ;
    
    KPiResonanceBase* clone() const {
        return new KPiResonanceLASS( *this );
    }
    
    std::string getType() const { return "LASS"; } 

private:
    double mean_;
    double width_;
    double radius_;
    double a_;
    double r_;
    double cutoff_;
};


inline unsigned int KPiResonanceLASS::spin() const {
    return 0;
}


#endif
