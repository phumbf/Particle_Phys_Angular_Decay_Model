#ifndef KPIRESONANCEBW_H 
#define KPIRESONANCEBW_H

#include "KPiResonanceBase.h" 
#include "FitParameter.h"

#include <complex>
#include <string>


//Kaon + pion Resonance Breit-Wigner class - inherits from Base Resonance Class

class KPiResonanceBW : public KPiResonanceBase {
public:
    KPiResonanceBW( const double mean,
                    const double width,
                    const double mB,
                    const double mpsi,
                    const double mkaon,
                    const double mpion,
                    const unsigned int LR,
                    const unsigned int LB,
                    const double radius_ ) ;
    
    KPiResonanceBW( const double mean,
                    const double width,
                    const double mB,
                    const double mpsi,
                    const double mkaon,
                    const double mpion,
                    const unsigned int LR,
                    const unsigned int LB,
                    const double radiusR_,
                    const double radiusB_ ) ;
    
    
    KPiResonanceBW( const KPiResonanceBW& other ) ;
    
    virtual ~KPiResonanceBW() {} ;
    
    std::complex<double> lineshape( const double m, const double pB, const double pR ) const ;
    
    unsigned int spin() const ;

    double targetfitfrac() const ;
    
    KPiResonanceBase* clone() const {
        return new KPiResonanceBW( *this );
    }

    std::string getType() const { return "RelativisticBW"; }
    
    
private:
    double mean_;
    double width_;
    double LR_;
    double LB_;
    double radiusR_;
    double radiusB_;
    double targetfitfrac_;
    
    
};

inline unsigned int KPiResonanceBW::spin() const {
    return LR_ ;
}


#endif
