#ifndef KPIRESONANCEBWNorm_H 
#define KPIRESONANCEBWNorm_H

#include "KPiResonanceBase.h" 
#include "FitParameter.h"

#include <complex>
#include <string>

//Similar Breit-Wigner Class with added Breit-Wigner normalisation term

class KPiResonanceBWNorm : public KPiResonanceBase {
public:
    KPiResonanceBWNorm( const double mean,
                    const double width,
                    const double mB,
                    const double mpsi,
                    const double mkaon,
                    const double mpion,
                    const unsigned int LR,
                    const unsigned int LB,
                    const double radius_ ) ;
    
    KPiResonanceBWNorm( const double mean,
                    const double width,
                    const double mB,
                    const double mpsi,
                    const double mkaon,
                    const double mpion,
                    const unsigned int LR,
                    const unsigned int LB,
                    const double radiusR_,
                    const double radiusB_ ) ;
    
    KPiResonanceBWNorm( const KPiResonanceBWNorm& other ) ;
    
    virtual ~KPiResonanceBWNorm() {} ;
    
    std::complex<double> lineshape( const double m, const double pB, const double pR ) const ;
    
    unsigned int spin() const ;

    
    KPiResonanceBase* clone() const {
        return new KPiResonanceBWNorm( *this );
    }
    
    std::string getType() const { return "NormRelativisticBW"; }

private:
    double mean_;
    double width_;
    double LR_;
    double LB_;
    double radiusR_;
    double radiusB_;
    
};

inline unsigned int KPiResonanceBWNorm::spin() const {
    return LR_ ;
}


#endif
