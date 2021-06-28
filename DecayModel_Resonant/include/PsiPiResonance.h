#ifndef PSIPIRESONANCE_H
#define PSIPIRESONANCE_H


#include <complex>
#include <string> 

#include "Event.h"


/*
 * Model for (psi pi) resonances using a 
 * Relativistic BW.
 * 
 */

class PsiPiResonance {
public:
    
    PsiPiResonance( const double mean,
                   const double width,
                   const double mB,
                   const double mpsi,
                   const double mkaon,
                   const double mpion,
                   const unsigned int LR ,
                   const unsigned int LB ,
                   const double radiusR,
                   const double radiusB );
    
    PsiPiResonance( const PsiPiResonance& other );

    PsiPiResonance* clone() const {  return new PsiPiResonance( *this ); }
    
    virtual ~PsiPiResonance() {} ;

    std::complex<double> evaluate( const Event& event, const int lambda ) const ;
    
    std::complex<double> lineshape( const double m, const double pB, const double pR ) const;
    
    void setAmplitudesPolar( const double magA0  , const double magAp, const double magAm,
			     const double phaseA0, const double phaseAp, const double phaseAm ) ;
    
    void setAmplitudesCartesian( const double reA0, const double imA0, 
				 const double reAp, const double imAp, 
				 const double reAm, const double imAm  );
    
    std::string name() const ;
    
    void setName( const std::string& name ) ;
    
    bool isNamed() const ;
    
    bool isNamed( const std::string& name ) const ;
    
    bool isNameIn( const std::string& name ) const ;
    
    double project( const double m ) const ;

    
protected:
    
    double mean_;
    double width_ ;
    
    double mB_  ;
    double mpsi_;
    double mkaon_;
    double mpion_;
    
    unsigned int LR_ ;
    unsigned int LB_ ;
    
    double radiusR_ ;
    double radiusB_ ;
    
    double min_  ;
    double max_  ;
    
    std::complex<double> A0_;
    std::complex<double> Ap_;
    std::complex<double> Am_; 

    std::string name_ ;
    
};



#endif


