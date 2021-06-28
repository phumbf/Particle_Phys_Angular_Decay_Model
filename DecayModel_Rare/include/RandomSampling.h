#ifndef RANDOMSAMPLING_H 
#define RANDOMSAMPLING_H

#include "TVectorD.h"
#include "TMatrixDSym.h"
#include "TRandom3.h"

class RandomSampling {
    
public:
    
    RandomSampling( const TVectorD& vec, const TMatrixDSym& cov, const TMatrixDSym& inv );
    
    RandomSampling( const TVectorD& vec, const TMatrixDSym& cov );
    
    RandomSampling( const TVectorD& vec, const TVectorD& errs);
    
    RandomSampling( const RandomSampling& other );
    
    ~RandomSampling() ;
    
    TVectorD parameters() const ;
    
    TVectorD random() ;
    
    TVectorD uncorrelated() ;
    
    void seed( const UInt_t num ) ;
    
protected:
    
    TVectorD    vec_ ;
    TVectorD    errs_ ;
    TMatrixDSym cov_ ;
    TMatrixDSym inv_ ;
    
};




class RandomBSZ15FormFactor : RandomSampling {
public:
    RandomBSZ15FormFactor( const TVectorD& vec, const TMatrixDSym& cov, const TMatrixDSym& inv, const double mothermass, const double resmass );
    
    RandomBSZ15FormFactor( const TVectorD& vec, const TMatrixDSym& cov, const double mothermass, const double resmass );
    
    RandomBSZ15FormFactor( const RandomBSZ15FormFactor& other );
    
    ~RandomBSZ15FormFactor() ;
    
    TVectorD parameters() const ;
    
    TVectorD random() ;
    
    TVectorD uncorrelated() ;
    
    const double m_mother;

    const double m_mass;
    
private:
    
    TVectorD insert( const TVectorD& vec ) const ;
    
};


#endif
