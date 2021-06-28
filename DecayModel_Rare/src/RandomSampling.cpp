#include "RandomSampling.h" 
#include "Parameters.h"

#include <cmath>
#include <cassert>
#include <iostream>

#include "TDecompSVD.h"
#include "TDecompChol.h"

#include "TRandom3.h"


RandomBSZ15FormFactor::RandomBSZ15FormFactor( const TVectorD& vec, const TMatrixDSym& cov, const double mothermass, const double resmass ) :
    m_mother(mothermass),
    m_mass(resmass),
    RandomSampling( vec, cov ) {}  

RandomBSZ15FormFactor::RandomBSZ15FormFactor( const TVectorD& vec, const TMatrixDSym& cov, const TMatrixDSym& inv, const double mothermass, const double resmass ) :
    m_mother(mothermass),
    m_mass(resmass),
    RandomSampling( vec, cov, inv ) {}  

RandomBSZ15FormFactor::RandomBSZ15FormFactor( const RandomBSZ15FormFactor& other ) :
    m_mother(other.m_mother),
    m_mass(other.m_mass),
    RandomSampling( other ) {} 


RandomBSZ15FormFactor::~RandomBSZ15FormFactor() { } 

TVectorD RandomBSZ15FormFactor::insert( const TVectorD& vec )  const {
  assert( 19 == vec.GetNrows() );
  
  TVectorD par( 21 );
  
  for ( unsigned int i = 0; i < 6; ++i ) {
    par(i) = vec(i);
  }
  
  const double mBsq  = std::pow( m_mother, 2 );
  const double mVsq  = std::pow( m_mass, 2 );
  
  
  const double scale = (mBsq - mVsq)/(8.*(m_mother)*(m_mass));
  
  par(6) = scale*vec(0);
  
  for ( unsigned int i = 7; i < 15; ++i ) {
    par(i) = vec(i-1);
  }
  
  par(15) = vec(11);
  
  for ( unsigned int i = 16; i < 21; ++i ){
    par(i) = vec(i-2);
  }
  
  return par;
}



TVectorD RandomBSZ15FormFactor::parameters() const {
  return insert( vec_ ) ;
}

TVectorD RandomBSZ15FormFactor::random() {
  TVectorD par = RandomSampling::random() ;
  return insert( par );
}

TVectorD RandomBSZ15FormFactor::uncorrelated(){
    TVectorD par = RandomSampling::uncorrelated() ;
    return insert( par );
}

RandomSampling::RandomSampling( const TVectorD& vec, const TVectorD& errs) : 
	vec_ ( vec ),
	errs_ ( errs){}

RandomSampling::RandomSampling( const TVectorD& vec, const TMatrixDSym& cov ) :
    vec_ ( vec ),
    cov_ ( cov ),
    inv_ ( cov ){
    
    inv_.Invert();
} 

RandomSampling::RandomSampling( const TVectorD& vec, const TMatrixDSym& cov, const TMatrixDSym& inv ) :
    vec_ ( vec ),
    cov_ ( cov ),
    inv_ ( inv ){} 


RandomSampling::RandomSampling( const RandomSampling& other ) :
    vec_ ( other.vec_ ),
    cov_ ( other.cov_ ),
    inv_ ( other.inv_ ){} 


RandomSampling::~RandomSampling() { } 

void RandomSampling::seed( const UInt_t num ) {
    gRandom->SetSeed( num );
}

TVectorD RandomSampling::parameters() const {
    return vec_;
}


TVectorD RandomSampling::uncorrelated() {
  TVectorD par( vec_.GetNrows() );
  
  for ( int i = 0; i < vec_.GetNrows(); ++i ) {
    par(i) = gRandom->Gaus( vec_(i), std::sqrt( cov_(i,i) ) );
  }
  
  return par;
}

TVectorD RandomSampling::random(){
    
  const int nrows = vec_.GetNrows();
  
  if ( 1 == nrows ){ 
    return uncorrelated();
  }

  TVectorD par( nrows );
  
  TVectorD tmp( nrows );
  
  for ( int i = 0; i < nrows; ++i ){
    tmp(i) = gRandom->Gaus( 0., 1. );
  }
  
  TDecompChol cholesky( cov_ );
  bool success = cholesky.Decompose();
  
  assert(success);
  
  TMatrixD U  = cholesky.GetU();
  TMatrixD UT( U );
  
  UT.Transpose( U );
  
  par = vec_ + UT*tmp;
  
  return par;
}
