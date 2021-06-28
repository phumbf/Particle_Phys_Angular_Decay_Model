#include "SplineModel.h"
#include "TSpline.h"
#include <algorithm>
#include <iterator>

SplineModelKnot::SplineModelKnot( const double m ) :
  m_ ( m ),
  re_( 0 ),
  im_( 0 ) {} ;

SplineModelKnot::SplineModelKnot() :
  m_ ( 0 ),
  re_( 0 ),
  im_( 0 ) {} ;

SplineModelKnot::SplineModelKnot( const double m, const double re, const double im ) :
  m_ ( m ),
  re_( re ),
  im_( im ) {} ;

SplineModelKnot::SplineModelKnot( const SplineModelKnot& other ) :
  m_ ( other.m_  ),
  re_( other.re_ ),
  im_( other.im_ ) {} ;

double SplineModelKnot::mass() const {
  return m_;
}

double SplineModelKnot::real() const {
  return re_;
}

double SplineModelKnot::imag() const {
  return im_;
}

SplineModel::SplineModel() {} ;

SplineModel::SplineModel( const unsigned int knots, const double min, const double max ) {
  for ( unsigned int i = 0; i < knots; ++i ){
    knots_.push_back( SplineModelKnot( min + i*(max - min )/(knots - 1) ) );
  }
}

void SplineModel::addKnot( const double m ){ 
  knots_.push_back( SplineModelKnot( m ) );
  return; 
}

SplineModel::SplineModel( const SplineModel& other ){
  // copy the knots
  std::copy ( other.knots_.begin(), other.knots_.end(), std::back_inserter(knots_) );
}

std::complex<double> SplineModel::evaluate( const double m ) const {
  
  const unsigned int nknots = knots_.size();
  
  if ( nknots  < 3 ) return std::complex<double>(0,0);
  
  double mass[ nknots ];
  double real[ nknots ];
  double imag[ nknots ];
  
  for ( unsigned int i = 0; i < nknots ; ++i ){
    mass[i] = knots_[i].mass();
    real[i] = knots_[i].real();
    imag[i] = knots_[i].imag();
  }
  
  TSpline3 respline("Re",mass,real,nknots,0,0,0);
  TSpline3 imspline("Im",mass,imag,nknots,0,0,0);
  
  double re = respline.Eval( m );
  double im = imspline.Eval( m );
  
  return std::complex<double>( re, im );
}
