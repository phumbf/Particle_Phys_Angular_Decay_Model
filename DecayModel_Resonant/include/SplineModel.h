#ifndef SPLINEMODEL_H
#define SPLINEMODEL_H 

#include <vector>
#include <complex>

class SplineModelKnot {
public:
    SplineModelKnot() ;
    
    SplineModelKnot( const double m ) ;
    
    SplineModelKnot( const double m, const double re, const double im ) ;
    
    SplineModelKnot( const SplineModelKnot& other );
    
    virtual ~SplineModelKnot() {} ;
    
    double mass() const ;
    
    double real() const ;
    
    double imag() const ;
    
private:
    double m_;
    double re_;
    double im_;
};


class SplineModel {
public:
    SplineModel() ;
    
    SplineModel( const unsigned int knots, const double min, const double max ) ;
    
    SplineModel( const SplineModel& other );
    
    virtual ~SplineModel() {} ;
    
    void addKnot( const double m ) ;
    
    std::complex<double> evaluate( const double m ) const ;
    
private:
    std::vector< SplineModelKnot > knots_;
};


#endif
