#ifndef HHLINESHAPEBASE_H
#define HHLINESHAPEBASE_H 1

#include <complex>

//Base class to describe the particle invariant mass lineshape
class HHLineShapeBase{

 public:
  
  virtual std::complex<double> lineshape( const double m ) const = 0;
  
  HHLineShapeBase() = default;
  
  HHLineShapeBase(const HHLineShapeBase&) = default;
  
  HHLineShapeBase(HHLineShapeBase&&) = default;
  
  HHLineShapeBase& operator=(const HHLineShapeBase&) = default;
  
  HHLineShapeBase& operator=(HHLineShapeBase&&) = default;
  
  virtual ~HHLineShapeBase() = default;

  double normalise(const double mother, const double daug1, const double daug2) const;

  double shape( const double mkpi ) const;

};
#endif
