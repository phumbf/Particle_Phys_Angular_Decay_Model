#ifndef FORMFACTORBASE_H
#define FORMFACTORBASE_H 1

#include "TVectorD.h"

//Classes for the form factors with the different scalar and vector terms 
class VectorFormFactorBase{ 
 public: 
  virtual ~VectorFormFactorBase() {} ; 

  virtual double A1( const double qsq ) const = 0; 
 
  virtual double A12( const double qsq ) const = 0; 
  
  virtual double A2( const double qsq ) const = 0; 
  
  virtual double A0( const double qsq ) const = 0; 
  
  virtual double V( const double qsq ) const = 0; 
  
  virtual double T1( const double qsq ) const = 0; 
  
  virtual double T2( const double qsq ) const = 0; 

  virtual double T3tilde( const double qsq ) const = 0; 

  virtual VectorFormFactorBase* clone() const = 0;
  
  virtual double T3( const double qsq ) const = 0; 
  
  virtual double T23( const double qsq ) const = 0; 

  virtual TVectorD getParams() const = 0;

  virtual void setParams( const TVectorD& params ) = 0; 
  
  virtual void setParams( const std::vector< double >& params ) = 0;
  
  virtual void print() const = 0;
  
  void setMass( const double value ) { m_mass  = value; } 

  void setMotherMass( const double value ){ m_motherMass = value; } 

  void setB0m( const double value ) { m_mB0m = value; } 

  void setB1m( const double value ) { m_mB1m = value; } 

  void setB1p( const double value ) { m_mB1p = value; } 
  
  void setBst( const double value ) { m_mBst = value; } 
  
  double m_mass{0};
  double m_motherMass{0};
  double m_mB0m{0};
  double m_mB1p{0};
  double m_mB1m{0};
  double m_mBst{0};
};


class ScalarFormFactorBase{ 
 public:
  
  virtual ~ScalarFormFactorBase() {} ; 

  virtual double fplus( const double qsq ) const = 0; 

  virtual double fminus( const double qsq ) const = 0;

  virtual double fzero( const double qsq ) const = 0;

  virtual double fT( const double qsq ) const = 0;

  virtual ScalarFormFactorBase* clone() const = 0; 

  virtual TVectorD getParams() const = 0;

  virtual void setParams( const TVectorD& params ) = 0; 
 
  virtual void setParams( const std::vector< double >& params ) = 0; 
  
  virtual void print() const = 0;
  
   void setMass( const int value ) { m_mass  = value; } 

  void setMotherMass( const int value ){ m_motherMass = value; } 

  void setB0m( const double value ) { m_mB0m = value; } 

  void setB1m( const double value ) { m_mB1m = value; } 

  void setB1p( const double value ) { m_mB1p = value; } 
  
  void setBst( const double value ) { m_mBst = value; } 
  
  double m_mass{0};
  double m_motherMass{0};
  double m_mB0m{0};
  double m_mB1p{0};
  double m_mB1m{0};
  double m_mBst{0};
};


#endif
