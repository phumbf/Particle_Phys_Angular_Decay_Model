#ifndef EVTGENFORMFACTOR_H
#define EVTGENFORMFACTOR_H

#include "FormFactorBase.h" 

//Class to calculate the Form Factor terms which contribute to the amplitudes
class EvtGenBsFormFactor : virtual public VectorFormFactorBase { 
 public: 

  VectorFormFactorBase* clone() const { 
    return  new EvtGenBsFormFactor();
  }
  
  double A1( const double qsq ) const override; 
  
  double A12( const double qsq ) const override { return 0;} ; 
  
  double A2( const double qsq ) const override; 
  
  double A0( const double qsq ) const override; 
  
  double V( const double qsq ) const override; 
  
  double T1( const double qsq ) const override; 
  
  double T2( const double qsq ) const override; 
  
  double T23( const double qsq ) const override { return 0; }; 

  double T3tilde( const double qsq ) const override; 
 
  double T3( const double qsq ) const override ;

  void setParams( const std::vector<double>& params ) override {} ;

  void setParams( const TVectorD& params ) override {} ; 
  
  TVectorD getParams() const { return TVectorD(); } 

  void print() const {} ;
  

 protected:

  double pole( const double qsq, const double r, const double m ) const ;
  
  double sqpole( const double qsq, const double r, const double m ) const;
}; 

class EvtGenBzFormFactor : virtual public VectorFormFactorBase { 
 public: 
  
  VectorFormFactorBase* clone() const { 
    return  new EvtGenBzFormFactor();
  }
  
  double A1( const double qsq ) const override; 
  
  double A12( const double qsq ) const override { return 0;} ; 
  
  double A2( const double qsq ) const override; 
  
  double A0( const double qsq ) const override; 
  
  double V( const double qsq ) const override; 
  
  double T1( const double qsq ) const override; 
  
  double T2( const double qsq ) const override; 
  
  double T23( const double qsq ) const override { return 0; }; 

  double T3tilde( const double qsq ) const override; 

  double T3( const double qsq ) const ;

  void setParams( const std::vector<double>& params ) override {} ;

  void setParams( const TVectorD& params ) override {} ; 
  
  TVectorD getParams() const { return TVectorD(); } 

  void print() const {} ;
  
 protected:
  double pole( const double qsq, const double r, const double m ) const ;
  
  double sqpole( const double qsq, const double r, const double m ) const;
    
}; 

#endif 
