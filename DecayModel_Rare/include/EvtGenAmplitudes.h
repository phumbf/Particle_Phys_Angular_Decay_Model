#ifndef EVTGENAMPLITUDES_H 
#define EVTGENAMPLITUDES_H 1

#include <complex>
#include <vector>
#include "FormFactorBase.h" 
#include "Amplitudes.h" 
#include "TObject.h"
#include "TVectorD.h"

//Amplitudes class to mirror the treatment found in the EvtGen 
//software package https://evtgen.hepforge.org/

class EvtGenAmplitudes {

public: 
    
  EvtGenAmplitudes( const bool isBs, const bool isbtod );
  
  EvtGenAmplitudes( const EvtGenAmplitudes& other ) ;
  
  virtual ~EvtGenAmplitudes( ); ///< Destructor
  
  double rate( const double s ) const;
  
  double integral( const double min, const double max ) const ;
  
  double differential( const AmplitudeSet& amps, 
		       const double s, 
		       const double ctl,
		       const double ctk, 
		       const double phi ) const ; 
    
  double differential( const double s, 
		       const double ctl, 
		       const double ctk, 
		       const double phi ) const ;
  
  double operator()( const double* x ); 
  
  double operator()( const std::vector<double> &x ); 
  
  void setDeltaC9( const double reC9, const double imC9 ){ 
    C9_  = std::complex<double>( reC9, imC9 ); 
  }

  void setDeltaC10( const double reC10, const double imC10 ){ 
    C10_ = std::complex<double>( reC10, imC10 ); 
  }
 
 protected:  
  
  std::complex<double> AZERO( const double s, const double LR, const int flavour = 1 ) const ;
  
  std::complex<double> APARA( const double s, const double LR, const int flavour = 1 ) const ;
  
  std::complex<double> APERP( const double s, const double LR, const int flavour = 1 ) const ;
  
  std::complex<double> ATIME( const double s, const int flavour = 1 ) const ;
  
  double norm( const double s ) const ;

  double I1s( const AmplitudeSet& amps, const double s ) const;
  double I1c( const AmplitudeSet& amps, const double s ) const;
  double I2s( const AmplitudeSet& amps, const double s ) const;
  double I2c( const AmplitudeSet& amps, const double s ) const;
  double I3 ( const AmplitudeSet& amps, const double s ) const;
  double I4 ( const AmplitudeSet& amps, const double s ) const;
  double I5 ( const AmplitudeSet& amps, const double s ) const;
  double I6s( const AmplitudeSet& amps, const double s ) const;
  double I7 ( const AmplitudeSet& amps, const double s ) const;
  double I8 ( const AmplitudeSet& amps, const double s ) const;
  double I9 ( const AmplitudeSet& amps, const double s ) const;

 public:
  bool btod_;
  
 private:
  
  std::complex<double> C9_ ; 
  std::complex<double> C10_;
    
  VectorFormFactorBase*  FF_;
  
};

#endif 
