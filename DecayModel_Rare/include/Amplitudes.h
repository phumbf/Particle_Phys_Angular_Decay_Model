#ifndef AMPLITUDES_H 
#define AMPLITUDES_H 1

#include <complex>
#include <vector>
#include <memory>

#include "FFcalculation.h"
#include "FormFactorBase.h"
#include "Wilson.h"
#include "ResonanceBase.h"
#include "AmplitudeSet.h"

#include "TObject.h"
#include "TVectorD.h"

//Amplitude class built using the full amplitude model 

class ResonanceBase;

class Amplitudes {

public: 
  /// Standard constructor
    
  Amplitudes();
  
  //Amplitudes( const Amplitudes& other ) ;
  
  virtual ~Amplitudes( ); ///< Destructor
  
  void addResonance( std::unique_ptr<ResonanceBase>& r );
  
  double differential( const double s  ,
		       const double mhh,
		       const double ctl,
		       const double ctk,
		       const double phi ) const ;
    
  double operator()( const double* x ) const ; 
  
  double operator()( const std::vector<double> x) const ; 
  
  double FL  ( const std::string name, const double s, const double mhh ) const;
  
  double AFB  ( const std::string name, const double s, const double mhh ) const;
  
  double AT2 ( const std::string name, const double s, const double mhh ) const;

  
  double rate( const std::string name, const double s, const double mhh ) const;

  double rate( const double s, const double mhh ) const ;

  double integral( const std::string name, 
		   const double qsq_min, 
		   const double qsq_max, 
		   const double mhh_min, 
		   const double mhh_max ) const ;

  double integral( const double qsq_min, 
		   const double qsq_max, 
		   const double mhh_min, 
		   const double mhh_max ) const ;
  
  double integralMC( const std::string name, 
		     const double qsq_min, 
		     const double qsq_max, 
		     const double mhh_min, 
		     const double mhh_max ) const ;
  

  
  double rate_fn( const ResonanceBase* res, const double* x ) const ;

  double total_rate_fn( const double* x ) const ;

  ResonanceBase* getResonance( const std::string name ) const ;
  
 protected:  
  
  double FL  ( const ResonanceBase* res, const double s, const double mhh ) const;
  
  double AFB  ( const ResonanceBase* res, const double s, const double mhh ) const;

  double AT2 ( const ResonanceBase* res, const double s, const double mhh ) const;
  
  double rate( const ResonanceBase* res, const double s, const double mhh ) const;
  
  double norm( const double s, const double mV, const double mB ) const ;

  double integral( const ResonanceBase* res, 
		   const double qsq_min, 
		   const double qsq_max, 
		   const double mhh_min, 
		   const double mhh_max ) const ;
  
  double integralMC( const ResonanceBase* res, 
		     const double qsq_min, 
		     const double qsq_max, 
		     const double mhh_min, 
		     const double mhh_max ) const ;
  
  double I1s(const double ctk,
	     const double mhh,
	     const double s ) const;

  double I1c(const double ctk,
	     const double mhh,
	     const double s ) const;

  double I2s(const double ctk,
	     const double mhh,
	     const double s ) const;

  double I2c(const double ctk,
	     const double mhh,
	     const double s ) const;

  double I3 (const double ctk,
	     const double mhh,
	     const double s ) const;

  double I4 (const double ctk,
	     const double mhh,
	     const double s ) const;

  double I5 (const double ctk,
	     const double mhh,
	     const double s ) const;

  double I6s(const double ctk,
	     const double mhh,
	     const double s ) const;

  double I7 (const double ctk,
	     const double mhh,
	     const double s ) const;

  double I8 (const double ctk,
	     const double mhh,
	     const double s ) const;

  double I9 (const double ctk,
	     const double mhh,
	     const double s ) const;

 private:
  std::vector< std::unique_ptr<ResonanceBase> > m_res; 
  
};

#endif // AMPLITUDES_VLLAMPLITUES_H
