#ifndef FFPARSER_HPP
#define FFPARSER_HPP

// std
#include <vector>
#include <set>

//Root
#include "TMatrixDSym.h"
#include "TVectorD.h"

// local
#include "FFcalculation.h" 

//Class to parse in the form factor information from json files
class FFparser {
 
 public:
  
  FFparser();
  
  ~FFparser();
  
  std::unique_ptr<YANG10FFCalculation> getYANG10VectorFormFactors( const std::string filename ) const ;
  
  std::unique_ptr<ESS18FFCalculation> getESS18VectorFormFactors( const std::string filename ) const ;
  
  std::unique_ptr<MS00FFCalculation> getMS00VectorFormFactors( const std::string filename ) const ;
  
  std::unique_ptr<ABHH99FFCalculation> getABHH99VectorFormFactors( const std::string filename ) const ;
  
  std::unique_ptr<BSZ15FFCalculation>  getBSZ15VectorFormFactors ( const std::string filename ) const;
  
  std::unique_ptr<BZ04FFCalculation>   getBZ04VectorFormFactors  ( const std::string filename ) const;

  std::unique_ptr<AAS07FFCalculation>  getAAS07ScalarFormFactors ( const std::string filename ) const ;

  std::unique_ptr<CFW10FFCalculation>  getCFW10ScalarFormFactors ( const std::string filename ) const ;
  
  std::vector<double> getFabFormFactors( const std::string filename) const ;
  
  std::vector<double> getCCCFormFactors( const std::string filename) const ;

  TMatrixDSym getBSZ15Covariance( const std::string filename ) const;
  
  TVectorD getfabUncertainty( const std::string filename ) const;

  TVectorD getaaaUncertainty( const std::string filename ) const;

  TMatrixDSym dropParams( const TMatrixDSym& M, const std::set< int >& droppedParams ) const;

  TVectorD    dropParams( const TVectorD& V, const std::set< int >& droppedParams ) const;

  
  void setDebug( const bool flag ){ m_debug = flag; } 
  
 private:
  bool m_debug{false};

};


#endif
