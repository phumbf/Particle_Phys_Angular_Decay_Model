#ifndef AMPLITUDESET_H
#define AMPLITUDESET_H

#include <complex>

//Class to store the amplitude polarisations
class AmplitudeSet { 
  
 public: 
  
  AmplitudeSet( const AmplitudeSet& other ) ; 
  
  AmplitudeSet( const std::complex<double> aZEROL,
		const std::complex<double> aZEROR,
		const std::complex<double> aPARAL,
		const std::complex<double> aPARAR,
		const std::complex<double> aPERPL,
		const std::complex<double> aPERPR,
		const std::complex<double> aTIME );
  
  AmplitudeSet() ; 
  
  ~AmplitudeSet() {} ; 
  
  void print() const ;

  void flipCP(const int spin);
  
 public:
  
  std::complex<double> AZEROL;
  std::complex<double> AZEROR;
  std::complex<double> APARAL;
  std::complex<double> APARAR;
  std::complex<double> APERPL;
  std::complex<double> APERPR;
  std::complex<double> ATIME ;  
}; 

#endif
