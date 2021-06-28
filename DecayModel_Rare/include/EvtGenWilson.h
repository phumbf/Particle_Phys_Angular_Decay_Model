#include <complex>

//To calculate Wilson Coefficients
namespace EvtGenWilson { 
  
  std::complex<double> GetC7Eff ( const double q2, const bool nnlo, const bool btod );
  
  std::complex<double> GetC9Eff ( const double q2, const bool nnlo, const bool btod );
  
  std::complex<double> GetC10Eff( const double q2, const bool nnlo, const bool btod );
} 
