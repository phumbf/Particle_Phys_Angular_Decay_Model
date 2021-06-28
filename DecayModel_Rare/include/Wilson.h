#ifndef WILSON_H 
#define WILSON_H 1

#include <complex> 
#include <array>
#include <vector>

//Class to store the Wilson Coefficients 
class ResonantState { 
 public:
  ResonantState( const double mpole,
		 const double gamma,
		 const double phase,
		 const double magnitude );
  
  ResonantState( const ResonantState& other ) ;
  
  double mass;
  double width;
  double arg;
  double scale;
};

class Wilson {
public: 
  
  /// Standard constructor
  enum Param { SM,varyC9,varyC9C10, Zero };  
  
  Wilson( ); 
  
  Wilson( const Wilson& other ) ; 

  virtual ~Wilson( ); ///< Destructor

  std::complex<double>  operator()(const unsigned int i) const;
  std::complex<double>& operator[](const unsigned int i) ;

  Wilson& operator=(const Wilson& other);

  void print() const ;

  void init( Wilson::Param param );

  std::complex<double> Y( const double s ) const;
  
  std::complex<double> Ysum( const double s, const std::vector< ResonantState >& res ) const;
  
protected:

  std::complex<double> h ( const double s, const double m ) const ;
  std::complex<double> h0( const double s ) const;

  std::complex<double> resonance( const double s,
				  const double mean,
				  const double gamma,
				  const double arg,
				  const double scale ) const; 

private:
  std::array< std::complex< double >, 11 > par_;

};
#endif // WILSON_H
