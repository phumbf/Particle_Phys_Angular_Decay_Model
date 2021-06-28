#include "Wilson.h"
#include "Parameters.h"
#include <cassert>
#include <iostream>
#include "TMath.h"

//Class to store the Wilson Coefficients 
Wilson::Wilson(){
}

Wilson::Wilson( const Wilson& other ){
  
  for ( unsigned i = 0; i < other.par_.size(); ++i ){ 
    par_[i] = other.par_[i];
  }
}

void Wilson::init( Wilson::Param param )
{
  if ( param == Wilson::SM ){
    par_[1]  = std::complex<double>( -0.257 ,0. );
    par_[2]  = std::complex<double>(  1.009 ,0. );
    par_[3]  = std::complex<double>( -0.005 ,0. );
    par_[4]  = std::complex<double>( -0.078 ,0. );
    par_[5]  = std::complex<double>(  0.000 ,0. );
    par_[6]  = std::complex<double>(  0.001 ,0. );
    par_[7]  = std::complex<double>( -0.304 ,0. );
    par_[8]  = std::complex<double>( -0.167 ,0. );
    par_[9]  = std::complex<double>(  4.211 ,0. );
    par_[10] = std::complex<double>( -4.103 ,0. );
  }
  else if (param == Wilson::varyC9){
    par_[1]  = std::complex<double>( -0.257 ,0. );
    par_[2]  = std::complex<double>(  1.009 ,0. );
    par_[3]  = std::complex<double>( -0.005 ,0. );
    par_[4]  = std::complex<double>( -0.078 ,0. );
    par_[5]  = std::complex<double>(  0.000 ,0. );
    par_[6]  = std::complex<double>(  0.001 ,0. );
    par_[7]  = std::complex<double>( -0.304 ,0. );
    par_[8]  = std::complex<double>( -0.167 ,0. );
    par_[9]  = std::complex<double>(  3.211 ,0. );
    par_[10] = std::complex<double>( -4.103 ,0. );
  }
  else if (param == Wilson::varyC9C10){
    par_[1]  = std::complex<double>( -0.257 ,0. );
    par_[2]  = std::complex<double>(  1.009 ,0. );
    par_[3]  = std::complex<double>( -0.005 ,0. );
    par_[4]  = std::complex<double>( -0.078 ,0. );
    par_[5]  = std::complex<double>(  0.000 ,0. );
    par_[6]  = std::complex<double>(  0.001 ,0. );
    par_[7]  = std::complex<double>( -0.304 ,0. );
    par_[8]  = std::complex<double>( -0.167 ,0. );
    par_[9]  = std::complex<double>(  3.511 ,0. );
    par_[10] = std::complex<double>( -3.403 ,0. );
  }
  else{
    par_[1]  = std::complex<double>( 0.,0. );
    par_[2]  = std::complex<double>( 0.,0. );
    par_[3]  = std::complex<double>( 0.,0. );
    par_[4]  = std::complex<double>( 0.,0. );
    par_[5]  = std::complex<double>( 0.,0. );
    par_[6]  = std::complex<double>( 0.,0. );
    par_[7]  = std::complex<double>( 0.,0. );
    par_[8]  = std::complex<double>( 0.,0. );
    par_[9]  = std::complex<double>( 0.,0. );
    par_[10] = std::complex<double>( 0.,0. );
  }  
}

Wilson::~Wilson() {} 

std::complex<double> Wilson::operator()( const unsigned int i ) const
{
  assert( i > 0 && i < 11 );

  return par_[i];
}

std::complex<double>& Wilson::operator[]( const unsigned int i ) 
{
  assert( i > 0 && i < 11 );

  return par_[i];
}

Wilson& Wilson::operator=(const Wilson& other){

	for ( unsigned int i = 0; i < other.par_.size(); ++i ){ 
	    
	    par_[i] = other.par_[i];
	}

	return *this;
}

void Wilson::print() const { 
  std::cout
    << "  C7 = " << par_[7] 
    << "  C9 = " << par_[9] 
    << " C10 = " << par_[10]
    << std::endl;
  return ;
} 

std::complex< double > Wilson::h0( const double s ) const
{
    double reh = (8./27.) + (8./9.)*TMath::Log(Parameters::mu) - (4./9.)*TMath::Log(s);
    
    double imh = (4./9.)*TMath::Pi() ;
  
    return std::complex<double>(reh,imh);
}

std::complex<double> Wilson::h( const double s, const double m )  const 
{
    const double z = 4*m*m/s;
    
    double reh = -(4./9.)*( TMath::Log(std::pow(m/Parameters::mu,2)) - (2./3.) - z );
    double imh = 0.0;

    double multiplier = -(4./9.)*(2. + z)*sqrt( fabs(z-1) );
  
    if ( z > 1.0 ){
        reh = reh + multiplier*atan(1.0/(sqrt(z-1)));
    }
    else {
        reh = reh + multiplier*TMath::Log((1. + sqrt( 1.0-z ))/sqrt(z));
        imh = -multiplier*TMath::PiOver2();
    }
  
  return std::complex<double>( reh, imh );
}

std::complex<double> Wilson::Y( const double s ) const 
{
  std::complex<double> calcY(0,0);

  calcY += 1.0*h(s,Parameters::mc)*( (4./3.)*par_[1] + par_[2] + 6.*par_[3] + 60.*par_[5] );
  
  calcY -= 0.5*h(s,Parameters::mb)*( (7.)*par_[3] + (4./3.)*par_[4] + (76.)*par_[5] + (64./3.)*par_[6] );
  
  calcY -= 0.5*h0(s)*( par_[3] + (4./3.)*par_[4] + (16.)*par_[5] + (64./3.)*par_[6] );
  
  calcY += (4./3.)*par_[3] + (64./9.)*par_[5] + (64./27.)*par_[6];
  
  return calcY;
}


std::complex<double> Wilson::Ysum( const double s, const std::vector< ResonantState >& res ) const {
  
  std::complex<double> result(0,0);
  
  for ( const auto & r: res ){ 
    result += resonance( s, r.mass, r.width, r.arg, r.scale );
  }
  
  return result ;
}


std::complex<double> Wilson::resonance( const double s,
					const double mean,
					const double gamma,
					const double arg,
					const double scale ) const
{
  double q = std::sqrt(s);
  
  double cotphi =  ( mean * mean - q * q ) /  ( mean * gamma ) ;
   
  static const std::complex<double> i( 0 , 1 );

  std::complex<double> amp = scale/( cotphi - i );

  std::complex<double> rotate( std::cos( arg ), std::sin( arg ) );

  return rotate*amp;  
}

ResonantState::ResonantState( const double mpole,
			      const double gamma,
			      const double phase,
			      const double magnitude ) : 
  mass( mpole ) , width( gamma ), arg( phase ), scale( magnitude ) {}


ResonantState::ResonantState( const ResonantState& other ) : 
  mass ( other.mass  ),
  width( other.width ),
  arg  ( other.arg   ), 
  scale( other.scale ) {} 
