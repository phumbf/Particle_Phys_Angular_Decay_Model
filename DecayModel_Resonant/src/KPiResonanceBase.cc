#include "KPiResonanceBase.h"
#include "Wignerd.h"

#include "HelperFunctions.h"
#include "Variables.h"

#include <cmath>
#include <iostream>

//The base class for the resonances, containing information 
//relevant for all the different resonances

KPiResonanceBase::KPiResonanceBase( const double mB,
                                    const double mpsi,
                                    const double mkaon,
                                    const double mpion ) :
    mB_     ( mB ) ,
    mpsi_   ( mpsi ),
    mkaon_  ( mkaon ),
    mpion_  ( mpion ),
    min_    ( mpion + mkaon ),
    max_    ( mB - mpsi ),
    polar_  ( true ), 
    offset_ ( 0 ),
    A0_     ( 0, 0 ),
    Ap_     ( 0, 0 ), 
    Am_     ( 0, 0 ),
    A0CP_   ( 0, 0 ),
    ApCP_   ( 0, 0 ), 
    AmCP_   ( 0, 0 ){} ; 


KPiResonanceBase::KPiResonanceBase( const KPiResonanceBase& other ) :
    mB_     ( other.mB_ ) ,
    mpsi_   ( other.mpsi_ ),
    mkaon_  ( other.mkaon_ ),
    mpion_  ( other.mpion_ ),
    min_    ( other.mpion_ + other.mkaon_ ),
    max_    ( other.mB_ - other.mpsi_ ),
    polar_  ( other.polar_ ),
    offset_ ( 0 ),
    A0_     ( other.A0_ ),
    Ap_     ( other.Ap_ ),
    Am_     ( other.Am_ ),
    A0CP_   ( other.A0CP_ ),
    ApCP_   ( other.ApCP_ ),
    AmCP_   ( other.AmCP_ ) {} ;


std::ostream& operator<<( std::ostream& stream, const KPiResonanceBase& resonance ){
    
    unsigned int j = resonance.spin();
    
    stream
        << " ==> INFO (KPiResonanceBase): Resonance with J = "
        << j
        << "\n";
    
    if ( resonance.isNamed() ) {
        stream
            << "     Name = "
            << resonance.name()
            <<  "\n";
    }
    
    stream
      << "       A0 = " << resonance.A0_ << "\n";
    
    if ( j > 0 ){
        stream
	  << "       A+ = " << resonance.Ap_  << "\n" 
	  << "       A- = " << resonance.Am_  << "\n" ;
    }

    stream << "     Type = " 
	   << resonance.getType() << "\n" ;
        
    return stream;
}

void KPiResonanceBase::setAmplitudes( const std::complex<double>& A0, 
				      const std::complex<double>& Ap, 
				      const std::complex<double>& Am ){
  A0_ = A0;
  Ap_ = Ap;
  Am_ = Am;
  
  Helper::cpBasisSwap( A0_, Ap_, Am_, A0CP_, ApCP_, AmCP_ , spin() );
}

void KPiResonanceBase::setAmplitudesPolar( const double magA0, const double phaseA0 ){
  
  A0_ = std::polar( magA0, phaseA0 ); 
  
  if ( spin() > 0 ){
    std::cout << " ==> WARNING (KPiResonanceBase) : J > 0, but only setting A0 " << std::endl;
  }
  
  Helper::cpBasisSwap( A0_, Ap_, Am_, A0CP_, ApCP_, AmCP_ , spin() );
  
  return ;
}

void KPiResonanceBase::setAmplitudesPolar( const double magA0, const double magAp, const double magAm,
					   const double phaseA0, const double phaseAp, const double phaseAm ){
    
  A0_ = std::polar( magA0, phaseA0 );
  
  if ( spin() > 0 ){
    Ap_ = std::polar( magAp, phaseAp );
    Am_ = std::polar( magAm, phaseAm ); 
  }
  
  Helper::cpBasisSwap( A0_, Ap_, Am_, A0CP_, ApCP_, AmCP_ , spin() );
  
  return ;
}

void KPiResonanceBase::setAmplitudesCartesian( const double reA0, const double imA0 ){
  
  A0_ = std::complex<double>( reA0, imA0 );
  
  if ( spin() > 0 ){
    std::cout << " ==> WARNING (KPiResonanceBase) : J > 0, but only setting A0 " << std::endl;
  }
    
  Helper::cpBasisSwap( A0_, Ap_, Am_, A0CP_, ApCP_, AmCP_ , spin() );
  
  return ;
}


void KPiResonanceBase::setAmplitudesCartesian( const double reA0, const double imA0, 
					       const double reAp, const double imAp, 
					       const double reAm, const double imAm ){
  
  A0_ = std::complex<double>( reA0, imA0 );
  
  if ( spin() > 0 ){ 
    Ap_ = std::complex<double>( reAp, imAp );
    Am_ = std::complex<double>( reAm, imAm );
  }
  
  Helper::cpBasisSwap( A0_, Ap_, Am_, A0CP_, ApCP_, AmCP_ , spin() );
  
  return ; 
}




std::complex<double> KPiResonanceBase::evaluate( const Event& event, const int lambda, const bool isCP ) const {

  static std::complex<double> i(0,1);
  
  const unsigned int j = spin();
  
  std::complex<double> result(0,0);

  result += ( ( isCP ? A0CP_ : A0_ )*
	      Wignerd::djmpm( event.cospsi_, 1, lambda, 0 )*
	      Wignerd::djmpm( event.coskpi_, j, 0, 0 ) );
  
  if ( j > 0 ){
    result += ( ( isCP ? ApCP_ : Ap_ )*
		Wignerd::djmpm( event.cospsi_, 1, lambda, 1 )*
		Wignerd::djmpm( event.coskpi_, j, 0, 1 )*
		std::exp( -i*event.phi_ ) );
    
    result += ( ( isCP ? AmCP_ : Am_ )*
		Wignerd::djmpm( event.cospsi_, 1, lambda, -1 )*
		Wignerd::djmpm( event.coskpi_, j, 0, -1 )*
		std::exp( i*event.phi_ ) );
  }
  result *= lineshape( event.mkpi_, event.pB_, event.pR_ ) ;
  return result;
}

std::complex<double> KPiResonanceBase::evaluate( const Event& event, const int lambda) const {

  static std::complex<double> i(0,1);
  
  const unsigned int j = spin();
  
  std::complex<double> result(0,0);
  
  result += ( A0_ *
	      Wignerd::djmpm( event.cospsi_, 1, lambda, 0 )*
	      Wignerd::djmpm( event.coskpi_, j, 0, 0 ) );
  
  if ( j > 0 ){
    result += ( Ap_ *
		Wignerd::djmpm( event.cospsi_, 1, lambda, 1 )*
		Wignerd::djmpm( event.coskpi_, j, 0, 1 )*
		std::exp( -i*event.phi_ ) );
    
    result += ( Am_ *
		Wignerd::djmpm( event.cospsi_, 1, lambda, -1 )*
		Wignerd::djmpm( event.coskpi_, j, 0, -1 )*
		std::exp( i*event.phi_ ) );
  }
  
  result *= lineshape( event.mkpi_, event.pB_, event.pR_ ) ;

  return result;
}

double KPiResonanceBase::project( const double m ) const {
    
    const double pB = Helper::momentum( mB_, m, mpsi_ );
    const double pR = Helper::momentum( m, mkaon_, mpion_ );
    
    std::complex<double> shape = lineshape( m, pB, pR  );
    
    return getNorm()*std::norm( shape );
}


std::complex<double> KPiResonanceBase::evaluate( const Event& event, const int lambdapsi, const int lambda, const bool isCP ) const {
    
    const unsigned int j = spin();
    
    static std::complex<double> i(0,1);
    
    std::complex<double> result(0,0);
    
    switch ( lambdapsi ){
    case -1:
      result = ( isCP ? AmCP_ : Am_ );
      break;
    case  0:
      result = ( isCP ? A0CP_ : A0_ );
      break;
    case  1:
      result = ( isCP ? ApCP_ : Ap_ );
      break;
    default:
      break;
    }
    
    // angular part
    result *= ( Wignerd::djmpm( event.coskpi_, 1, lambda, lambdapsi )*
		Wignerd::djmpm( event.cospsi_, j, 0, lambdapsi )*
		std::exp( -event.phi_*lambdapsi*i )  );
    
    // mass lineshape
    result *= lineshape( event.mkpi_, event.pB_, event.pR_ ) ;
   
    return result;
}


std::string KPiResonanceBase::name() const {
    return name_ ;
}

bool KPiResonanceBase::isNamed() const {
    return ( name_.length() > 0 );
}

bool KPiResonanceBase::isNamed( const std::string& name ) const {
    return ( name_ == name );
}

bool KPiResonanceBase::isNameIn( const std::string& name ) const { 
  return ( name.find( name_ ) != std::string::npos );
}

void KPiResonanceBase::setName( const std::string& name ){
    name_ = name;
}


void KPiResonanceBase::usePolar() { 
  polar_ = true;
  return ;
}

void KPiResonanceBase::useCartesian() { 
  polar_ = false;
  return ;
}


bool KPiResonanceBase::forbidden( const double m ) const { 
  return ( ( m < ( mpion_ + mkaon_ ) ) || 
	   ( m > ( mB_    - mpsi_  ) ) );
}


std::complex<double> KPiResonanceBase::getA0() const { 
  return A0_;
}

std::complex<double> KPiResonanceBase::getAp() const { 
  return Ap_;
}

std::complex<double> KPiResonanceBase::getAm() const { 
  return Am_;
}

double KPiResonanceBase::getNorm() const {
  return ( std::norm( A0_ ) + 
	   std::norm( Ap_ ) + 
	   std::norm( Am_ ) );
}


double KPiResonanceBase::getA0Fraction() const { 
  return std::norm( A0_ )/getNorm(); 
}

double KPiResonanceBase::getApFraction() const { 
  return std::norm( Ap_ )/getNorm(); 
}

double KPiResonanceBase::getAmFraction() const { 
  return std::norm( Am_ )/getNorm(); 
}

double KPiResonanceBase::getAparaFraction() const { 
  
  std::complex<double> Azero, Apara, Aperp; 
  
  Helper::convertToTransversityBasis( A0_, Ap_, Am_, Azero, Apara, Aperp );
  
  return std::norm( Apara )/getNorm(); 
}

double KPiResonanceBase::getAperpFraction() const { 
  
  std::complex<double> Azero, Apara, Aperp; 
  
  Helper::convertToTransversityBasis( A0_, Ap_, Am_, Azero, Apara, Aperp );
  
  return std::norm( Aperp )/getNorm(); 
}


void KPiResonanceBase::setParameters( const std::vector<double>& par ) { 
  
  if ( polar_ ){ 
    setAmplitudesPolar( par[ offset_ ], par[ offset_ + 1 ] );
  }
  else { 
    setAmplitudesCartesian( par[ offset_ ], par[ offset_ + 1 ] );
  }
  
  
  if ( spin() > 0 ){
    if ( polar_ ){ 
      setAmplitudesPolar( par[ offset_ ]    , par[ offset_ + 2 ],  par[ offset_ + 4 ], 
			  par[ offset_ + 1 ], par[ offset_ + 2 ],  par[ offset_ + 4 ] );
    }
    else { 
      setAmplitudesCartesian( par[ offset_ ], par[ offset_ + 1 ],  
			      par[ offset_ + 2 ], par[ offset_ + 3 ],  
			      par[ offset_ + 4 ], par[ offset_ + 5 ] );
    }
  }

  return;
}

unsigned KPiResonanceBase::getParameters( ROOT::Minuit2::MnUserParameters& par )   {
  
  offset_ = par.Parameters().size() ;

  if ( polar_ ){ 
    par.Add( name() + "_magA0", std::abs(A0_) , 0.1, 0, 100*std::abs(A0_) );
    par.Add( name() + "_argA0", std::arg(A0_) , 0.1, -TMath::Pi(), TMath::Pi() );
  } else { 
    par.Add( name() + "_reA0", std::real(A0_) , 0.1, 0, 100*std::abs(A0_) );
    par.Add( name() + "_imA0", std::imag(A0_) , 0.1, 0, 100*std::abs(A0_));
  }
  
  if ( spin() == 0 ) return 2;
  
  if ( polar_ ) { 
    par.Add( name() + "_magAp", std::abs(Ap_) , 0.1, 0, 100*std::abs(Ap_) );
    par.Add( name() + "_argAp", std::arg(Ap_) , 0.1, -TMath::Pi(), TMath::Pi() );
    par.Add( name() + "_magAm", std::abs(Am_) , 0.1, 0, 100*std::abs(Am_) );
    par.Add( name() + "_argAm", std::arg(Am_) , 0.1 , -TMath::Pi(), TMath::Pi() );
  } else { 
    par.Add( name() + "_reAp", std::real(Ap_) , 0.1, 0, 100*std::abs(Ap_) );
    par.Add( name() + "_imAp", std::imag(Ap_) , 0.1, 0, 100*std::abs(Ap_));
    par.Add( name() + "_reAm", std::real(Am_) , 0.1, 0, 100*std::abs(Am_) );
    par.Add( name() + "_imAm", std::imag(Am_) , 0.1, 0, 100*std::abs(Am_));
  }

  return 6;
}

void KPiResonanceBase::Rescale(const double scalefactor){


	A0_ *= scalefactor;
	A0CP_ *= scalefactor;
	Ap_ *= scalefactor;
	ApCP_ *= scalefactor;
	Am_ *= scalefactor;
	AmCP_ *= scalefactor;


	return;
}
