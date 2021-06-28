// Local 

#include "ResonanceBase.h"
#include "FFparser.h"
#include "Wilson.h"
#include "Parameters.h"
#include "Kinematic.h"


//The base class for a resonant object
double ResonanceBase::getMotherMass() const{
  return m_mothermass;
}

double ResonanceBase::getMass() const{
  return m_mass;
}

void ResonanceBase::setMotherMass( const double mothermass ){
  m_mothermass = mothermass;
}

void ResonanceBase::setMass( const double mass ){
  m_mass = mass;
}

void ResonanceBase::setBFfac( const double BFfac){

	m_branching = BFfac;
}

void ResonanceBase::setisCP( const bool isCP){
	m_isCP = isCP;
}

bool ResonanceBase::isCP(){
	if(m_isCP){
		return true;
	}
	else{
		return false;
	}
}

void ResonanceBase::setisEvtGenWC( const bool isEvtGenWC){
    m_isEvtGenWC = isEvtGenWC;
}

bool ResonanceBase::isEvtGenWC(){

    if(m_isEvtGenWC){
		return true;
	}
	else{
		return false;
	}
}

void ResonanceBase::setisbtod( const bool isbtod){
    m_isbtod = isbtod;
}

bool ResonanceBase::isbtod(){

    if(m_isbtod){
		return true;
	}
	else{
		return false;
	}
}

int ResonanceBase::getL() const{
  return m_L;
}

void ResonanceBase::setL( const int l ){
  m_L = l;
}

std::string ResonanceBase::getName() const{
  return m_name;
}

void ResonanceBase::setName( const std::string name ){
  m_name = name;
}

bool ResonanceBase::hasName( const std::string name ) const { 
  return ( m_name == name );
}

void ResonanceBase::setDaughterMasses( const double mone, const double mtwo ){
  m_mone = mone;
  m_mtwo = mtwo;
}

void ResonanceBase::print() const{
  std::cout 
    << "Resonance: " << getName() << " (L = " << getL() << ")\n " 
    << "\t Mother mass = " << getMotherMass() 
    << std::endl;

}

void ResonanceBase::setShape( std::unique_ptr<HHLineShapeBase>& shape ){

  m_lineshape = std::move(shape);
}

void ResonanceBase::getAmplitudes( const double s   ,
				   const double mhh ,
				   AmplitudeSet& amp) const {


  const std::complex<double> lineshape = m_lineshape->lineshape( mhh );

  // add in the SphericalHarmonics
  amp.AZEROL = AZERO( s, mhh,  1 )*lineshape; 
  amp.AZEROR = AZERO( s, mhh, -1 )*lineshape; 
  amp.APARAL = APARA( s, mhh,  1 )*lineshape; 
  amp.APARAR = APARA( s, mhh, -1 )*lineshape; 
  amp.APERPL = APERP( s, mhh,  1 )*lineshape; 
  amp.APERPR = APERP( s, mhh, -1 )*lineshape; 
  amp.ATIME  = ATIME( s, mhh )*lineshape; 

}

void ResonanceBase::getAmplitudesCPswap( const double s,
                                         const double mhh,
                                         AmplitudeSet& amp) const{


    const std::complex<double> lineshape = m_lineshape->lineshape( mhh );

    // add in the SphericalHarmonics
    amp.AZEROL = AZERO( s, mhh,  1 )*lineshape; 
    amp.AZEROR = AZERO( s, mhh, -1 )*lineshape; 
    amp.APARAL = APARA( s, mhh,  1 )*lineshape; 
    amp.APARAR = APARA( s, mhh, -1 )*lineshape; 
    amp.APERPL = APERP( s, mhh,  1 )*lineshape; 
    amp.APERPR = APERP( s, mhh, -1 )*lineshape; 
    amp.ATIME  = ATIME( s, mhh )*lineshape; 

    //swap signs of amplitudes depending upon the spin of the resonance
    if(this->getL() % 2 == 0){

		amp.AZEROR = -amp.AZEROR;
		amp.AZEROL = -amp.AZEROL;
		amp.APARAR = -amp.APARAR;
		amp.APARAL = -amp.APARAL;
	}
	else{
		amp.APERPL = -amp.APERPL;
		amp.APERPR = -amp.APERPR;
	}

}

double ResonanceBase::normalisation( const double s, const double m ) const { 
  
  const double beta   = Kinematic::beta( s, Parameters::ml );
  const double lambda = Kinematic::lambda( s, m, m_mothermass );
  
  double Nsq  = s*std::sqrt( lambda )*beta;
  
  return std::sqrt( Nsq );                                                                                                                
}


void ResonanceBase::setCorrections( const std::complex<double>& cZERO, 
				    const std::complex<double>& cPARA,
				    const std::complex<double>& cPERP ) {
  
  m_cPERP = cPERP; 
  m_cPARA = cPARA;
  m_cZERO = cZERO;
  
  return;
}

void ResonanceBase::setCorrections( const std::complex<double>& cZERO ){
  m_cZERO = cZERO;
}

void ResonanceBase::setWilsons(Wilson &wc, 
		                 Wilson &wp){


	m_WC = wc;
	m_WP = wp;
}

void ResonanceBase::setPolarResFactor(const double alpha,
				 const double phi){
	m_resfactor = std::polar(alpha,phi);
}
void ResonanceBase::setResFactor(const double alpha,
				 const double phi){

    std::complex<double> num(alpha,phi);
    m_resfactor = num;
}
double ResonanceBase::prefactor( const int mproj,
			         const double mass,
				 const double lambdaM,
				 const double lambdaH) const{

  const int l = getL(); 
  const int m = std::abs( mproj );

  //lambdaH normalized term
  const double normalizedterm = Kinematic::lambda(getMass()*getMass(),m_mone,m_mtwo);
  //const double normalizedterm = 1.0;

  //Factors of mass: come from P = lambda^1/2 argument:
  //i.e. in jpsi code we had a P/m term so if have 
  //l =0 -> have mass^1/2, l=1 -> have mass^3/2 etc.`
  if( l == 0){
	return std::sqrt(m_branching*std::sqrt(lambdaH/normalizedterm)/(mass/getMass()));
  }
  else if (l == 1){
	return std::sqrt(m_branching*std::sqrt(std::pow(lambdaH/normalizedterm,3))/std::pow(mass/getMass(),3));
  }
  else if ( l == 2 ){ 
    if ( m == 0 ){
	return std::sqrt( m_branching*lambdaM*std::sqrt(std::pow(lambdaH/normalizedterm,5)) / 6. /std::pow(mass/getMass(),5)) / (m_mothermass*mass);
    }
    if ( m == 1 ){
        return std::sqrt( m_branching*lambdaM*std::sqrt(std::pow(lambdaH/normalizedterm,5)) / 8. /std::pow(mass/getMass(),5)) / (m_mothermass*mass);    
    }
  }
  else{
	std::cout << "ERROR: Spin of resonance not supported" << std::endl;
	return 0;
  }

}
		                 
