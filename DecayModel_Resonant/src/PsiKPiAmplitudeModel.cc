/*
 * Amplitude model for B --> psi K pi, including
 * sum of conventional (K pi) resonances and exotic
 * (psi pi) states. 
 *
 * Model is configured using an external json file in 
 * PsiKPIAmplitudeModel::initialise 
 * 
 */


// Class header
#include "PsiKPiAmplitudeModel.h"

// std 
#include <complex>
#include <cmath>

// Local
#include "KPiResonanceBase.h"
#include "PsiPiResonance.h" 
#include "Variables.h"
#include "Parser.h"

PsiKPiAmplitudeModel::PsiKPiAmplitudeModel() {} ;

PsiKPiAmplitudeModel::PsiKPiAmplitudeModel( const std::string& json ){
  initialise( json );
} ;


PsiKPiAmplitudeModel::PsiKPiAmplitudeModel( const PsiKPiAmplitudeModel& other ){
  
  std::vector< KPiResonanceBase* >::const_iterator it;
  std::vector< PsiPiResonance*   >::const_iterator jt ; 
  
  for ( it  = other.kpi_resonances_.begin();
	it != other.kpi_resonances_.end(); ++it ) {
    kpi_resonances_.push_back( (*it)->clone() );
  }

  for ( jt  = other.psipi_resonances_.begin();
	jt != other.psipi_resonances_.end(); ++jt ) {
    psipi_resonances_.push_back( (*jt)->clone() );
  }
}


PsiKPiAmplitudeModel::~PsiKPiAmplitudeModel() {
    
  std::vector< KPiResonanceBase* >::iterator it;
  std::vector< PsiPiResonance*   >::iterator jt ; 

  for ( it  = kpi_resonances_.begin();
	it != kpi_resonances_.end(); ++it ) {
    delete (*it);
  }
  
  for ( jt  = psipi_resonances_.begin();
	jt != psipi_resonances_.end(); ++jt ) {
    delete (*jt);
  }
}



void PsiKPiAmplitudeModel::clearResonanceList() { 
  /*
   * Clear the list of resonances
   */ 
  
  std::vector< KPiResonanceBase* >::iterator it;
  std::vector< PsiPiResonance*   >::iterator jt ; 

  for ( it  = kpi_resonances_.begin();
	it != kpi_resonances_.end(); ++it ) {
    delete (*it);
  }
  
  for ( jt  = psipi_resonances_.begin();
	jt != psipi_resonances_.end(); ++jt ) {
    delete (*jt);
  }

  kpi_resonances_.clear();
  psipi_resonances_.clear(); 
  
  return;
}


void PsiKPiAmplitudeModel::setParameters( const std::vector<double>& par ) {
  /* 
   * Propogate fit parameters to the model 
   */
  std::vector< KPiResonanceBase* >::iterator it;
  
  for ( it  = kpi_resonances_.begin();
	it != kpi_resonances_.end(); ++it ) {
    (*it)->setParameters( par ); 
  }
  
  return ;
}

unsigned PsiKPiAmplitudeModel::getParameters( ROOT::Minuit2::MnUserParameters& par ){ 
  /* 
   * Get set of fit parameters for the model 
   */
  
  std::vector< KPiResonanceBase* >::iterator it;
  
  unsigned int nparam = 0;
  
  for ( it  = kpi_resonances_.begin();
	it != kpi_resonances_.end(); ++it ) {
    nparam += (*it)->getParameters( par ); 
  }
  
  return nparam; 
}
 

double PsiKPiAmplitudeModel::operator()( const double* x ) const {
  /* 
   * operator for use in event generation and fitting 
   */
  
  Event evt( 511, x[0], x[1], x[2], x[3] );
  
  return evaluate( evt );
}

void PsiKPiAmplitudeModel::initialise( const std::string& json ){
  /*
   *  Read list of resonances from JSON input file
   */
  clearResonanceList(); 
    
  // Parse resonances from JSON file (see Parser.{h,cc})
 
  parse_resonance_list( json, kpi_resonances_ , psipi_resonances_, false );

  return ;
}

double PsiKPiAmplitudeModel::evaluate( const Event& event ) const {
  double result = 0;
  

  // loop over lambda = -1, +1 for the psi decay
  for ( int lambda = -1; lambda <= 1; lambda += 2 ){
    result += std::norm( getAmplitude( lambda, event ) );
  }

  // multiplied by Kpi momentum in B-frame and K momentum in Kpi frame
  double kinematic = event.pB_*event.pR_;
  if(kpi_resonances_.empty() && psipi_resonances_.empty()){
	std::cout << "Kinematic is: " << kinematic << std::endl;
	return kinematic;
  }

  result *= kinematic;
  return result;
}

double PsiKPiAmplitudeModel::evaluateCP( const Event& event ) const {
  
  double result = 0;

  // loop over lambda = -1, +1 for the psi decay
  for ( int lambda = -1; lambda <= 1; lambda += 2 ){

   std::complex<double> A    = getAmplitude( lambda, event, false ) ;
   std::complex<double> Abar = getAmplitude( lambda, event, true  ) ;
    //add coherently |A + Abar|^2 = |A|^2 + |Abar|^2 + crossterms.
    result += std::norm( A );
    result += std::norm( Abar );
    double D = Event::D;
    std::complex<double> qoverp = std::polar(1.0,-0.763);
    result -= 2.*D*std::real( qoverp * std::conj( A ) * Abar );
  }
  
  // multiplied by Kpi momentum in B-frame and K momentum in Kpi frame
  double kinematic = event.pB_*event.pR_;
  
  if(kpi_resonances_.empty() && psipi_resonances_.empty()){
	return kinematic;
  }
  result *= kinematic;
  
  return result;
}

void PsiKPiAmplitudeModel::print(){
	
 std::vector< KPiResonanceBase* >::const_iterator it;
  
  for ( it = kpi_resonances_.begin(); it != kpi_resonances_.end(); ++it ){
    std::cout << **it << std::endl;
  }
}


double PsiKPiAmplitudeModel::evaluateCP( const std::string& name, const Event& event ) const {
  
  double result = 0;

  // loop over lambda = -1, +1 for the psi decay
  for ( int lambda = -1; lambda <= 1; lambda += 2 ){

   std::complex<double> A    = getAmplitude( name,lambda, event, false ) ;
   std::complex<double> Abar = getAmplitude( name,lambda, event, true  ) ;
    //add coherently |A + Abar|^2 = |A|^2 + |Abar|^2 + crossterms.
    result += std::norm( A );
    result += std::norm( Abar );
    double D = 0.01;
    std::complex<double> qoverp = std::polar(1.0,-0.763);
    result -= 2.*D*std::real( qoverp * std::conj( A ) * Abar );
  }
  
  // multiplied by Kpi momentum in B-frame and K momentum in Kpi frame
  result *= event.pB_*event.pR_;
  
  return result;
}
double PsiKPiAmplitudeModel::evaluate( const std::string& name, const Event& event ) const {
    
  /*
   *  Evaluate model for resonances appearing in name 
   *  e.g. if name = "K*(892),K0(800)" only the
   *  "K*(892)" and the "K0(800)" will be considered
   */

  double result = 0;
  
  // loop over lambda = -1, +1 for the psi decay
  for ( int lambda = -1; lambda <= 1; lambda += 2 ){
    result += std::norm( getAmplitude( name, lambda, event ) ); 
  }
  // multiplied by Kpi momentum in B-frame and K momentum in Kpi frame
  result *= event.pB_*event.pR_;
  
  return result;
}


double PsiKPiAmplitudeModel::evaluate( const unsigned int j, const Event& event ) const {
    
  double result = 0;
  
  for ( int lambda = -1; lambda <= 1; lambda += 2 ){
    result += std::norm( getAmplitude( j, lambda, event ) ); 
  }
  
  // multiplied by Kpi momentum in B-frame and K momentum in Kpi frame
  result *= event.pB_*event.pR_;
  
  return result;
}

const KPiResonanceBase* PsiKPiAmplitudeModel::getResonance( const std::string& name ) const { 
  
  std::vector< KPiResonanceBase* >::const_iterator it;
  
  for ( it = kpi_resonances_.begin(); it != kpi_resonances_.end(); ++it ){
    if ( (*it)->isNamed( name ) ) return (*it);
  }
  
  return 0;
}

std::complex<double>  PsiKPiAmplitudeModel::getAmplitude( const int lambda, const Event& event, const bool isCP ) const{

  std::vector< KPiResonanceBase* >::const_iterator it;
  std::vector< PsiPiResonance*   >::const_iterator jt;
  
  std::complex<double> fn(0,0);
  // sum over K pi resonances
  for ( it = kpi_resonances_.begin(); it != kpi_resonances_.end(); ++it ){

    fn += (*it)->evaluate( event, lambda, isCP ); // summing Kpi helicities
  }

  // sum over psi pi resonances
  for ( jt = psipi_resonances_.begin(); jt != psipi_resonances_.end(); ++jt ){ 
    fn += (*jt)->evaluate( event, lambda );
  }
  
  return fn;

}

std::complex<double>  PsiKPiAmplitudeModel::getAmplitude( const std::string& name, const int lambda, const Event& event, const bool isCP ) const{

    std::complex<double> fn(0,0) ;
    
    std::vector< KPiResonanceBase* >::const_iterator it;
    std::vector< PsiPiResonance*   >::const_iterator jt;

    // sum over k pi resonances
    for ( it = kpi_resonances_.begin(); it != kpi_resonances_.end(); ++it ){
      if ( (*it)->isNameIn( name ) ){ 
	fn += (*it)->evaluate( event, lambda, isCP );
      }
    }

    for ( jt = psipi_resonances_.begin(); jt != psipi_resonances_.end(); ++jt ){ 
      if ( (*jt)->isNameIn( name ) ) { 
	fn += (*jt)->evaluate( event, lambda );
      }
    }

    return fn;
}


std::complex<double>  PsiKPiAmplitudeModel::getAmplitude( const unsigned int j, const int lambda, const Event& event, const bool isCP ) const{
    
  std::complex<double> fn(0,0) ;
    
  std::vector< KPiResonanceBase* >::const_iterator it;

  for ( it = kpi_resonances_.begin(); it != kpi_resonances_.end(); ++it ){
    if ( (*it)->spin() == j ) {
      fn += (*it)->evaluate( event, lambda, isCP );
    }
  }

  return fn;
}


//normal evaluate procedure

std::complex<double>  PsiKPiAmplitudeModel::getAmplitude( const int lambda, const Event& event) const{

  
  std::vector< KPiResonanceBase* >::const_iterator it;
  std::vector< PsiPiResonance*   >::const_iterator jt;
  
  std::complex<double> fn(0,0);
  
  // sum over K pi resonances
  
  for ( it = kpi_resonances_.begin(); it != kpi_resonances_.end(); ++it ){
    fn += (*it)->evaluate( event, lambda); // summing Kpi helicities
  }
  
  
  // sum over psi pi resonances
  for ( jt = psipi_resonances_.begin(); jt != psipi_resonances_.end(); ++jt ){ 
    fn += (*jt)->evaluate( event, lambda );
  }
  
  return fn;

}

std::complex<double>  PsiKPiAmplitudeModel::getAmplitude( const std::string& name, const int lambda, const Event& event) const{

    std::complex<double> fn(0,0) ;
    
    std::vector< KPiResonanceBase* >::const_iterator it;
    std::vector< PsiPiResonance*   >::const_iterator jt;

    // sum over k pi resonances
    for ( it = kpi_resonances_.begin(); it != kpi_resonances_.end(); ++it ){
      if ( (*it)->isNameIn( name ) ){ 
	fn += (*it)->evaluate( event, lambda);
      }
    }

    for ( jt = psipi_resonances_.begin(); jt != psipi_resonances_.end(); ++jt ){ 
      if ( (*jt)->isNameIn( name ) ) { 
	fn += (*jt)->evaluate( event, lambda );
      }
    }

    return fn;
}


std::complex<double>  PsiKPiAmplitudeModel::getAmplitude( const unsigned int j, const int lambda, const Event& event) const{
    
  std::complex<double> fn(0,0) ;
    
  std::vector< KPiResonanceBase* >::const_iterator it;

  for ( it = kpi_resonances_.begin(); it != kpi_resonances_.end(); ++it ){
    if ( (*it)->spin() == j ) {
      fn += (*it)->evaluate( event, lambda);
    }
  }

  return fn;
}

std::map<std::string,double> PsiKPiAmplitudeModel::fitfractions(const bool isCP){

	std::vector< KPiResonanceBase* >::const_iterator it;
	std::vector< PsiPiResonance*   >::const_iterator jt;
	
	std::map<std::string,double> fracmap;

	const unsigned int nGenerated = 100000; 

	const std::string def{""};
	//loop through resonances and fill fitfractions
	double fitfracdenom = integrate(def,nGenerated,isCP);

	for ( it = kpi_resonances_.begin(); it != kpi_resonances_.end(); ++it ){
		fracmap[(*it)->name()] = (integrate((*it)->name(),nGenerated,isCP)) / fitfracdenom;
	}
	for ( jt = psipi_resonances_.begin(); jt != psipi_resonances_.end(); ++jt ){
		fracmap[(*jt)->name()] = (integrate((*jt)->name(),nGenerated,isCP)) / fitfracdenom;
	}

	return fracmap;
}

void PsiKPiAmplitudeModel::rescaleres(const std::string& refres, std::map<std::string,double>& fitfracs, std::map<std::string,double>& targetfracs){

	//Change to iterators over const iterators
	std::vector< KPiResonanceBase* >::iterator it;
	std::vector< PsiPiResonance*   >::iterator jt;
	
	double scalefactor{0},desired{0},current{0};

	for ( it = kpi_resonances_.begin(); it != kpi_resonances_.end(); ++it ){
		desired = targetfracs[(*it)->name()] / targetfracs[ refres ];
		current = fitfracs[(*it)->name()] / fitfracs[refres];
		scalefactor = sqrt(desired / current);
		(*it)->Rescale(scalefactor);
	}

	return;
}

double PsiKPiAmplitudeModel::integrate(std::string name, const unsigned int numevts, const bool isCP){
	
	double result{0};

	double mkpi{0};
	double cospsi{0};
	double coskpi{0};
	double phi{0};

	double mkpi_max = Event::mB - Event::mPsi;
	double mkpi_min = Event::mKaon + Event::mPion;
	double cos_max = 1.0, cos_min = -1.0, phi_max = TMath::Pi(), phi_min = -TMath::Pi();
	double volume = (mkpi_max-mkpi_min)*(cos_max-cos_min)*(cos_max-cos_min)*(phi_max-phi_min);
	
	for(unsigned int i=0; i<numevts; i++){

		mkpi   = gRandom->Uniform( mkpi_min, mkpi_max);
		coskpi = gRandom->Uniform( cos_min, cos_max );
		cospsi = gRandom->Uniform( cos_min, cos_max );
		phi    = gRandom->Uniform( phi_min, phi_max );

		Event evt(511, mkpi,coskpi,cospsi,phi);

		if(name.empty()){
		
			if(isCP){
				result += this->evaluateCP(evt);
			}
			else{
				result += this->evaluate(evt);
			}
		}

		else{
			if(isCP){
				result += this->evaluateCP(name,evt);
			}
			else{
				result += this->evaluate(name,evt);
			}
		}
	}
	

	double toRet = (volume/numevts)*result;
	return toRet;	
}

std::vector<std::string> PsiKPiAmplitudeModel::resonancelist(){

	std::vector< KPiResonanceBase* >::iterator it;
	std::vector< PsiPiResonance*   >::iterator jt;

	std::vector<std::string> resvector;

	for ( it = kpi_resonances_.begin(); it != kpi_resonances_.end(); ++it ){
		resvector.push_back((*it)->name());
	}
	for ( jt = psipi_resonances_.begin(); jt != psipi_resonances_.end(); ++jt ){
		resvector.push_back((*jt)->name());

	}
	return resvector;

}
