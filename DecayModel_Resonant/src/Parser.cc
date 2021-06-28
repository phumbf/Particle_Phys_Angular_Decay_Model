#include "Parser.h"

#include "Variables.h"
#include "KPiResonanceBase.h"
#include "KPiResonanceBW.h"
#include "KPiResonanceBWNorm.h"
#include "KPiResonanceFlatte.h"
#include "KPiResonanceLASS.h" 
#include "PsiPiResonance.h" 
#include "HelperFunctions.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <string>
#include <vector>


//Class to Parse in resonances from a json list
using boost::property_tree::ptree;


std::string string_from_tree( const std::string& s ) {
    return "\"" + s + "\"";
}

template <typename T> std::vector<T>
    list_from_tree( const ptree& pt, const std::string& key ){
    
    std::vector<T> result;
    
    for ( const auto& i : pt.get_child( key ) ){
        result.push_back( i.second.get_value<T>() );
    }
    
    return result;
}



void print_tree( const ptree& pt, int level )
{
    const std::string sep( 2 * level, ' ' );
    
    for ( const auto& i : pt) {
        std::cout
            << sep
            << string_from_tree( i.first ) 
	    << " : " 
	    << string_from_tree( i.second.data() ) 
	    << "\n";
        
        print_tree( i.second, level + 1 );
    }
    return;
}


std::map<std::string,double> get_target_fractions( const std::string json){ 

  std::map<std::string,double> fraction;

  ptree pt;
  
  read_json( json ,  pt );
  
  for ( const auto& i : pt.get_child("resonance") ){
    boost::optional< double > f =  i.second.get_optional<double>("targetfitfrac") ;
    
    if ( f ) {
      fraction[ i.first ] = f.get();
    }
  }
  
  return fraction;
}


void update_json( const std::string json, 
		      std::map< std::string, double >& scalemap ){
  
    ptree pt;
   
    read_json( json ,  pt );
    
    for ( auto& i : pt.get_child("resonance") ){
      if ( scalemap.find( i.first ) != scalemap.end() ) {
	i.second.put< double >("Scaler",scalemap[i.first]);
      }
    }

    write_json( std::cout, pt );
} 
 
void parse_resonance( const std::string& name,
                      std::vector< KPiResonanceBase* >& resonances,
                      const ptree& pt,
		      const std::string& mode ){

  //Read from the Json file and create the resonance base
  std::cout << "Debug: mode is : " << mode << std::endl;
  
  std::vector<double> amp = list_from_tree<double>( pt, "amplitude" );
  int LR = pt.get<int>("LR");
  const std::string basis = pt.get<std::string>("basis");
  const std::string resonance = pt.get<std::string>("type");
  KPiResonanceBase* res;
  
  std::complex<double> H0(0,0);
  std::complex<double> Hp(0,0);
  std::complex<double> Hm(0,0);
  
  if(mode == "cartesian"){
    
    if(LR == 0){
      H0 = std::complex<double>(amp[0],amp[1]);
    }
    else{
      if(basis == "transversity"){
	std::cout << "Debug: HERE " << std::endl;
	std::complex<double> Azero(amp[0],amp[3]);
	std::complex<double> Apara(amp[1],amp[4]);
	std::complex<double> Aperp(amp[2],amp[5]);
	std::cout << "Debug: cart,transversity: A0,Apara,Aperp" << Azero << " " << Apara << " " << Aperp << std::endl;
	Helper::convertToHelicityBasis(Azero,Apara,Aperp,H0,Hp,Hm);		    
	std::cout << "Debug: cart,transversity: H0,Hp,Hm" << H0 << " " << Hp << " " << Hm << std::endl;
	
      }
      else{
	H0 = std::complex<double>(amp[0],amp[3]);
	Hp = std::complex<double>(amp[1],amp[4]);
	Hm = std::complex<double>(amp[2],amp[5]);
      }
    }
  }
  else if(mode == "polar"){
    
    if( LR == 0 ){
      H0 = std::polar(amp[0],amp[1]);
    }
    else{
      
      if(basis == "transversity"){
	std::complex<double> Azero = std::polar(amp[0],amp[3]);
	std::complex<double> Apara = std::polar(amp[1],amp[4]);
	std::complex<double> Aperp = std::polar(amp[2],amp[5]);
	Helper::convertToHelicityBasis(Azero,Apara,Aperp,H0,Hp,Hm);		    
      }
      else{
	H0 = std::polar(amp[0],amp[3]);
	Hp = std::polar(amp[1],amp[4]);
	Hm = std::polar(amp[2],amp[5]);
      }
    }
  }
  else if(mode == "fractions"){
    
    if(LR == 0){
      H0 = std::polar(amp[0],amp[1]);
    }
    
    else{
      std::vector<double> transfracs = list_from_tree<double>( pt, "transversityfracs" );
      
      //Scale with the  transversity fractions.
      std::complex<double> Azero = std::polar(amp[0],amp[3]);
      std::complex<double> Apara = std::polar(amp[1],amp[4]);
      std::complex<double> Aperp = std::polar(amp[2],amp[5]);
      
      std::cout << "Debug: Azero Apara Aperp " << Azero << " " << Apara << " " << Aperp << std::endl;
      
      Helper::transFracScale(Azero,Apara,Aperp,transfracs[1],transfracs[2]);
      
      std::cout << "Debug: for reson: A0,Ap,Aper " << name << " " << Azero<< " " <<Apara << " " <<Aperp << std::endl;
      
      Helper::convertToHelicityBasis(Azero,Apara,Aperp,H0,Hp,Hm);
      
    }
  }
  
  if ( resonance == "BW" || resonance == "BWNorm" ){
    
    double mass   = pt.get<double>("mean");
    double width  = pt.get<double>("width");
    
    boost::optional<double> opt_radius =
      pt.get_optional<double>("radius");
    
    double radius = ( opt_radius ) ? opt_radius.get() : 1.6 ;
    
    boost::optional<int> opt_LB = pt.get_optional<int>("LB");
        
    int LB  = ( opt_LB ) ? opt_LB.get() : abs( LR - 1 );
    
    if ( resonance == "BWNorm" ){
      res = new KPiResonanceBWNorm( mass, width, Event::mB, Event::mPsi, Event::mKaon,  Event::mPion, LR, LB,  radius );  
    }
    else {
      res = new KPiResonanceBW( mass, width,  Event::mB,  Event::mPsi,  Event::mKaon,  Event::mPion, LR, LB,  radius );
    }
    
    res->setName( name );
    std::cout << "Debug: for reson: H0,Hp,Hm " << name << " " << H0 << " " <<Hp << " " <<Hm << std::endl;
    std::cout << "Debug: mass, width" << mass << " " << width <<  std::endl;
    
    res->setAmplitudes( H0, Hp, Hm );
    resonances.push_back( res );
  }
  
  if (resonance == "FLATTE"){
    
    double mass  = pt.get<double>("mean");
    double gpipi = pt.get<double>("gpipi");
    double gkk   = pt.get<double>("gkk");
    double alpha = pt.get<double>("alpha");
    
    boost::optional<double> opt_radius =
      pt.get_optional<double>("radius");
    
    double radius = ( opt_radius ) ? opt_radius.get() : 1.6 ;
    
    res = new KPiResonanceFlatte(mass,gkk,gpipi,alpha,Event::mB,Event::mPsi,Event::mKaon,Event::mPion, 1, radius);
    res->setName(name);
    res->setAmplitudes( H0, Hp, Hm );
    resonances.push_back( res );
  }
  
  if ( resonance == "LASS" ){ 
    
    res = new KPiResonanceLASS(  Event::mB,  Event::mPsi,  Event::mKaon,  Event::mPion );
    res->setName( "LASS" );
    res->setAmplitudesPolar( amp[0], amp[1] );
    
    resonances.push_back( res );
  }
  
  return ;
}


void parse_exotic( const std::string& name,
		   std::vector< PsiPiResonance* >& resonances,
		   const ptree& pt ){
    
  PsiPiResonance* res;
  
  double mass   = pt.get<double>("mean");
  double width  = pt.get<double>("width");
  
  boost::optional<double> opt_radius =
    pt.get_optional<double>("radius");
  
  double radius = ( opt_radius ) ? opt_radius.get() : 1.6 ;
  
  int J = pt.get<int>("J");
  std::string parity = pt.get<std::string>("P");
  
  std::vector<double> amp = list_from_tree<double>( pt, "amplitude" );
  
  res = new PsiPiResonance( mass, width, Event::mB, Event::mPsi, Event::mKaon, Event::mPion, J, 0, radius, radius );
  res->setName( name );
  
  //
  if ( J > 0 ){
    res->setAmplitudesPolar( amp[0], amp[1], amp[2], amp[3], amp[4], amp[5] );
  }
  else {
    res->setAmplitudesPolar( amp[0], amp[1], 0, 0, 0, 0 );
  }
  
  resonances.push_back( res );
  
  return ;
}


void parse_masses( const ptree& pt ){ 

  boost::optional< double > massB = pt.get_optional<double>("mMother");
  if ( massB ) Event::mB = massB.get();
  
  boost::optional< double > massPsi = pt.get_optional<double>("mPsi");
  if ( massPsi ) Event::mPsi = massPsi.get();
  
  boost::optional< double > massh1 = pt.get_optional<double>("mh1");
  if ( massh1 ) Event::mKaon = massh1.get();
  
  boost::optional< double > massh2 = pt.get_optional<double>("mh2");
  if ( massh2 ) Event::mPion = massh2.get();
  
  return ; 
}


void parse_resonance_list( const std::string json,
                           std::vector<KPiResonanceBase*>& kpi_resonances,
			   std::vector<PsiPiResonance*>& psipi_resonances,
			   const bool debug ){

    //Used to store the Azero for the first resonance to scale other resonances against (fit fracs)
    ptree pt;
   
    read_json( json ,  pt );
    
    if ( debug ) {
      print_tree( pt, 0 );
    }
    
    //Mode 
    std::string mode{""};
    boost::optional<std::string> ismode = pt.get_optional<std::string>("mode");
    mode =  ismode ? ismode.get() : "cartesian";

    //D factor for B/Bs and bar interfering
    double D{0};
    boost::optional<double> isD = pt.get_optional<double>("Dfactor");
    D = isD ? isD.get() : 0.01;
    Event::D = D;
    std::cout << "DEBUG: D factor is" << Event::D << std::endl;

    // This comes first to ensure correct daughter masses are picked up
    Event::mB    = Variables::mB;
    Event::mPsi  = Variables::mPsi ;
    Event::mKaon = Variables::mKaon;
    Event::mPion = Variables::mPion;
    
    if ( pt.get_child_optional( "masses" ) ){
      parse_masses( pt.get_child( "masses" ) );
    }

    ///Run the frac version of parse_resonance and parse_exotic
    for ( const auto& i : pt.get_child("resonance") ){
	    parse_resonance( i.first, kpi_resonances, i.second,mode);
    }

    if ( pt.get_child_optional( "exotic" ) ){
	    for ( const auto& i : pt.get_child("exotic") ){
		    parse_exotic( i.first, psipi_resonances, i.second);
	    }
    }

    return;
}

void parse_amplitudes( const unsigned int LR, 
		       const ptree& pt, 
		       std::complex<double>& H0, 
		       std::complex<double>& Hp, 
		       std::complex<double>& Hm ){ 
  
  std::vector<double> amp = list_from_tree<double>( pt, "amplitude" );
  
  const std::string basis = pt.get<std::string>("basis");
  const std::string mode  = pt.get<std::string>("mode" );
  
  if ( LR == 0 ) { 
    if ( mode == "polar" ) { 
      H0 = std::polar(amp[0],amp[1]);  
    } else {
      H0 = std::complex<double>( amp[0], amp[1] );
    }
  }
  else { 
    if ( basis == "transversity" ){ 
      if ( mode == "polar" ){ 
	std::complex<double> Azero = std::polar(amp[0],amp[3]);
	std::complex<double> Apara = std::polar(amp[1],amp[4]);
	std::complex<double> Aperp = std::polar(amp[2],amp[5]);

	Helper::convertToHelicityBasis(Azero,Apara,Aperp,H0,Hp,Hm);
      } else { 
	std::complex<double> Azero(amp[0],amp[3]);
	std::complex<double> Apara(amp[1],amp[4]);
	std::complex<double> Aperp(amp[2],amp[5]);

	Helper::convertToHelicityBasis(Azero,Apara,Aperp,H0,Hp,Hm);
      }
      
    }
    else {
      if ( mode == "polar" ){ 
	H0 = std::polar(amp[0],amp[3]);
	Hp = std::polar(amp[1],amp[4]);
	Hm = std::polar(amp[2],amp[5]);
      } else { 
	H0 = std::complex<double>(amp[0],amp[3]);
	Hp = std::complex<double>(amp[1],amp[4]);
	Hm = std::complex<double>(amp[2],amp[5]);
      }
    }
  }
  
  return ; 
}
