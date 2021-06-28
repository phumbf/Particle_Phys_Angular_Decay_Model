// Boost
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

// std
#include <iostream>

//Root
#include "TMatrixDSym.h"

//Local
#include "ResonanceBase.h"
#include "ScalarResonance.h"
#include "VectorResonance.h"
#include "ResonanceParser.h"
#include "FFparser.h"
#include "Parameters.h" 
#include "Wilson.h"

//Class to allow the resonance parser to be read in 
using boost::property_tree::ptree;

Resparser::Resparser() {}

Resparser::~Resparser() {}

std::vector< std::unique_ptr<ResonanceBase> > 
Resparser::getResonances( const std::string filename ) const{
  
  if ( m_debug ){ 
    std::cout 
      << "Parsing file " << filename 
      << std::endl;
  }
  
  ptree pt;
  read_json( filename, pt );

  std::vector< std::unique_ptr<ResonanceBase> > list;

  const double mothermass = pt.get<double>("MotherM");

  //Wilson coeffs
  const std::string WCstr = pt.get<std::string>("WC");
  const std::string WPstr = pt.get<std::string>("WP");
  
  //Daughter masses
  const double daugM1 = pt.get<double>("Daughter1M");
  const double daugM2 = pt.get<double>("Daughter2M");

  Wilson WC, WP;
  if(WCstr == "SM"){
    WC.init(Wilson::SM);
  }
  else if(WCstr == "varyC9"){
    WC.init(Wilson::varyC9);
  }
  else if(WCstr == "varyC9C10"){
    WC.init(Wilson::varyC9C10);
  }
  else{
    WC.init(Wilson::Zero);
  }

  WPstr == "SM" ? WP.init(Wilson::SM) : WP.init(Wilson::Zero);

  for ( const auto& i : pt.get_child("resonance") ){
    
    std::unique_ptr< ResonanceBase > res;
    
    const std::string name = i.first; 
    const double mass =  i.second.get<double>("mesonM") ;
    const int L       =  i.second.get<int>("mesonL") ;
    
    std::string lineshape      = i.second.get<std::string>("lineshape") ;
    std::string formfactorjson = i.second.get<std::string>("formfactors") ;
    
    const double BFfac =  i.second.get<double>("BFfac") ;
    const bool isCP = i.second.get<bool>("isCP");
    const bool isEvtGenWC = i.second.get<bool>("isEvtGenWC");
    const bool isbtod = i.second.get<bool>("isbtod");
    
    // Resonance masses for form-factors
    boost::optional<double> B0m = i.second.get_optional<double>("B0m");
    boost::optional<double> B1p = i.second.get_optional<double>("B1p");
    boost::optional<double> B1m = i.second.get_optional<double>("B1m");
    boost::optional<double> Bst = i.second.get_optional<double>("Bst");
    
    boost::optional<double> cR = i.second.get_optional<double>("couplingReal");
    boost::optional<double> cI = i.second.get_optional<double>("couplingImag");
    boost::optional<double> cV = i.second.get_optional<double>("couplingval");
    boost::optional<double> cA = i.second.get_optional<double>("couplingpha");

    double mB0m = B0m ? B0m.get() : 0 ;
    double mB1m = B1m ? B1m.get() : 0 ;
    double mB1p = B1p ? B1p.get() : 0 ;
    double mBst = Bst ? Bst.get() : 0 ;
    
    if( L == 0 ){
      auto r = std::make_unique<ScalarResonanceBase>();
      
      std::unique_ptr<ScalarFormFactorBase> ff = 
	getScalarFFs(formfactorjson); 
      
      ff->setMass( mass );
      ff->setMotherMass( mothermass );
      ff->setB0m( mB0m );
      ff->setB1m( mB1m );
      ff->setB1p( mB1p );
      ff->setBst( mBst );
      
      r->setFFs( ff );
      r->setBFfac(BFfac);
      r->setisCP(isCP);
      r->setisEvtGenWC(isEvtGenWC);
      r->setisbtod(isbtod);
      if(cR && cI){
        r->setResFactor(cR.get(),cI.get());
      }
      else if(cV && cA){
        r->setPolarResFactor(cV.get(),cA.get());
      }
      
      res = std::move( r );
    }
    else{
      auto r = std::make_unique<VectorResonanceBase>();
      
      std::unique_ptr<VectorFormFactorBase> ff = 
	getVectorFFs(formfactorjson);

      ff->setMass( mass );
      ff->setMotherMass( mothermass );
      ff->setB0m( mB0m );
      ff->setB1m( mB1m );
      ff->setB1p( mB1p );
      ff->setBst( mBst );
      
      r->setFFs( ff );
      r->setBFfac(BFfac);
      r->setisCP(isCP);
      r->setisEvtGenWC(isEvtGenWC);
      r->setisbtod(isbtod);
      if(cR && cI){
          r->setResFactor(cR.get(),cI.get());
      }
      else if(cV && cA){
          r->setPolarResFactor(cV.get(),cA.get());
      }  
      res = std::move( r );
    }
    
    // set resonance properties
    res->setWilsons(WC,WP);
    res->setName( name );
    res->setMotherMass( mothermass );
    res->setMass( mass );
    res->setL( L );
    res->setDaughterMasses( daugM1, daugM2 );
    
    if(lineshape == "BW"){
	std::cout 
	  << "Building BW lineshape" 
	  << std::endl;
      
      const double width = i.second.get<double>("mesonW") ;
	    
      std::unique_ptr<HHLineShapeBase> shape =  
	std::make_unique<BWShape>( mass, width, mothermass, daugM1, daugM2 );

      res->setShape( shape );
    }
    else if(lineshape == "CoupledBW"){
      if ( m_debug ){ 
	std::cout 
	  << "Building Coupled BW lineshape" 
	  << std::endl;
      }
	    
      const double mK = i.second.get<double>("mK");
      const double gK = i.second.get<double>("gK");
      const double mR = i.second.get<double>("mR");
      const double gR = i.second.get<double>("gR");
      const double cR = i.second.get<double>("couplingReal");
      const double cI = i.second.get<double>("couplingImag");

      std::complex<double> coupling = std::polar( cR, cI ); 
      
      std::unique_ptr<HHLineShapeBase> shape = 
	std::make_unique<BWCoupledShape>(mK,gK,mR,gR,coupling,mothermass,daugM1,daugM2);

      res->setShape( shape );
    }
    else if(lineshape == "3CoupledBW"){
      if ( m_debug ){ 
	std::cout 
	  << "Building 3Coupled BW lineshape" 
	  << std::endl;
      }
	    
      const double m1 = i.second.get<double>("m1");
      const double g1 = i.second.get<double>("g1");
      const double m2 = i.second.get<double>("m2");
      const double g2 = i.second.get<double>("g2");
      const double m3 = i.second.get<double>("m3");
      const double g3 = i.second.get<double>("g3");
      const double c12real = i.second.get<double>("c12real");
      const double c12imag = i.second.get<double>("c12imag");
      const double c13real = i.second.get<double>("c13real");
      const double c13imag = i.second.get<double>("c13imag");

      std::complex<double> coupling12 = std::polar( c12real, c12imag ); 
      std::complex<double> coupling13 = std::polar( c13real, c13imag ); 
      
      std::unique_ptr<HHLineShapeBase> shape = 
	std::make_unique<BW3CoupledShape>(m1,g1,m2,g2,m3,g3,coupling12,coupling13,mothermass,daugM1,daugM2);

      res->setShape( shape );
    }

    else{
      std::cout << "ERROR: Lineshape is not recognised" << std::endl;
    }
    
    list.push_back( std::move(res) );
  }
  
  return list;
}

std::unique_ptr<VectorFormFactorBase> Resparser::getVectorFFs( const std::string filename ) const{
  
  FFparser parser;
  parser.setDebug( m_debug );

  std::unique_ptr<VectorFormFactorBase> ff;
  
  if ( filename.find("BZ04") != std::string::npos ){
    
    std::unique_ptr<BZ04FFCalculation> bz( parser.getBZ04VectorFormFactors(filename) );
    ff = std::move(bz);
  }
  else if ( filename.find("BSZ15") != std::string::npos ) {
    
    std::unique_ptr<BSZ15FFCalculation> bsz(parser.getBSZ15VectorFormFactors( filename ));
    ff = std::move(bsz);
  }
  else if (filename.find("YANG") != std::string::npos ) {
    std::unique_ptr<YANG10FFCalculation> yang(parser.getYANG10VectorFormFactors( filename ));
    
    ff = std::move(yang);
  }
  else if (filename.find("ESS18") != std::string::npos ) {
    std::unique_ptr<ESS18FFCalculation> yang(parser.getESS18VectorFormFactors( filename ));
    
    ff = std::move(yang);
  }
  else if (filename.find("MS00") != std::string::npos ){
    std::unique_ptr<MS00FFCalculation> abhh(parser.getMS00VectorFormFactors( filename ));
    
    ff = std::move(abhh);
  }

  else if (filename.find("ABHH99") != std::string::npos ){
    std::unique_ptr<ABHH99FFCalculation> abhh(parser.getABHH99VectorFormFactors( filename ));
    
    ff = std::move(abhh);
  }
  return ff;
}

std::unique_ptr<ScalarFormFactorBase> Resparser::getScalarFFs( const std::string filename ) const{
  
  FFparser parser;
  parser.setDebug( m_debug );
  
  
  std::unique_ptr<ScalarFormFactorBase> ff;
  
  if(filename.find( "AAS07") != std::string::npos ){
    
    std::unique_ptr<AAS07FFCalculation> aas( parser.getAAS07ScalarFormFactors( filename ) );
    ff = std::move(aas);
  }
  
  else{
    std::unique_ptr<CFW10FFCalculation> cfw( parser.getCFW10ScalarFormFactors( filename ) );
    ff = std::move(cfw);
  }
  
  return ff;
}

