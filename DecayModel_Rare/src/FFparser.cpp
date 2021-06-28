// Boost
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

// std
#include <iostream>

// Local
#include "FFparser.h"
#include "FFcalculation.h"

//Root
#include "TMatrixDSym.h"

//Class to parse in the form factor information from json files
using boost::property_tree::ptree;

FFparser::FFparser() {} 

FFparser::~FFparser(){}


std::unique_ptr<BZ04FFCalculation>
FFparser::getBZ04VectorFormFactors( const std::string filename) const { 

// define ptree and read json file
ptree pt;
read_json( filename, pt );
  
std::vector<double> params;

for ( const auto& i : pt.get_child("central") ){
    
    const std::string ffname = i.first;
    
    if ( m_debug ){
      std::cout << "Reading parameters for " << ffname << std::endl;
    }
    
    boost::optional<double> r1  = i.second.get_optional<double>("r1");
    boost::optional<double> r2  = i.second.get_optional<double>("r2");
    boost::optional<double> r   = i.second.get_optional<double>("r");
    boost::optional<double> msq = i.second.get_optional<double>("msq");
    
    if( r1 ) { params.push_back( r1.get() );  }
    if( r2 ) { params.push_back( r2.get() );  }
    if( r  ) { params.push_back( r.get()  );  }
    if( msq ){ params.push_back( msq.get() ); }
  }
  
  std::unique_ptr<BZ04FFCalculation> ff = 
    std::make_unique<BZ04FFCalculation>( params );
  
  //ff->setParams( params );
  
  return ff;
}

std::unique_ptr<BSZ15FFCalculation>
FFparser::getBSZ15VectorFormFactors( const std::string filename ) const { 

  // define ptree and read json file
  ptree pt;
  read_json( filename, pt );

  std::vector<double> params;

  for ( const auto& i : pt.get_child("central") ){
    
    const std::string ffname = i.first;
    
    if ( m_debug ){ 
      std::cout << "Reading parameters for " << ffname << std::endl;
    }
    
    const double a0 = i.second.get<double>("a0");
    const double a1 = i.second.get<double>("a1");
    const double a2 = i.second.get<double>("a2");
    
    params.push_back( a0 );
    params.push_back( a1 );
    params.push_back( a2 );
  }
  
  std::unique_ptr<BSZ15FFCalculation> ff = 
    std::make_unique<BSZ15FFCalculation>( params );

  return ff;
}


std::unique_ptr<YANG10FFCalculation> 
FFparser::getYANG10VectorFormFactors( const std::string filename ) const { 
  
  std::vector< double > params = getFabFormFactors( filename );

  std::unique_ptr<YANG10FFCalculation> ff = 
    std::make_unique<YANG10FFCalculation>( params ); 
  //ff->setParams( params );
  
  return ff; 
}

std::unique_ptr<ESS18FFCalculation> 
FFparser::getESS18VectorFormFactors( const std::string filename ) const { 
  
  std::vector< double > params = getFabFormFactors( filename );

  std::unique_ptr<ESS18FFCalculation> ff = 
    std::make_unique<ESS18FFCalculation>( params ); 
  //ff->setParams( params );
  
  return ff; 
}

std::unique_ptr<MS00FFCalculation> 
FFparser::getMS00VectorFormFactors( const std::string filename ) const { 
  
  std::vector< double > params = getFabFormFactors( filename );

  std::unique_ptr<MS00FFCalculation> ff = 
    std::make_unique<MS00FFCalculation>( params ); 
  //ff->setParams( params );
  
  return ff; 
}

std::unique_ptr<ABHH99FFCalculation> 
FFparser::getABHH99VectorFormFactors( const std::string filename ) const { 
  
  std::vector< double > params = getCCCFormFactors( filename );

  std::unique_ptr<ABHH99FFCalculation> ff = 
    std::make_unique<ABHH99FFCalculation>( params ); 
  //ff->setParams( params );
  
  return ff; 
}

std::unique_ptr<AAS07FFCalculation>  
FFparser::getAAS07ScalarFormFactors ( const std::string filename ) const {
  
  std::vector< double > params = getFabFormFactors( filename );
  
  std::unique_ptr<AAS07FFCalculation> ff = 
    std::make_unique<AAS07FFCalculation>( params );
  
  return ff;
}

std::unique_ptr<CFW10FFCalculation>  
FFparser::getCFW10ScalarFormFactors ( const std::string filename ) const {
  
  std::vector< double > params = getFabFormFactors( filename );

  std::unique_ptr<CFW10FFCalculation> ff = 
    std::make_unique<CFW10FFCalculation>( params ); 
  
  return ff;
}


std::vector<double> FFparser::getFabFormFactors( const std::string filename ) const { 

  ptree pt;
  read_json( filename , pt );
  
  std::vector<double> params;

  for ( const auto& i : pt.get_child("central") ){
    
    const std::string ffname = i.first;
    
    if ( m_debug ){ 
      std::cout << "Reading parameters for " << ffname << std::endl;
    }
    
    const double f = i.second.get<double>("f0");
    const double a = i.second.get<double>("a");
    const double b = i.second.get<double>("b");
  
    params.push_back( f );
    params.push_back( a );
    params.push_back( b );
  }

  return params;
}

std::vector<double> FFparser::getCCCFormFactors( const std::string filename ) const { 

  ptree pt;
  read_json( filename , pt );
  
  std::vector<double> params;

  for ( const auto& i : pt.get_child("central") ){
    
    const std::string ffname = i.first;
    
    if ( m_debug ){ 
      std::cout << "Reading parameters for " << ffname << std::endl;
    }
    
    const double f = i.second.get<double>("f0");
    const double c1 = i.second.get<double>("c1");
    const double c2 = i.second.get<double>("c2");
    const double c3 = i.second.get<double>("c3");
  
    params.push_back( f );
    params.push_back( c1 );
    params.push_back( c2 );
    params.push_back( c3 );
  }

  return params;
}

TMatrixDSym FFparser::dropParams( const TMatrixDSym& M, const std::set< int >& droppedParams ) const { 
  
  for ( auto& p: droppedParams ){ 
    assert( p >= 0 && p < M.GetNrows() );
  }

  const int noriginal = M.GetNrows();
  const int nfinal    = noriginal - droppedParams.size();

  TMatrixDSym result( nfinal );
  
  // rest counter
  int x = 0;
  int y = 0;
  
  for ( int i = 0; i < noriginal; i++ ){
    
    // check against drop list
    if ( droppedParams.count( i ) ){ continue;  }
    
    for ( int j = 0; j < noriginal; j++ ){
      
      if ( droppedParams.count( j ) ){ continue; }
      
      result( x, y ) = M( i, j );
      
      y++;
    }
    
    y = 0;
    x++;
  }
  
  return result;
}

TVectorD FFparser::dropParams( const TVectorD& V, const std::set< int >& droppedParams ) const { 
  for ( auto& p: droppedParams ){ 
    assert( p >= 0 && p < V.GetNrows() );
  }

  const int noriginal = V.GetNrows();
  const int nfinal    = noriginal - droppedParams.size();

  TVectorD result( nfinal );
  
  int n = 0;

  for ( int i = 0; i < noriginal; i++ ){
    if ( droppedParams.count( i ) ){ continue; } 
    
    result( n ) = V( i ); n++;
  }
  
  return result;
}

TVectorD FFparser::getfabUncertainty( const std::string filename ) const {
  
	ptree pt;
	read_json( filename , pt );
	
	std::vector<double> vec;

	for ( const auto& i : pt.get_child("uncertainty") ){

		vec.push_back( i.second.get<double>("f0"));
		vec.push_back( i.second.get<double>("a"));
		vec.push_back( i.second.get<double>("b"));
	}

	TVectorD retvec(vec.size(),&vec[0]);

	return retvec;
}

TVectorD FFparser::getaaaUncertainty( const std::string filename ) const {
  
	ptree pt;
	read_json( filename , pt );
	
	std::vector<double> vec;

	for ( const auto& i : pt.get_child("uncertainty") ){

		vec.push_back( i.second.get<double>("a0"));
		vec.push_back( i.second.get<double>("a1"));
		vec.push_back( i.second.get<double>("a2"));
	}

	TVectorD retvec(vec.size(),&vec[0]);

	return retvec;
}

TMatrixDSym FFparser::getBSZ15Covariance( const std::string filename ) const {
	
  ptree pt;
  read_json( filename , pt );

  const int npar = 21; 
  int x{0};
  int y{0};

  TMatrixDSym M( npar );
  
  // Read raw covariance matrix
  for ( const auto& i : pt.get_child("covariance") ){
    
    const std::string ffname = i.first;
    
    std::cout 
      << "Reading covariance for combination " 
      << ffname 
      << std::endl;
    
    const double a0a0 = i.second.get<double>("a0a0");
    const double a0a1 = i.second.get<double>("a0a1");
    const double a0a2 = i.second.get<double>("a0a2");
    const double a1a0 = i.second.get<double>("a1a0");
    const double a1a1 = i.second.get<double>("a1a1");
    const double a1a2 = i.second.get<double>("a1a2");
    const double a2a0 = i.second.get<double>("a2a0");
    const double a2a1 = i.second.get<double>("a2a1");
    const double a2a2 = i.second.get<double>("a2a2");
    
    if ( y >= npar ){
      y  = 0;
      x += 3;
    }

    M( x  , y   ) = a0a0;
    M( x  , y+1 ) = a0a1;
    M( x  , y+2 ) = a0a2;
    M( x+1, y   ) = a1a0;
    M( x+1, y+1 ) = a1a1;
    M( x+1, y+2 ) = a1a2;
    M( x+2, y   ) = a2a0;
    M( x+2, y+1 ) = a2a1;
    M( x+2, y+2 ) = a2a2;
    
    y += 3; 
  }
  
  return M;
}
