#include "OptionsParser.h" 

// ROOT 
#include "Minuit2/MinuitParameter.h" 


OptionsParser::OptionsParser() : 
  nsample   ( 1000000 ), 
  appname   ( "unknown" ), 
  filename  ( "generated.root" ), 
  modelname ( "BdKstarMuMu/Amplitude/config/model.json" ){};

OptionsParser::OptionsParser( const OptionsParser& other ) : 
  nsample   ( other.nsample ), 
  appname   ( other.appname ), 
  filename  ( other.filename ), 
  modelname ( other.modelname ){};

OptionsParser::OptionsParser( int argc, char* argv[] ) : 
  nsample   ( 1000000 ) , 
  appname   ( "" ),
  filename  ( "generated.root" ),
  modelname ( "BdKstarMuMu/Amplitude/config/model.json" )
{
  appname = std::string( argv[0] );
  
  for ( int i=1 ; i < argc; ++i ){
    std::string arg = argv[i];
    
    if ( 0 == arg.compare("-help") ){ 
      printUsageAndExit();
    }
    
    if ( 0 == arg.compare("-output") ){ 
      filename = std::string( argv[++i] ); 
    }
    if ( 0 == arg.compare("-model") ) {
      modelname = std::string( argv[++i] ); 
    }
    if ( 0 == arg.compare("-samples") ){ 
      nsample = std::stoul( std::string( argv[++i] ) );
    }
    if ( 0 == arg.compare("-resonance")){
	resonance = std::string(argv[++i]);
    }
  }
};


void OptionsParser::printUsageAndExit(){ 
  std::cout 
    << " ==> Usage of " << appname << "\n\n"
    << "     -help    -- Prints this message \n" 
    << "     -output  -- Specify the name of the output ROOT file (default: generated.root) \n" 
    << "     -model   -- Specify the model to use                 (default: config/model.json) \n" 
    << "     -samples -- Specify the number of points to sample   (default: one million) \n" 
    << std::endl;
  
  exit(0);
}

std::ostream& operator<<( std::ostream& stream, const OptionsParser& opt ){
  
  stream << " ==> External configuration: " << "\n" ;
  stream << "     Application           = " << opt.appname   << "\n";
  stream << "     Output file           = " << opt.filename  << "\n";
  stream << "     Model                 = " << opt.modelname << "\n";
  stream << "     Sample size           = " << opt.nsample   << "\n";

  return stream;
}

void print_configuration( const ROOT::Minuit2::MnUserParameters& params ){
  const std::vector< ROOT::Minuit2::MinuitParameter >& par =
    params.Parameters();

  std::cout << " ==> INFO: Minuit2 input parameters \n" << std::endl;

  for ( unsigned int i = 0; i < par.size(); ++i ){
    std::cout 
      <<  "       (" << i << ") " 
      << par[i].Name()  
      << " : value = "   
      << par[i].Value() << "\n" ;
  }

  std::cout << std::endl;

  return ;
}
