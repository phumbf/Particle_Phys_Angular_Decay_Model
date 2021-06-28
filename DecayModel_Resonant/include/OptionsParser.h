#ifndef OPTIONSPARSER_H
#define OPTIONSPARSER_H

#include <string> 
#include <iostream> 

// ROOT 
#include "Minuit2/MnUserParameters.h" 

//Class to handle the parsing in of options from the different .json files
class OptionsParser { 
  
 public: 
  OptionsParser() ; 
  
  OptionsParser( int argc, char* argv[] ) ;

  OptionsParser( const OptionsParser& other ); 
  
  ~OptionsParser(){} ;
  
  friend std::ostream& operator<<( std::ostream& stream, const OptionsParser& opt );

 public:
  unsigned int nsample ; 
  
  std::string appname  ; 
  
  std::string filename ;

  std::string modelname; 

  std::string resonance;

  bool normalise;

 protected:
  void printUsageAndExit(); 
};

void parse_configuration( const ROOT::Minuit2::MnUserParameters& params ) ; 

#endif 
