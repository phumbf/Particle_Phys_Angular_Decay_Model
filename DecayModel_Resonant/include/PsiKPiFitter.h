#ifndef PSIKPIFITTER_H
#define PSIKPIFITTER_H

#include "Event.h" 
#include "PsiKPiAmplitudeModel.h" 

// std
#include <vector> 
#include <string> 

// ROOT 
#include "TH1.h" 
#include "TTree.h" 

// Minuit2 
#include "Minuit2/MnUserParameters.h" 

class PsiKPiFitter { 
 public:
  
  PsiKPiFitter() ; 

  PsiKPiFitter( const std::string& json) ; 
  
  PsiKPiFitter( const PsiKPiFitter& other );

  ~PsiKPiFitter() ;

  double operator()( const std::vector< double >& par ) ;
  
 public:
  
  TH1* legendre( const char* name, const unsigned int order, const std::vector< Event >& collection ) const ;
  
 public: 
  
  void addMC  ( const Long64_t entries );
    
  void addMC  ( TTree* tree,  const Long64_t entries );

  void addData( TTree* tree,  const Long64_t entries );

  unsigned getParameters( ROOT::Minuit2::MnUserParameters& par ) ; 
  
 private:
  
  void addToCollection( std::vector< Event >& collection, TTree* tree, const Long64_t entries );
  
  double integrate( const std::vector< Event >& collection ) const ;

 private:

  // Model 
  PsiKPiAmplitudeModel model_; 
  
  // Data and MC stores
  std::vector< Event > data_;
  std::vector< Event > mc_  ;
}; 

#endif 
