
// Local 
#include "Variables.h" 
#include "PsiKPiFitter.h" 

// std
#include <iostream> 

// ROOT 
#include "TH1D.h" 
#include "TMath.h"

// ROOT Math
#include "Math/QuasiRandom.h"
#include "Math/SpecFunc.h"


PsiKPiFitter::PsiKPiFitter(){} ; 

PsiKPiFitter::PsiKPiFitter( const std::string& json ) : 
  model_( json ) {} ; 


PsiKPiFitter::PsiKPiFitter( const PsiKPiFitter& other ) : 
  model_( other.model_ ), 
  data_ ( other.data_ ), 
  mc_   ( other.mc_ ){} ;


PsiKPiFitter::~PsiKPiFitter(){ 
  data_.clear(); 
  mc_.clear();
}


double PsiKPiFitter::operator()( const std::vector<double>& par ) {
  /* 
   * operator for use in fitting 
   */
  
  double result = 0;
  double prob   = 0;
  
  // Set the parameters of the model 
  model_.setParameters( par ) ; 
  
  // compute the integral over phasespace using MC integration
  double integral = integrate( mc_ );
  
  // sum over the events and evaluate the likelihood
  for ( auto evt: data_ ){ 
    
    prob = model_.evaluate( evt );
    
    if ( prob > 0 ){
      result -= 2.*log( prob/integral );
    }
    else {
      result -= 100;
    }
  }
  
  return result;
}

double PsiKPiFitter::integrate( const std::vector< Event >& collection ) const {
  /* 
   * Integrate the model PDF by MC integration 
   */
  
  const double V = 8*TMath::Pi()*( Variables::mB    - Variables::mPsi -
				   Variables::mKaon - Variables::mPion );
  
  double result = 0;
  
  for ( auto evt: collection ) { 
    result += model_.evaluate( evt );
  }
  
  return V*result/double(collection.size());
}


void PsiKPiFitter::addToCollection( std::vector< Event >& collection, TTree* tree, const Long64_t entries ) { 
  /* 
   * Add events from TTree to either the MC or data store 
   */
  
  Double_t kpi, ctl, ctk, phi;
  Int_t    id;
  
  if ( tree ) {
    Long64_t maxentries = 
      ( entries > tree->GetEntries() ? tree->GetEntries() : entries );
    
    if ( maxentries < entries ){
      std::cout 
	<< " ==> WARNING (PsiKPiFitter) tree size is smaller than requested events"
	<< std::endl;
    }
    
    tree->SetBranchAddress( "mKpi" , &kpi );
    tree->SetBranchAddress( "costhetal", &ctl );
    tree->SetBranchAddress( "costhetak", &ctk );
    tree->SetBranchAddress( "phi", &phi );
    tree->SetBranchAddress( "ID", &id );
    
    for ( Long64_t i = 0; i < maxentries; ++i ){
      tree->GetEntry( i );
      collection.push_back( Event( id, kpi, ctk, ctl, phi  ) );
    }
  }
  else {
    std::cout 
      << " ==> WARNING (PsiKPiFitter) NULL pointer for TTree"
      << std::endl;
  }
  
  return ;
}

void PsiKPiFitter::addData( TTree* tree , const Long64_t entries ) {
   /* 
   * Add events to the data store 
   */
  
  data_.clear();
    
  addToCollection( data_, tree, entries );
  
  std::cout 
    << " ==> INFO (PsiKPiFitter) added "
    << data_.size() 
    << " events to data container"
    << std::endl;
  
  return;
}

void PsiKPiFitter::addMC( TTree* tree , const Long64_t entries ) {
  /* 
   * Add events to the MC store 
   */
  
  mc_.clear();
  
  addToCollection( mc_, tree, entries );
  
  std::cout 
    << " ==> INFO (PsiKPiFitter) added "
    << data_.size() 
    << " events to MC container"
    << std::endl;
  
  return;
}


void PsiKPiFitter::addMC( const Long64_t entries ) {
  /* 
   * Generate pseudo-random MC points and 
   * add the resulting events to the MC store 
   */
  
  mc_.clear();
  
  double kpi, ctl, ctk, phi;
  
  ROOT::Math::QuasiRandomSobol rndm(4);
  double x[4];
  
  const double min_kpi = Variables::mKaon + Variables::mPion;
  const double max_kpi = Variables::mB - Variables::mPsi;
  
  for ( Long64_t i = 0; i < entries ; ++i ){
    
    rndm.Next(x);
    
    kpi =  (max_kpi - min_kpi)*x[0] + min_kpi;
    ctk =  2.0*x[1] - 1.0 ;
    ctl =  2.0*x[2] - 1.0 ;
    phi =  TMath::TwoPi()*x[3] - TMath::Pi();
    
    mc_.push_back( Event( 511, kpi, ctk, ctl, phi  ) );
  }
  
  return ;
}

TH1* PsiKPiFitter::legendre( const char* name, const unsigned int order, const std::vector< Event >& collection ) const {
  
  /*
   *  Evaluate Legendre polynomial in cos theta_kpi 
   */
  
  const unsigned int nbins = 20;
  
  TH1D* hist = new TH1D( name, name, nbins,
			 Variables::mKaon + Variables::mPion,
			 Variables::mB - Variables::mPsi );
  hist->Sumw2();
  
  for ( auto evt: collection ){
    hist->Fill( evt.mkpi_, ROOT::Math::legendre( order, evt.coskpi_ ) );
  }
  
  return hist;
}

unsigned PsiKPiFitter::getParameters( ROOT::Minuit2::MnUserParameters& par ) { 
  
  unsigned int nparam = model_.getParameters( par ); 
  
  return nparam;  
}
