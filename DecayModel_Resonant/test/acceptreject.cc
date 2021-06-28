#include "Event.h"
#include "PsiKPiAmplitudeModel.h"
#include "Variables.h" 
#include "OptionsParser.h" 


#include <iostream>
#include <string> 

#include "TGraph.h"
#include "TFile.h"
#include "TStopwatch.h"

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"

//A simple accept-reject test to see whether the model produces the correct distributions from the literature

double get_maximum( PsiKPiAmplitudeModel& model, const unsigned int nburn ){ 
  
  double result = 0; 
  
  double mkpi, cospsi, coskpi, phi ;

  for ( unsigned int i =0 ; i < nburn; ++i ){
    
    mkpi   = gRandom->Uniform( Event::mKaon + Event::mPion, Event::mB - Event::mPsi );
    coskpi = gRandom->Uniform( -1, 1 );
    cospsi = gRandom->Uniform( -1, 1 );
    phi    = gRandom->Uniform( -TMath::Pi(), TMath::Pi() );
    
    Event evt( 511, mkpi, coskpi, cospsi, phi );
    
    double value = model.evaluate( evt );
    
    if ( value > result ) { 
      result = value; 
    }
  }
  
  return 1.2*result;
}


int main( int argc, char* argv[] ) {
  
  OptionsParser parser( argc, argv );
  std::cout << parser << std::endl;
  
  PsiKPiAmplitudeModel model;
  model.initialise( parser.modelname );
	

  model.print();

  TFile* file = new TFile( parser.filename.c_str() ,"RECREATE");
  file->cd();
  
  TTree* tree = new TTree("Events","Events");
  
  double wgt; 
  double mkpi, mz;
  double cospsi;
  double coskpi;
  double phi;
  
  tree->Branch("mkpi",&mkpi,"mkpi/D");
  tree->Branch("coskpi",&coskpi,"coskpi/D");
  tree->Branch("cospsi",&cospsi,"cospsi/D");
  tree->Branch("phi",&phi,"phi/D");
  tree->Branch("mz", &mz , "mz/D");
    
  TStopwatch timer;
  timer.Start();
  
  
  double max_wgt = get_maximum( model, 100000 );
  
  for ( unsigned int i = 0; i < parser.nsample ; ++i ){
    
    do { 
      mkpi   = gRandom->Uniform( Event::mKaon + Event::mPion, 
				 Event::mB - Event::mPsi );
      coskpi = gRandom->Uniform( -1, 1 );
      cospsi = gRandom->Uniform( -1, 1 );
      phi    = gRandom->Uniform( -TMath::Pi(), TMath::Pi() );
      
      Event evt(511, mkpi,coskpi,cospsi,phi);
      
      mz     = evt.mz_;
      wgt    = model.evaluate( evt ) ;
   

    } while ( gRandom->Rndm() > wgt/max_wgt ) ;
    
    tree->Fill();
  }
  
  timer.Stop();
  
  std::cout << " time elapsed = " << timer.CpuTime() << std::endl;
  
  
  tree->Write();
  file->Close();
  
  return 0;
}
