#include "Event.h"
#include "PsiKPiAmplitudeModel.h"
#include "Event.h" 
#include "OptionsParser.h" 
#include "Parser.h"

#include <iostream>
#include <string> 

#include "TGraph.h"
#include "TFile.h"
#include "TStopwatch.h"

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TH1D.h"

//Test to determine the mass distribution for a Bs mother
double get_maximum( PsiKPiAmplitudeModel& model, const unsigned int nburn ){ 
  
  double result = 0; 
  
  double mpipi, cospsi, cospipi, phi ;

  for ( unsigned int i =0 ; i < nburn; ++i ){
    
    mpipi   = gRandom->Uniform( Event::mPion + Event::mPion, Event::mB - Event::mPsi );
    cospipi = gRandom->Uniform( -1, 1 );
    cospsi = gRandom->Uniform( -1, 1 );
    phi    = gRandom->Uniform( -TMath::Pi(), TMath::Pi() );
   
    Event evt( 511, mpipi, cospipi, cospsi, phi );
    
    double value = model.evaluateCP( evt );
    
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
  
  //std::vector<std::string> resonancevector = {"f0(980)","f0(1500)","f2(1270)","f0(1790)","f2'(1525)"};
  std::vector<std::string> resonancevector = model.resonancelist();

//Scale amplitudes - iterate a number of times
  std::map<std::string,double> targetfracmap = get_target_fractions(parser.modelname);
  std::map<std::string,double> fitfracmap = model.fitfractions(true);
  model.print();
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

  model.rescaleres(resonancevector[0],fitfracmap,targetfracmap);
  model.print();
  std::map<std::string,double> newfitfracmap = model.fitfractions(true);

  double sumtarg{0},sumsec{0};
  for(auto z : resonancevector){


	  double targ = targetfracmap[z] / targetfracmap[resonancevector[0]];
	  double sec = newfitfracmap[z] / newfitfracmap[resonancevector[0]];

          sumtarg += targetfracmap[z];
          sumsec += newfitfracmap[z];


          std::cout << "For resonance " << z << "The ratios targ and sec are " << targ  << " "<<sec <<  std::endl;
  }
  std::cout << "Sumtarg, Sumsec" << sumtarg << " " << sumsec << std::endl;

  TFile* file = new TFile( parser.filename.c_str() ,"RECREATE");
  file->cd();
  
  double wgt; 
  double mpipi, mz;
  double cospsi;
  double cospipi;
  double phi;
  TH1D* mpipiHist = new TH1D("mpipiHist","",100,Event::mPion + Event::mPion,Event::mB - Event::mPsi);
  TH1D* cospsiHist = new TH1D("cospsiHist","",100,-1,1);
  TH1D* cospipiHist = new TH1D("cospipiHist","",100,-1,1);
  TH1D* phiHist = new TH1D("phiHist","",100,-TMath::Pi(),TMath::Pi());

  TStopwatch timer;
  timer.Start();
 
  
  double max_wgt = get_maximum( model, 100000 );
  
  for ( unsigned int i = 0; i < parser.nsample ; ++i ){
    
    do { 
      mpipi   = gRandom->Uniform( Event::mPion + Event::mPion, 
				 Event::mB - Event::mPsi );
      cospipi = gRandom->Uniform( -1, 1 );
      cospsi = gRandom->Uniform( -1, 1 );
      phi    = gRandom->Uniform( -TMath::Pi(), TMath::Pi() );
      
      Event evt(511, mpipi,cospipi,cospsi,phi);
      
      mz     = evt.mz_;

      if(parser.resonance.length() > 0){
	wgt    = model.evaluateCP( parser.resonance,evt ) ;

      }
      else{
        wgt    = model.evaluateCP( evt ) ;

      }
   

    } while ( gRandom->Rndm() > wgt/max_wgt ) ;
    
    
    mpipiHist->Fill(mpipi);
    cospsiHist->Fill(cospsi);
    cospipiHist->Fill(cospipi);
    phiHist->Fill(phi);


  }
  
  timer.Stop();
  
  std::cout << " time elapsed = " << timer.CpuTime() << std::endl;
 
if(parser.resonance.length() > 0){

	mpipiHist->Scale(newfitfracmap[parser.resonance]);
	phiHist->Scale(newfitfracmap[parser.resonance]);
	cospsiHist->Scale(newfitfracmap[parser.resonance]);
	cospipiHist->Scale(newfitfracmap[parser.resonance]);


}

  mpipiHist->Write(); 
  cospsiHist->Write();
  cospipiHist->Write();
  phiHist->Write();

  file->Close();

  return 0;
}
