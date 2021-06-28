//Reweight the K*Jpsi MC using the ratio of the bellemodel with that amplitude and the full bellemodel amplitude.
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

int main( int argc, char* argv[] ) {
  
  TString inputfile = argv[1];
  TString inputtreepath = argv[2];
  TString resonance  = argv[3];
  
  PsiKPiAmplitudeModel bellemodel,Decmodel;
  bellemodel.initialise( "../config/bellemodel.json");
  Decmodel.initialise( "../config/DecFileKpi.json");
  
  TFile *DataFile = TFile::Open(inputfile);
  TTree *DataTree = dynamic_cast<TTree*>( DataFile->Get(inputtreepath) );
  TFile *newDataFile = new TFile(inputfile + "w.root","RECREATE");
  TTree *newDataTree = DataTree->CloneTree(0); //0 argument means don't copy events
  
  double wgt; 
  double wgtkstjpsi;
  double wgtratio;
  double rawratio;
  double mpipi;
  double cospsi;
  double cospipi;
  double phi;

  DataTree->SetBranchAddress("Cos_kaon",&cospipi);
  DataTree->SetBranchAddress("Cos_muon1",&cospsi);
  DataTree->SetBranchAddress("Phi",&phi);
  DataTree->SetBranchAddress(resonance,&mpipi);

  TStopwatch timer;
  timer.Start();
 
  double wgtratiosum{0};

  for ( Int_t i=0; i<DataTree->GetEntries(); i++){

      DataTree->GetEntry(i); 
      mpipi = mpipi / 1000;

      Event evt(511,mpipi,cospipi,cospsi,phi);
      wgt    = bellemodel.evaluate( evt ) ;
      wgtkstjpsi  = Decmodel.evaluate( evt ) ;

      if(wgtkstjpsi ==0 || wgt == 0 ){
	rawratio = 0;
      }
      else{
      	rawratio = wgt / wgtkstjpsi;
      }
      wgtratiosum += rawratio;
    
  }

  double normfactor = wgtratiosum / (DataTree->GetEntries());

  TBranch * newb = newDataTree->Branch("wgtratio",&wgtratio,"wgtratio/D");

  for ( Int_t i=0; i<DataTree->GetEntries(); i++){

      DataTree->GetEntry(i); 
      mpipi = mpipi / 1000;

      Event evt(511,mpipi,cospipi,cospsi,phi);
      
      wgt    = bellemodel.evaluate( evt ) ;
      wgtkstjpsi  = Decmodel.evaluate( evt ) ;
      if(wgtkstjpsi ==0 || wgt == 0 ){
	wgtratio = 0;
      }
      else{
      	wgtratio = wgt/wgtkstjpsi;
      	wgtratio /= normfactor;
      }
    
      newDataTree->Fill(); 
  }

  newDataTree->AutoSave();
  
  timer.Stop();
  
  std::cout << " time elapsed = " << timer.CpuTime() << std::endl;

  DataFile->Close();
  newDataFile->Close();

  return 0;
}
