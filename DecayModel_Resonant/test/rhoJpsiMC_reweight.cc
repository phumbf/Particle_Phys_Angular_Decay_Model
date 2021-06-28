//Reweight the rhoJpsi MC using the ratio of the pipimodel with that amplitude and the full pipimodel amplitude.
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
  
  PsiKPiAmplitudeModel pipimodel,Decmodel;
  pipimodel.initialise("../config/B2rhoJpsi_from_2014_012.json");
  Decmodel.initialise("../config/DecFilerho.json");
  
  TFile *DataFile = TFile::Open(inputfile);
  TTree *DataTree = dynamic_cast<TTree*>( DataFile->Get(inputtreepath) );
  TFile *newDataFile = new TFile(inputfile + "w.root","RECREATE");
  TTree *newDataTree = DataTree->CloneTree(0); //0 argument means don't copy events
  
  double wgt; 
  double wgtrhojpsi;
  double wgtratio;
  double mpipi;
  double cospsi;
  double cospipi;
  double phi;

  DataTree->SetBranchAddress("Cos_pion1",&cospipi);
  DataTree->SetBranchAddress("Cos_muon1",&cospsi);
  DataTree->SetBranchAddress("Phi",&phi);
  DataTree->SetBranchAddress(resonance,&mpipi);

  TBranch * newb = newDataTree->Branch("wgtratio",&wgtratio,"wgtratio/D");

  TStopwatch timer;
  timer.Start();
 
  double wgtratiosum{0};

  for ( Int_t i=0; i<DataTree->GetEntries(); i++){

      DataTree->GetEntry(i); 
      mpipi = mpipi / 1000;

      std::cout << "mpipi, cospipi, cospsi, phi" << mpipi << " " << cospipi << " " << cospsi << " " << phi << std::endl;

      Event evt(511,mpipi,cospipi,cospsi,phi);
      
      wgt    = pipimodel.evaluate( evt ) ;
      std::cout << "wgt is : " << wgt << std::endl;
      wgtrhojpsi  = Decmodel.evaluate( evt ) ;
      std::cout << "wgtrhojpsi is : " << wgtrhojpsi << std::endl;
      wgtratio = wgt / wgtrhojpsi;
      if(wgt ==0 && wgtrhojpsi ==0){
		wgtratio = 0;
      }
      wgtratiosum += wgtratio;
      std::cout << "wgtratio is : " << wgtratio << std::endl;
    
  }
  std::cout << "--------------------------------------------------FIN---------------------------------" << std::endl;
  double normfactor = wgtratiosum / (DataTree->GetEntries());

  std::cout << "Norm factor is:  " << normfactor << std::endl;
  std::cout << "wgtratiosum:  " << wgtratiosum << std::endl;
  std::cout << "DataTree->GetEntries:  " << DataTree->GetEntries() << std::endl;

  for ( Int_t i=0; i<DataTree->GetEntries(); i++){

      DataTree->GetEntry(i); 
      mpipi = mpipi / 1000;

      std::cout << "mpipi, cospipi, cospsi, phi" << mpipi << " " << cospipi << " " << cospsi << " " << phi << std::endl;

      Event evt(511,mpipi,cospipi,cospsi,phi);
      
      wgt    = pipimodel.evaluate( evt ) ;
      wgtrhojpsi  = Decmodel.evaluate( evt ) ;
      wgtratio = (wgt / wgtrhojpsi) / normfactor;
      if(wgt ==0 && wgtrhojpsi ==0){
		wgtratio = 0;
      }
    
      std::cout << "wgtratio is : " << wgtratio << std::endl;
    newDataTree->Fill(); 
  }

  newDataTree->AutoSave();
  
  timer.Stop();
  
  std::cout << " time elapsed = " << timer.CpuTime() << std::endl;

  DataFile->Close();
  newDataFile->Close();

  return 0;
}
