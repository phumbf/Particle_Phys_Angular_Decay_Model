/*
 *Test script to reweight MC using the full rare mode decay model
 * */

#include "../include/ResonanceBase.h"
#include "../include/VectorResonance.h"
#include "../include/ScalarResonance.h"
#include "../include/ResonanceParser.h"
#include "../include/Amplitudes.h"
#include "../include/Parameters.h" 
#include "../include/FFcalculation.h" 
#include "../include/FFparser.h" 
#include "../include/StatTools.h" 
#include "../include/RandomSampling.h"
#include "../include/generator.h"

#include <memory> 

#include "TGraph.h" 
#include "TMath.h" 
#include "TCanvas.h" 
#include "TAxis.h" 
#include "TH1D.h" 
#include "TFile.h" 
#include "TTree.h" 
#include "TLorentzVector.h" 

#include <iostream>
#include <fstream>

std::map< std::pair<Int_t,Long64_t>, std::vector<double> >  gettruthvals(const TString truthFile,
									 const TString track1,
									 const TString track2){

	std::map< std::pair<Int_t,Long64_t>, std::vector<double> > truthmap;
	//Vector should include mhh,qsq,cosl,cosh,phi
	std::vector<double> phasespacevec;
	phasespacevec.resize(10);

	TFile *truthfile = TFile::Open(truthFile);
	TTree *truthtree = dynamic_cast<TTree*>(truthfile->Get("MCDecayTree"));

	Double_t tr1PE{0},tr1PX{0},tr1PY{0},tr1PZ{0};
	Double_t tr2PE{0},tr2PX{0},tr2PY{0},tr2PZ{0};
	Double_t tr3PE{0},tr3PX{0},tr3PY{0},tr3PZ{0};
	Double_t tr4PE{0},tr4PX{0},tr4PY{0},tr4PZ{0};
	UInt_t runNumber{0};
	ULong64_t eventNumber{0};

	Double_t Cos_h{0}, Cos_l{0}, Phi{0};

	Double_t t1ID{0}, t2ID{0}, t3ID{0}, t4ID{0};
	Double_t BID{0};

	truthtree->SetBranchAddress("runNumber",&runNumber);
	truthtree->SetBranchAddress("eventNumber",&eventNumber);
	truthtree->SetBranchAddress(track1 + "_TRUEP_E",&tr1PE);
	truthtree->SetBranchAddress(track1 + "_TRUEP_X",&tr1PX);
	truthtree->SetBranchAddress(track1 + "_TRUEP_Y",&tr1PY);
	truthtree->SetBranchAddress(track1 + "_TRUEP_Z",&tr1PZ);
	truthtree->SetBranchAddress(track2 + "_TRUEP_E",&tr2PE);
	truthtree->SetBranchAddress(track2 + "_TRUEP_X",&tr2PX);
	truthtree->SetBranchAddress(track2 + "_TRUEP_Y",&tr2PY);
	truthtree->SetBranchAddress(track2 + "_TRUEP_Z",&tr2PZ);
	truthtree->SetBranchAddress("muplus_TRUEP_E",&tr3PE);
	truthtree->SetBranchAddress("muplus_TRUEP_X",&tr3PX);
	truthtree->SetBranchAddress("muplus_TRUEP_Y",&tr3PY);
	truthtree->SetBranchAddress("muplus_TRUEP_Z",&tr3PZ);
	truthtree->SetBranchAddress("muminus_TRUEP_E",&tr4PE);
	truthtree->SetBranchAddress("muminus_TRUEP_X",&tr4PX);
	truthtree->SetBranchAddress("muminus_TRUEP_Y",&tr4PY);
	truthtree->SetBranchAddress("muminus_TRUEP_Z",&tr4PZ);
	truthtree->SetBranchAddress("Cos_h",&Cos_h);
	truthtree->SetBranchAddress("Cos_l",&Cos_l);
	truthtree->SetBranchAddress("Phi",&Phi);
	truthtree->SetBranchAddress(track1 + "_MCID",&t1ID);
	truthtree->SetBranchAddress(track2 + "_MCID",&t2ID);
	truthtree->SetBranchAddress("muplus_MCID",&t3ID);
	truthtree->SetBranchAddress("muminus_MCID",&t4ID);

	TLorentzVector tr1, tr2, tr3, tr4, hh, mumu, B;

	for(Int_t i=0; i<truthtree->GetEntries(); i++){

		truthtree->GetEntry(i);

		tr1.SetPxPyPzE(tr1PX,tr1PY,tr1PZ,tr1PE);
		tr2.SetPxPyPzE(tr2PX,tr2PY,tr2PZ,tr2PE);
		tr3.SetPxPyPzE(tr3PX,tr3PY,tr3PZ,tr3PE);
		tr4.SetPxPyPzE(tr4PX,tr4PY,tr4PZ,tr4PE);

		hh = tr1 + tr2;
		mumu = tr3 + tr4;
		B = hh + mumu;

		double qsq = mumu.M()*mumu.M()*1E-6;

		phasespacevec[0] = qsq;
		phasespacevec[1] = hh.M()*1E-3;
	    	phasespacevec[2] = Cos_h;	
		phasespacevec[3] = Cos_l;
		phasespacevec[4] = Phi;
		phasespacevec[5] = t1ID;
		phasespacevec[6] = t2ID;
		phasespacevec[7] = t3ID;
		phasespacevec[8] = t4ID;
		phasespacevec[9] = BID;

		truthmap[std::make_pair(runNumber,eventNumber)] = phasespacevec;
	}

	return truthmap;
}

void reweightMC(const TString MCfilename,
		const std::string fullname,
		const std::string decname,
		const TString truthfile,
		const TString track1,
		const TString track2,
		const TString truthtrack1,
		const TString truthtrack2,
		Double_t normalisation){

  	TFile *MCFile = TFile::Open(MCfilename + ".root");
	TTree *MCTree = dynamic_cast<TTree*>(MCFile->Get("DecayTree"));
  	TFile *newMCFile = new TFile(MCfilename + "w.root","RECREATE");
  	TTree *newMCTree = MCTree->CloneTree(0); 

	Long64_t num = MCTree->GetEntries();

	Double_t mhh{0}, cosh{0}, cosl{0}, phi{0}, jpsim{0};
	Int_t t1ID{0}, t2ID{0}, t3ID{0}, t4ID{0};
	Int_t t1TID{0}, t2TID{0}, t3TID{0}, t4TID{0};
	Int_t BID{0}, BTID{0};

	Bool_t pipisMuon{0}, pimisMuon{0};

	UInt_t runNumber{0};
	ULong64_t eventNumber{0};
	MCTree->SetBranchAddress("hh_M",&mhh);
	MCTree->SetBranchAddress("Cos_h",&cosh);
	MCTree->SetBranchAddress("Cos_l",&cosl);
	MCTree->SetBranchAddress("Phi",&phi);
	MCTree->SetBranchAddress("Jpsi_M",&jpsim);
	MCTree->SetBranchAddress("runNumber",&runNumber);
	MCTree->SetBranchAddress("eventNumber",&eventNumber);
	MCTree->SetBranchAddress(track1 + "_ID",&t1ID);
	MCTree->SetBranchAddress(track2 + "_ID",&t2ID);
	MCTree->SetBranchAddress("muplus_ID",&t3ID);
	MCTree->SetBranchAddress("muminus_ID",&t4ID);
	MCTree->SetBranchAddress(track1 + "_TRUEID",&t1TID);
	MCTree->SetBranchAddress(track2 + "_TRUEID",&t2TID);
	MCTree->SetBranchAddress("muplus_TRUEID",&t3TID);
	MCTree->SetBranchAddress("muminus_TRUEID",&t4TID);
	MCTree->SetBranchAddress("B0_ID",&BID);
	MCTree->SetBranchAddress("B0_TRUEID",&BTID);
	MCTree->SetBranchAddress(track1 + "_isMuon",&pipisMuon);
	MCTree->SetBranchAddress(track2 + "_isMuon",&pimisMuon);

	Double_t weight{0};
	newMCTree->Branch("wgtratio",&weight,"wgtratio/D");

	//Amplitudes - fullmodel
	Resparser fullparser;	
	std::vector<std::unique_ptr<ResonanceBase>> fulllist = fullparser.getResonances(fullname);
	Amplitudes fullamp;
	for(auto &r : fulllist){
		fullamp.addResonance(r);
	}
	//Amplitudes - DecFilemodel
	Resparser decparser;
	std::vector<std::unique_ptr<ResonanceBase>> declist = decparser.getResonances(decname);
	Amplitudes decamp;
	for(auto &r : declist){
		decamp.addResonance(r);
	}
	
	//Obtain truthmap
	std::map< std::pair<Int_t,Long64_t>, std::vector<double> > truthmap =  gettruthvals(truthfile+".root",truthtrack1,truthtrack2);

    //Needs to happen twice: Once for normalisation 
    double weightsum{0};
    
    //Loop through the MC 
    for(Long64_t i=0; i<num; i++){

        MCTree->GetEntry(i);

        std::vector<double> evt;

        auto it = truthmap.find( std::make_pair( runNumber, eventNumber ) );
        if ( it == truthmap.end() ){ 
            std::cout << "Empty vec" << std::endl;
            std::cout << "RUNNUMBER: " << runNumber << " EVTNUMBER: " << eventNumber << std::endl;
            continue;
        }

        else{
            evt = it->second;
        }

        //Calculate wgt for each model
        double fullvalue{0};
        double decvalue{0};

        //Check for pionisMuon and also effectively ignore events with qsq > 19.0GeV
        if(evt[0] > 19.0){
            weight = 1;
        }
        else{

            fullvalue = fullamp(evt);
            decvalue = decamp(evt);
            weight = fullvalue/decvalue;

            if(decvalue == 0 || isnan(weight)){
                weight = 1;
            }
        }
        weight/=normalisation;
        newMCTree->Fill();
    }
    	newMCTree->AutoSave();

	MCFile->Close();
	newMCFile->Close();

	return;
}

void reweight(const TString MCfilename,
              const TString Truthfilename,
            const std::string fullname,
            const std::string decname,
            const TString track1,
            const TString track2,
            const TString truthtrack1,
            const TString truthtrack2){
        

    double normalisation = reweightTRUTH(Truthfilename,
            fullname,
            decname,
            truthtrack1,
            truthtrack2);


    reweightMC(MCfilename,
            fullname,
            decname,
            Truthfilename,
            track1,
            track2,
            truthtrack1,
            truthtrack2,
            normalisation);

    return; 

}
