/* Add helicity angles to file for use in Helicity Amplitude program*/

#include "TVector3.h"

void add(const TString& file, const TString& treename, const TString& track1, const TString& track2, const TString& track3, const TString& track4, const TString& hadres, const TString& psires){

	TFile *myFile = TFile::Open(file + ".root");
	TTree *myTree = dynamic_cast<TTree*>( myFile->Get(treename));

	TFile *newmyFile = new TFile(file + "H.root","RECREATE");
	TTree *newmyTree = myTree->CloneTree(0); //0 argument means don't copy events

	TString kaon_Branch = track1;
	TString pion_Branch = track2;
	TString muon1_Branch = track3;
	TString muon2_Branch = track4;
	TString Jpsi_Branch = psires;
	TString Kst_Branch = hadres;

	Double_t Pion_PX{0}, Pion_PY{0}, Pion_PZ{0}, Pion_PE{0};
	Double_t Kaon_PX{0}, Kaon_PY{0}, Kaon_PZ{0}, Kaon_PE{0};
	Double_t Jpsi_PX{0}, Jpsi_PY{0}, Jpsi_PZ{0}, Jpsi_PE{0};
	Double_t Kst_PX{0}, Kst_PY{0}, Kst_PZ{0}, Kst_PE{0};
	Double_t muon1_PX{0}, muon1_PY{0}, muon1_PZ{0},	muon1_PE{0};	
	Double_t muon2_PX{0}, muon2_PY{0}, muon2_PZ{0}, muon2_PE{0};
	Int_t Kaon_ID{0};

	myTree->SetBranchAddress(pion_Branch + "_PX",&Pion_PX);
	myTree->SetBranchAddress(pion_Branch + "_PY",&Pion_PY);
	myTree->SetBranchAddress(pion_Branch + "_PZ",&Pion_PZ);
	myTree->SetBranchAddress(pion_Branch + "_PE",&Pion_PE);
	myTree->SetBranchAddress(kaon_Branch + "_PX",&Kaon_PX);
	myTree->SetBranchAddress(kaon_Branch + "_PY",&Kaon_PY);
	myTree->SetBranchAddress(kaon_Branch + "_PZ",&Kaon_PZ);
	myTree->SetBranchAddress(kaon_Branch + "_PE",&Kaon_PE);
	myTree->SetBranchAddress(Jpsi_Branch + "_PX",&Jpsi_PX);
	myTree->SetBranchAddress(Jpsi_Branch + "_PY",&Jpsi_PY);
	myTree->SetBranchAddress(Jpsi_Branch + "_PZ",&Jpsi_PZ);
	myTree->SetBranchAddress(Jpsi_Branch + "_PE",&Jpsi_PE);
	myTree->SetBranchAddress(Kst_Branch + "_PX",&Kst_PX);
	myTree->SetBranchAddress(Kst_Branch + "_PY",&Kst_PY);
	myTree->SetBranchAddress(Kst_Branch + "_PZ",&Kst_PZ);
	myTree->SetBranchAddress(Kst_Branch + "_PE",&Kst_PE);
	myTree->SetBranchAddress(muon1_Branch + "_PX",&muon1_PX);
	myTree->SetBranchAddress(muon1_Branch + "_PY",&muon1_PY);	
	myTree->SetBranchAddress(muon1_Branch + "_PZ",&muon1_PZ);	
	myTree->SetBranchAddress(muon1_Branch + "_PE",&muon1_PE);	
	myTree->SetBranchAddress(muon2_Branch + "_PX",&muon2_PX);
	myTree->SetBranchAddress(muon2_Branch + "_PY",&muon2_PY);
	myTree->SetBranchAddress(muon2_Branch + "_PZ",&muon2_PZ);
	myTree->SetBranchAddress(muon2_Branch + "_PE",&muon2_PE);
	myTree->SetBranchAddress("Kaon_ID",&Kaon_ID);

	Double_t Cos_kaon{0}, Cos_muon1{0}, Phi{0}, kaonpion_M{0}, Psipion_M{0};

	TBranch * newb = newmyTree->Branch("Cos_kaon",&Cos_kaon,"Cos_kaon/D");
	TBranch * newb2 = newmyTree->Branch("Cos_muon1",&Cos_muon1,"Cos_muon1/D");
	TBranch * newb3 = newmyTree->Branch("Phi",&Phi,"Phi/D");
	TBranch * newb4 = newmyTree->Branch("kaonpion_M",&kaonpion_M,"kaonpion_M/D");
	TBranch * newb5 = newmyTree->Branch("Psipion_M",&Psipion_M,"Psipion_M/D");

	TLorentzVector Pion, Kaon, Jpsi, Kst, muon1, muon2, B;
	TLorentzVector kaonpion, Psipion;

	Int_t num = myTree->GetEntries();

	TH1D *coskaonhist = new TH1D("coskaonhist","CosKaon",100,-1,1);	
	TH1D *cosmuonhist = new TH1D("cosmuonhist","CosMuon",100,-1,1);	
	TH1D *phihist = new TH1D("phihist","Phi",100,-4,4);	

	for(Int_t i=0; i<num; i++){

		myTree->GetEntry(i);

		Pion.SetPxPyPzE(Pion_PX,Pion_PY,Pion_PZ,Pion_PE);
		Kaon.SetPxPyPzE(Kaon_PX,Kaon_PY,Kaon_PZ,Kaon_PE);
		Jpsi.SetPxPyPzE(Jpsi_PX,Jpsi_PY,Jpsi_PZ,Jpsi_PE);
		Kst.SetPxPyPzE(Kst_PX,Kst_PY,Kst_PZ,Kst_PE);
		muon1.SetPxPyPzE(muon1_PX,muon1_PY,muon1_PZ,muon1_PE);
		muon2.SetPxPyPzE(muon2_PX,muon2_PY,muon2_PZ,muon2_PE);
		
		B = Kst + Jpsi;
		kaonpion = Kaon + Pion;
		kaonpion_M = kaonpion.M();
		Psipion = Jpsi + Pion;
		Psipion_M = Psipion.M();

		/////// Calculation of helicity angles ////////
		//--------------------Kaon Angle-------------//	
		//B frame boost vector
		const TVector3 B_frameBoost = B.BoostVector();

		//Boost Jpsi into B frame
		TLorentzVector Jpsi_Bframe = Jpsi;
		Jpsi_Bframe.Boost(-1*B_frameBoost);

		//Boost Kst into B frame
		TLorentzVector Kst_Bframe = Kst;
		Kst_Bframe.Boost(-1*B_frameBoost);

		//Direction of Kst in the B frame
		const TVector3 Kst_frameBoost = Kst_Bframe.BoostVector();

		//Boost the Kaon into the Kst frame
		TLorentzVector Kaon_Kstframe = Kaon;
		Kaon_Kstframe.Boost(-1*B_frameBoost);
		Kaon_Kstframe.Boost(-1*Kst_frameBoost);
		//Get the Kaon vector and find first angle
		TVector3 Kaon_KstframeVec = Kaon_Kstframe.Vect();
		TVector3 Kst_BframeVec = Kst_Bframe.Vect();
		Double_t angle_kaon = Kst_BframeVec.Angle(Kaon_KstframeVec);
		//--------------------------------------------//

		//--------------------Muon Angle--------------//	
		//Direction of Jpsi in the B frame
		const TVector3 Jpsi_frameBoost = Jpsi_Bframe.BoostVector();

		//Boost the muon1 into the Jpsi frame
		TLorentzVector muon1_Jpsiframe = muon1;
		muon1_Jpsiframe.Boost(-1*B_frameBoost);
		muon1_Jpsiframe.Boost(-1*Jpsi_frameBoost);
		//Get the muon1 vector and find next angle
		TVector3 muon1_JpsiframeVec = muon1_Jpsiframe.Vect();
		TVector3 Jpsi_BframeVec = Jpsi_Bframe.Vect();
		Double_t angle_muon1 = Jpsi_BframeVec.Angle(muon1_JpsiframeVec);
		//--------------------------------------------//

		//Take cosines of Kaon and Muon angle
		Cos_kaon = TMath::Cos(angle_kaon);
		Cos_muon1 = TMath::Cos(angle_muon1);

		coskaonhist->Fill(Cos_kaon);
		cosmuonhist->Fill(Cos_muon1);

		//--------------------Phi (angle between Kst and Jpsi decay planes in B frame) --------------//	
		//Determine normal vectors of the two planes in the B frame
		TLorentzVector Kaon_Bframe = Kaon;
		Kaon_Bframe.Boost(-1*B_frameBoost);
		TLorentzVector Pion_Bframe = Pion;
		Pion_Bframe.Boost(-1*B_frameBoost);
		TVector3 Kaon_BframeVec = Kaon_Bframe.Vect();
		TVector3 Pion_BframeVec = Pion_Bframe.Vect();
		TVector3 Kst_Plane = Kaon_BframeVec.Cross(Pion_BframeVec);

		TLorentzVector muon1_Bframe = muon1;
		muon1_Bframe.Boost(-1*B_frameBoost);
		TLorentzVector muon2_Bframe = muon2;
		muon2_Bframe.Boost(-1*B_frameBoost);
		TVector3 muon1_BframeVec = muon1_Bframe.Vect();
		TVector3 muon2_BframeVec = muon2_Bframe.Vect();
		TVector3 Jpsi_Plane = muon1_BframeVec.Cross(muon2_BframeVec);

		Kst_Plane = Kst_Plane.Unit();
		Jpsi_Plane = Jpsi_Plane.Unit();

		//Determine angle between the two normal vectors
		Double_t cosphi = Kst_Plane.Dot(Jpsi_Plane);
		Double_t sinphi = ((Kst_Plane.Cross(Jpsi_Plane))).Dot(Kst_BframeVec.Unit());

		Phi = TMath::ATan2(sinphi,cosphi);
		if(Kaon_ID < 0){
			Phi = -Phi;
		}

		phihist->Fill(Phi);

		newmyTree->Fill();
	}

	newmyTree->AutoSave();
	myFile->Close();
	newmyFile->Close();
}

