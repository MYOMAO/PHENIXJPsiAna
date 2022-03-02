#include "TROOT.h"
#include "TH1.h"
#include "TTree.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TSystem.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TVector3.h"






#include "TRandom.h"
#include <iostream>
#include <fstream>


using namespace std;

using std::cout;
using std::endl;


void SkimMU(int Opt){

	TString FileName = "MU_ana_pp.root";

	TFile * fin = new TFile(FileName.Data());
	fin->cd();

	TTree * anaT2 = (TTree *) fin->Get("anaT2"); 


	TString RunFile;

	if(Opt == 0)  RunFile = "NorthRuns.txt";
	if(Opt == 1)  RunFile = "SouthRuns.txt";
	

	TString OutFileName;


	if(Opt == 0)  OutFileName = "NorthMU.root";
	if(Opt == 1)  OutFileName = "SouthMU.root";


	ifstream RunList(RunFile.Data());

	std::vector<int> GoodRun; 
	int RunID;


    while (!RunList.eof()){
		RunList >>	RunID;
		GoodRun.push_back(RunID);
	
//		cout << "RunID = " << RunID << endl;
	}


	int GoodRunSize = GoodRun.size();

	int RunNumber;
	int BBC_Rate;
	int Evt_Mult_FVTXN;
	int Evt_Mult_FVTXS;
	int Evt_Mult_FVTXN2;
	int Evt_Mult_FVTXS2;
	int Evt_Mult_SVX;
	int Evt_Mult_SVX2;

	int MB0_Scale_Down;
	int MB1_Scale_Down;
	int MB2_Scale_Down;
	int Evt_bbcZ;
	int Evt_vtxZ2;

	anaT2->SetBranchStatus("*",0);

	anaT2->SetBranchStatus("BBC_Rate",1);
	anaT2->SetBranchStatus("Evt_Mult_FVTXN",1);
	anaT2->SetBranchStatus("Evt_Mult_FVTXS",1);
	anaT2->SetBranchStatus("Evt_Mult_FVTXN2",1);
	anaT2->SetBranchStatus("Evt_Mult_FVTXS2",1);
	anaT2->SetBranchStatus("Evt_Mult_SVX",1);
	anaT2->SetBranchStatus("Evt_Mult_SVX2",1);
	anaT2->SetBranchStatus("Evt_Mult_SVX3",1);	
	anaT2->SetBranchStatus("Evt_Mult_FVTXN3",1);
	anaT2->SetBranchStatus("Evt_Mult_FVTXS3",1);
	anaT2->SetBranchStatus("Evt_bbcZ",1);
	anaT2->SetBranchStatus("MuID_N2D_Scale_Down",1);
	anaT2->SetBranchStatus("MuID_S2D_Scale_Down",1);
	anaT2->SetBranchStatus("dimuon_mass",1);
	anaT2->SetBranchStatus("dimuon_charge",1);
	anaT2->SetBranchStatus("dimuon_mass_fvtx",1);
	anaT2->SetBranchStatus("dimuon_mass_fvtxmutr",1);
	anaT2->SetBranchStatus("dimuon_dca_r",1);
	anaT2->SetBranchStatus("dimuon_dca_z",1);	
	anaT2->SetBranchStatus("dimuon_XYZ0",1);
	anaT2->SetBranchStatus("dimuon_XYZ0_fvtxmutr",1);
	anaT2->SetBranchStatus("dimuon_pXYZ",1);
	anaT2->SetBranchStatus("dimuon_pXYZ_fvtxmutr",1);
	anaT2->SetBranchStatus("dimuon_rapidity",1);
	
	anaT2->SetBranchStatus("RunNumber",1);

	
	anaT2->SetBranchAddress("RunNumber",&RunNumber);

	int NEvents = anaT2->GetEntries();
	
	int Pass;

	TFile * fout = new TFile(OutFileName.Data(),"RECREATE");

	fout->cd();

	TTree * anaT2_new = anaT2->CloneTree(0);


	for(int i = 0; i < NEvents; i++){

		if(i%10000==0) cout << "Now Working at Event " << i << endl;

		anaT2->GetEntry(i);
		
//		cout << "RunNumber = " << RunNumber << endl;
	
		Pass = 0;

		for(int j = 0; j < GoodRunSize; j++){


			if(RunNumber == GoodRun[j]) Pass = 1;
				

		}

		if(Pass == 1){

			anaT2_new->Fill();

		}
	

	}

	anaT2_new->Write();

	fout->Close();


}
