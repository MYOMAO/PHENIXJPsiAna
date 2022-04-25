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

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TGraph.h"

#include "TLatex.h"
#include "TStyle.h"

#include "TLegend.h"



void PreScale(int Opt){


	TString infile;


	if(Opt == 0) infile = "NorthMB.root";
	if(Opt == 1) infile = "SouthMB.root";


	TFile * fin = new TFile(infile.Data());
	fin->cd();

	TTree * anaT0 = (TTree *)  fin->Get("anaT0");


	int MB0_Scale_Down;
	int MB1_Scale_Down;
	int MB2_Scale_Down;
	float Evt_bbcZ;

	Float_t Evt_Mult_FVTXN;
	Float_t Evt_Mult_FVTXS;
	Float_t Evt_Mult_SVX;


	anaT0->SetBranchAddress("MB0_Scale_Down",&MB0_Scale_Down);
	anaT0->SetBranchAddress("MB1_Scale_Down",&MB1_Scale_Down);
	anaT0->SetBranchAddress("MB2_Scale_Down",&MB2_Scale_Down);
	anaT0->SetBranchAddress("Evt_bbcZ",&Evt_bbcZ);

	anaT0->SetBranchAddress("Evt_Mult_FVTXN",&Evt_Mult_FVTXN);
	anaT0->SetBranchAddress("Evt_Mult_FVTXS",&Evt_Mult_FVTXS);
	anaT0->SetBranchAddress("Evt_Mult_SVX",&Evt_Mult_SVX);



	/*
	   const int NBins = 5;
	   int MultBin[NBins + 1] = {0,2,5,8,12,19};
	   double MultBinHis[NBins + 1] = {0,2,5,8,12,19};
	   */

	const int NBins = 10;
	int MultBin[NBins + 1] = {0,1,2,3,4,5,6,8,10,12,19};
	double MultBinHis[NBins + 1] = {0,1,2,3,4,5,6,8,10,12,19};

	TH1D * Evt_Mult_FVTXN_His = new TH1D("Evt_Mult_FVTXN_His","",NBins,MultBinHis);
	TH1D * Evt_Mult_FVTXS_His = new TH1D("Evt_Mult_FVTXS_His","",NBins,MultBinHis);
	TH1D * Evt_Mult_SVX_His = new TH1D("Evt_Mult_SVX_His","",NBins,MultBinHis);

	int NEvents = anaT0->GetEntries();

	for(int i = 0; i < NEvents; i++){

		if(i%1000000 == 0)	cout << "Now Working on " << i << endl;

		anaT0->GetEntry(i);

		//cout << "Evt_bbcZ = " << Evt_bbcZ << "   abs(Evt_bbcZ)  = " << abs(Evt_bbcZ)  << endl;

		if(abs(Evt_bbcZ) < 10.0){
			Evt_Mult_FVTXN_His->Fill(Evt_Mult_FVTXN,MB2_Scale_Down + 1);
			Evt_Mult_FVTXS_His->Fill(Evt_Mult_FVTXS,MB2_Scale_Down + 1);
			Evt_Mult_SVX_His->Fill(Evt_Mult_SVX,MB2_Scale_Down + 1);

		}


	}



	TString OutFile;


	if(Opt == 0) OutFile = "NorthMBEvtMult.root";
	if(Opt == 1) OutFile = "SouthMBEvtMult.root";


	TFile * fout = new TFile(OutFile.Data(),"RECREATE");
	fout->cd();

	Evt_Mult_FVTXN_His->Write();
	Evt_Mult_FVTXS_His->Write();
	Evt_Mult_SVX_His->Write();

	fout->Close();




}
