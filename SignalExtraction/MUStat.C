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
//#include "tnp_weight_lowptPbPb.h"



//#include "his.h"
using namespace std;

using std::cout;
using std::endl;

void MUStat(int Opt, int MultOpt){

	gStyle->SetOptStat(0);
	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();

	
	TString inputdata;

	if(Opt == 0) inputdata = "CommonFiles/NorthMU.root";
	if(Opt == 1) inputdata = "CommonFiles/SouthMU.root";


	TString Var;
	
	if(MultOpt == 0) Var = "Evt_Mult_FVTXN";
	if(MultOpt == 1) Var = "Evt_Mult_FVTXS";
	if(MultOpt == 2) Var = "Evt_Mult_SVX";


	TFile * fin = new TFile(inputdata.Data());
	
	fin->cd();

	TTree * anaT2 = (TTree *) fin->Get("anaT2");

	
	const int NBins = 10;
	double MultBin[NBins + 1] = {0,1,2,3,4,5,6,8,10,12,19};

	TH1D * MUStatHis = new TH1D("MUStatHis","",NBins,MultBin);
	
	anaT2->Project("MUStatHis",Var.Data());
	
	
	MUStatHis->SetMarkerStyle(20);
	MUStatHis->SetMarkerColor(1);
	MUStatHis->SetLineColor(1);
	MUStatHis->SetMarkerSize(1);

	MUStatHis->SetTitle(Form("%s Multiplicity Distribution for MU Data",Var.Data()));
	MUStatHis->GetYaxis()->SetTitle("Counts");
	MUStatHis->GetYaxis()->SetTitleOffset(1.5);
	MUStatHis->GetXaxis()->SetTitle(Var.Data());
	MUStatHis->GetYaxis()->CenterTitle();
	MUStatHis->GetXaxis()->CenterTitle();
	
	MUStatHis->Draw("ep");

	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.042);
	lat->DrawLatex(0.50,0.75,"MU (Dimuon) Sample");

	c->SetLogy();
	c->SaveAs(Form("MUStatCheck/MUStat_%d_%d.png",Opt,MultOpt));


}
