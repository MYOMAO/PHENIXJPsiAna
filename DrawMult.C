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




//#include "his.h"
using namespace std;

using std::cout;
using std::endl;

void DrawMult(){


	gStyle->SetOptStat(0);


	const int NBins = 5;
	int MultBin[NBins + 1] = {0,2,5,8,12,19};

	TCanvas *c = new TCanvas("c","c",600,600);
	c->cd();

	c->SetLogy();

	TFile * fin = new TFile("merged_anaT2_run15pp_MU_v14.root");
	fin->cd();

	TTree * anaT2 = (TTree * ) fin->Get("anaT2");

	TH1D * Evt_Mult_FVTXN = new TH1D("Evt_Mult_FVTXN","",21,-0.5,20.5);
	Evt_Mult_FVTXN->GetXaxis()->SetTitle("Evt_Mult_FVTXN");
	Evt_Mult_FVTXN->GetYaxis()->SetTitle("Counts");
	Evt_Mult_FVTXN->GetYaxis()->SetTitleOffset(1.4);
	Evt_Mult_FVTXN->GetXaxis()->CenterTitle();
	Evt_Mult_FVTXN->GetYaxis()->CenterTitle();
	Evt_Mult_FVTXN->SetMaximum(20000);
	Evt_Mult_FVTXN->SetLineColor(kBlue);
	Evt_Mult_FVTXN->SetLineWidth(2);
	Evt_Mult_FVTXN->SetTitle("FVTX North Multiplicity Distribution");

	TH1D * Evt_Mult_FVTXS = new TH1D("Evt_Mult_FVTXS","",21,-0.5,20.5);
	Evt_Mult_FVTXS->GetXaxis()->SetTitle("Evt_Mult_FVTXS");
	Evt_Mult_FVTXS->GetYaxis()->SetTitle("Counts");
	Evt_Mult_FVTXS->GetYaxis()->SetTitleOffset(1.4);
	Evt_Mult_FVTXS->GetXaxis()->CenterTitle();
	Evt_Mult_FVTXS->GetYaxis()->CenterTitle();
	Evt_Mult_FVTXS->SetMaximum(20000);
	Evt_Mult_FVTXS->SetLineColor(kBlue);
	Evt_Mult_FVTXS->SetLineWidth(2);
	Evt_Mult_FVTXS->SetTitle("FVTX South Multiplicity Distribution");



	anaT2->Project("Evt_Mult_FVTXN","Evt_Mult_FVTXN");



	Evt_Mult_FVTXN->Draw("hist");


	TLine * l[NBins];

	for(int i = 0; i < NBins; i ++){
		l[i] =  new TLine(MultBin[i+1],1,MultBin[i+1],20000);
		l[i]->SetLineWidth(2);
		l[i]->SetLineStyle(2);
		l[i]->SetLineColor(kRed);
		l[i]->Draw("SAME");
	}


	c->SaveAs("MultPlots/Evt_Mult_FVTXN.png");


	anaT2->Project("Evt_Mult_FVTXS","Evt_Mult_FVTXS");
	Evt_Mult_FVTXS->Draw("hist");

	for(int i = 0; i < NBins; i++){
		l[i]->Draw("SAME");
	}


	c->SaveAs("MultPlots/Evt_Mult_FVTXS.png");

}



