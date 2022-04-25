#ifndef __CINT__
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


#include "sPhenixStyle.C"
#include "sPhenixStyle.h"

using namespace std;

using std::cout;
using std::endl;
#endif



void DrawTogetherCross(){



	gStyle->SetOptStat(0);
	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();



	const int NTypes = 2;
	const int NDir = 2;

	TString infile;
	TString infile2;

	TString OutFileName[NTypes]={"FVTXN","FVTXS"};
	TString DataName[NDir]={"FVTX North J/#psi","FVTX South J/#psi"};

	TLegend* leg[NTypes];

	float Min[NTypes] = {0,0};
	float Max[NTypes] = {18,18};


	TF1 * func = new TF1("Func","x",0,19);
	func->SetLineStyle(2);

	for(int i = 0; i < NTypes; i ++){


		leg[i] = new TLegend(0.17,0.68,0.50,0.88,NULL,"brNDC");
		leg[i]->SetBorderSize(0);
		leg[i]->SetTextSize(0.040);
		leg[i]->SetTextFont(42);
		leg[i]->SetFillStyle(0);
		leg[i]->SetLineWidth(3);






		infile = Form("/sphenix/user/zshi/PHENIX/JPsiAnaCodes/FinalResultsUnCorr/OutFiles/JPsiR_%d_%d.root",i,i);
		infile2 = Form("/sphenix/user/zshi/PHENIX/JPsiAnaCodes/FinalResults/OutFiles/JPsiR_%d_%d.root",i,i);

		TFile * fin = new TFile(infile.Data());
		TFile * fin2 = new TFile(infile2.Data());

		TH1D * RatioScaledCorr = (TH1D * ) fin->Get("RatioScaledCorr");
		TF1 * FitFunc = (TF1 * ) fin->Get("FitFunc");


		TH1D * RatioScaledCorr2 = (TH1D * ) fin2->Get("RatioScaledCorr");
		TF1 * FitFunc2 = (TF1 * ) fin2->Get("FitFunc");

		TGraphAsymmErrors *Syst = (TGraphAsymmErrors *) fin->Get("Graph");
		TGraphAsymmErrors *Syst2 = (TGraphAsymmErrors *) fin2->Get("Graph");

		RatioScaledCorr->SetMarkerColor(kGreen);
		RatioScaledCorr->SetLineColor(kGreen);
		FitFunc->SetLineColor(kGreen);


		Syst->SetFillColorAlpha(kGreen-9,0.5);
		Syst->SetLineColor(kGreen-9);


		RatioScaledCorr2->SetMarkerColor(kBlue);
		RatioScaledCorr2->SetLineColor(kBlue);
		FitFunc2->SetLineColor(kBlue);


		Syst2->SetFillColorAlpha(kBlue-9,0.5);
		Syst2->SetLineColor(kBlue-9);


		RatioScaledCorr->SetMinimum(Min[i]);
		RatioScaledCorr->SetMaximum(Max[i]);				
		RatioScaledCorr->SetTitle(DataName[i].Data());
		RatioScaledCorr->Draw("ep");


		RatioScaledCorr2->Draw("epSAME");


		FitFunc->Draw("SAME");
		FitFunc2->Draw("SAME");


		float p0 = FitFunc->GetParameter(0);
		float p1 = FitFunc->GetParameter(1);
		float p02 = FitFunc2->GetParameter(0);
		float p12 = FitFunc2->GetParameter(1);


		TString FuncForm = Form("Fit Function: y = %.2f + %.2f x",p0,p1);
		TString FuncForm2 = Form("Fit Function: y = %.2f + %.2f x",p02,p12);

		leg[i]->AddEntry(RatioScaledCorr,"Before Dimuon Subtraction","LP");
		leg[i]->AddEntry(FitFunc,FuncForm.Data(),"LP");
		leg[i]->AddEntry(RatioScaledCorr2,"After Dimuon Subtraction","LP");
		leg[i]->AddEntry(FitFunc2,FuncForm2.Data(),"LP");
		Syst->Draw("5SAME");
		Syst2->Draw("5SAME");




		leg[i]->Draw("SAME");
		func->Draw("SAME");

		c->SaveAs(Form("CrossPlots/%s.png",OutFileName[i].Data()));

	}








}




