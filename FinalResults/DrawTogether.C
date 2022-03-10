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



void DrawTogether(){



	gStyle->SetOptStat(0);
	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();



	const int NTypes = 3;
	const int NDir = 2;

	TString infile;

	TString OutFileName[NTypes]={"FVTXN","FVTXS","SVX"};
	TString DataName[NDir]={"FVTX North J/#psi","FVTX South J/#psi"};

	TLegend* leg[NTypes];

	float Min[NTypes] = {0,0,0};
	float Max[NTypes] = {18,24,12};

	
	TF1 * func = new TF1("Func","x",0,19);
	func->SetLineStyle(2);

	for(int i = 0; i < NTypes; i ++){


		leg[i] = new TLegend(0.17,0.68,0.50,0.88,NULL,"brNDC");
		leg[i]->SetBorderSize(0);
		leg[i]->SetTextSize(0.040);
		leg[i]->SetTextFont(42);
		leg[i]->SetFillStyle(0);
		leg[i]->SetLineWidth(3);


		for(int j = 0; j < NDir; j++){



			infile = Form("OutFiles/JPsiR_%d_%d.root",j,i);

			TFile * fin = new TFile(infile.Data());

			TH1D * RatioScaledCorr = (TH1D * ) fin->Get("RatioScaledCorr");
			TF1 * FitFunc = (TF1 * ) fin->Get("FitFunc");
				

			TGraphAsymmErrors *Syst = (TGraphAsymmErrors *) fin->Get("Graph");

			RatioScaledCorr->SetMarkerColor(j+3);
			RatioScaledCorr->SetLineColor(j+3);
			FitFunc->SetLineColor(j+3);
	
			if(j == 0){
				Syst->SetFillColorAlpha(kGreen-9,0.5);
				Syst->SetLineColor(kGreen-9);
			}
			if(j == 1){
				Syst->SetFillColorAlpha(kBlue-9,0.5);
				Syst->SetLineColor(kBlue-9);
			}

			if(j==0){
				RatioScaledCorr->SetMinimum(Min[i]);
				RatioScaledCorr->SetMaximum(Max[i]);				
				RatioScaledCorr->SetTitle("");
				RatioScaledCorr->Draw("ep");
			}
			
			if(j>0) RatioScaledCorr->Draw("epSAME");


			FitFunc->Draw("SAME");
			
	
			float p0 = FitFunc->GetParameter(0);
			float p1 = FitFunc->GetParameter(1);

			TString FuncForm = Form("Fit Function: y = %.2f + %.2f x",p0,p1);

			leg[i]->AddEntry(RatioScaledCorr,DataName[j].Data(),"LP");
			leg[i]->AddEntry(FitFunc,FuncForm.Data(),"LP");
			Syst->Draw("5SAME");

		}
		


		leg[i]->Draw("SAME");
		func->Draw("SAME");

		c->SaveAs(Form("FinalPlots/Together/%s.png",OutFileName[i].Data()));

	}








}




