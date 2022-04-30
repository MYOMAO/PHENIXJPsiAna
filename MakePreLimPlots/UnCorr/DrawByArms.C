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



void DrawByArms(){

	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.034);

	gStyle->SetOptStat(0);
	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();



	const int NTypes = 3;
	const int NDir = 2;

	TString infile;

	TString OutFileName[NTypes]={"FVTXN","FVTXS","SVX"};
	TString DataName[NDir]={"1.2 < y_{J/#psi} < 2.2","-2.2 < y_{J/#psi} < -1.2"};

	TLegend* leg[NTypes];

	float Min[NTypes] = {0,0,0};
	float Max[NTypes] = {20,20,20};

	int NUse = 2;





	leg[0] = new TLegend(0.12,0.65,0.50,0.80,NULL,"brNDC");
	leg[0]->SetBorderSize(0);
	leg[0]->SetTextSize(0.040);
	leg[0]->SetTextFont(42);
	leg[0]->SetFillStyle(0);
	leg[0]->SetLineWidth(3);




	for(int i = 0; i < NUse; i ++){






		infile = Form("OutFiles/JPsiR_%d_%d.root",i,i);

		TFile * fin = new TFile(infile.Data());

		TH1D * RatioScaledCorr = (TH1D * ) fin->Get("RatioScaledCorr");



		TGraphAsymmErrors *Syst = (TGraphAsymmErrors *) fin->Get("Graph");


		RatioScaledCorr->SetMarkerColor(i+3);
		RatioScaledCorr->SetLineColor(i+3);

		if(i==0){
			Syst->SetFillColorAlpha(kGreen-9,0.5);
			Syst->SetLineColor(kGreen-9);
		}
		if(i==1){
			Syst->SetFillColorAlpha(kBlue-9,0.5);
			Syst->SetLineColor(kBlue-9);
		}

		if(i==0){
			RatioScaledCorr->SetMinimum(Min[i]);
			RatioScaledCorr->SetMaximum(Max[i]);				
			RatioScaledCorr->SetTitle("");
			RatioScaledCorr->Draw("ep");
			RatioScaledCorr->GetXaxis()->SetTitle("N_{Ch}/<N_{Ch}>");
			
		}

		if(i==1) RatioScaledCorr->Draw("epSAME");


		leg[0]->AddEntry(RatioScaledCorr,DataName[i].Data(),"LP");
		Syst->Draw("5SAME");





	}


	leg[0]->Draw("SAME");


	lat->DrawLatex(0.15,0.61,"pp #sqrt{s} = 200 GeV");
	lat->DrawLatex(0.15,0.55,"N_{ch} and J/#psi: |#Delta y| = 0");
	lat->DrawLatex(0.15,0.49,"N_{ch}: 1.2 < |#eta| < 2.4");

	//		if(i == 2) lat->DrawLatex(0.15,0.57,"SVX: -1.0 < #eta < 1.0");

	c->SaveAs("Plots/ByArms/png/Same.png");
	c->SaveAs("Plots/ByArms/pdf/Same.pdf");


	leg[1] = new TLegend(0.12,0.65,0.50,0.80,NULL,"brNDC");
	leg[1]->SetBorderSize(0);
	leg[1]->SetTextSize(0.034);
	leg[1]->SetTextFont(42);
	leg[1]->SetFillStyle(0);
	leg[1]->SetLineWidth(3);


	for(int i = 0; i < NUse; i ++){





		infile = Form("OutFiles/JPsiR_%d_%d.root",i,1-i);

		TFile * fin = new TFile(infile.Data());

		TH1D * RatioScaledCorr = (TH1D * ) fin->Get("RatioScaledCorr");



		TGraphAsymmErrors *Syst = (TGraphAsymmErrors *) fin->Get("Graph");

		RatioScaledCorr->SetMarkerColor(i+3);
		RatioScaledCorr->SetLineColor(i+3);


		if(i==0){
			Syst->SetFillColorAlpha(kGreen-9,0.5);
			Syst->SetLineColor(kGreen-9);
		}
		if(i==1){
			Syst->SetFillColorAlpha(kBlue-9,0.5);
			Syst->SetLineColor(kBlue-9);
		}

		if(i==0){
			RatioScaledCorr->SetMinimum(Min[i]);
			RatioScaledCorr->SetMaximum(Max[i]);				
			RatioScaledCorr->SetTitle("");
			RatioScaledCorr->Draw("ep");
			RatioScaledCorr->GetXaxis()->SetTitle("N_{Ch}/<N_{Ch}>");
			
		}

		if(i==1) RatioScaledCorr->Draw("epSAME");


		leg[1]->AddEntry(RatioScaledCorr,DataName[i].Data(),"LP");
		Syst->Draw("5SAME");




	}


	leg[1]->Draw("SAME");


	lat->DrawLatex(0.15,0.61,"pp #sqrt{s} = 200 GeV");
	lat->DrawLatex(0.15,0.55,"N_{ch} and J/#psi: |#Delta y| > 1");
	lat->DrawLatex(0.15,0.49,"N_{ch}: 1.2 < |#eta| < 2.4");

	//		if(i == 2) lat->DrawLatex(0.15,0.57,"SVX: -1.0 < #eta < 1.0");

	c->SaveAs("Plots/ByArms/png/Opposite.png");
	c->SaveAs("Plots/ByArms/pdf/Opposite.pdf");



}




