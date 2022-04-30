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

	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.034);

	gStyle->SetOptStat(0);
	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();

	c->SetBottomMargin(0.12);


	const int NTypes = 3;
	const int NDir = 2;

	TString infile;

	TString OutFileName[NTypes]={"FVTXN","FVTXS","SVX"};
	TString DataName[NDir]={"1.2 < y_{J/#psi} < 2.2","-2.2 < y_{J/#psi} < -1.2"};
	TString Together = "(J/#psi #rightarrow #mu^{+} #mu^{-} subtracted)";


	TLegend* leg[NTypes];

	float Min[NTypes] = {0,0,0};
	float Max[NTypes] = {20,20,20};


	for(int i = 0; i < NTypes; i ++){


		leg[i] = new TLegend(0.12,0.65,0.50,0.80,NULL,"brNDC");
		leg[i]->SetBorderSize(0);
		leg[i]->SetTextSize(0.040);
		leg[i]->SetTextFont(42);
		leg[i]->SetFillStyle(0);
		leg[i]->SetLineWidth(3);


		for(int j = 0; j < NDir; j++){



			infile = Form("OutFiles/JPsiR_%d_%d.root",j,i);

			TFile * fin = new TFile(infile.Data());

			TH1D * RatioScaledCorr = (TH1D * ) fin->Get("RatioScaledCorr");

				

			TGraphAsymmErrors *Syst = (TGraphAsymmErrors *) fin->Get("Graph");

			RatioScaledCorr->SetMarkerColor(j+3);
			RatioScaledCorr->SetLineColor(j+3);


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

				if(i == 0) RatioScaledCorr->GetXaxis()->SetTitle("N_{ch}^{N}/<N_{ch}^{N}>, #tilde{N}_{ch}^{N}/<#tilde{N}_{ch}^{N}>");
				if(i == 1) RatioScaledCorr->GetXaxis()->SetTitle("N_{ch}^{S}/<N_{ch}^{S}>, #tilde{N}_{ch}^{S}/<#tilde{N}_{ch}^{S}>");
				RatioScaledCorr->SetTitle("");
				RatioScaledCorr->Draw("ep");

			}
			
			if(j>0) RatioScaledCorr->Draw("epSAME");


			//if(i != j) leg[i]->AddEntry(RatioScaledCorr,DataName[j].Data(),"LP");
			//if(i == j) leg[i]->AddEntry(RatioScaledCorr,Form("%s %s",DataName[j].Data(),Together.Data()),"LP");
			leg[i]->AddEntry(RatioScaledCorr,DataName[j].Data(),"LP");	
			
			Syst->Draw("5SAME");

		}



		leg[i]->Draw("SAME");


		lat->DrawLatex(0.15,0.61,"pp #sqrt{s} = 200 GeV");
		if(i == 0) lat->DrawLatex(0.15,0.55,"N_{ch}^{N}: 1.2 < #eta < 2.4");
		if(i == 1) lat->DrawLatex(0.15,0.55,"N_{ch}^{S}: -2.4 < #eta < -1.2");
		if(i == 2) lat->DrawLatex(0.15,0.55,"N_{ch}^{M}: -1.0 < #eta < 1.0");

		if(i == 0) lat->DrawLatex(0.15,0.49,"#tilde{N}_{ch}^{N}: 1.2 < #eta < 2.4, J/#psi #rightarrow #mu^{+} #mu^{-} subtracted");
		if(i == 1) lat->DrawLatex(0.15,0.49,"#tilde{N}_{ch}^{S}: -2.4 < #eta < -1.2, J/#psi #rightarrow #mu^{+} #mu^{-} subtracted");



		c->SaveAs(Form("Plots/Together/png/%s.png",OutFileName[i].Data()));
		c->SaveAs(Form("Plots/Together/pdf/%s.pdf",OutFileName[i].Data()));

	}








}




