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



	const int NTypes = 3;
	const int NDir = 2;

	TString infile;

	TString OutFileName[NTypes]={"FVTXN","FVTXS","SVX"};
	TString DataName[NDir]={"1.2 < y_{J/#psi} < 2.2","-2.2 < y_{J/#psi} < -1.2"};

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
				RatioScaledCorr->SetTitle("");
				RatioScaledCorr->Draw("ep");
			}
			
			if(j>0) RatioScaledCorr->Draw("epSAME");


			leg[i]->AddEntry(RatioScaledCorr,DataName[j].Data(),"LP");
			Syst->Draw("5SAME");

		}



		leg[i]->Draw("SAME");

	
		lat->DrawLatex(0.15,0.61,"pp #sqrt{s} = 200 GeV");
		if(i == 0) lat->DrawLatex(0.15,0.55,"N_{ch}^{N}: 1.2 < #eta < 2.4");
		if(i == 1) lat->DrawLatex(0.15,0.55,"N_{ch}^{S}: -2.4 < #eta < -1.2");
		if(i == 2) lat->DrawLatex(0.15,0.55,"N_{ch}^{M}: -1.0 < #eta < 1.0");

		c->SaveAs(Form("Plots/Together/png/%s.png",OutFileName[i].Data()));
		c->SaveAs(Form("Plots/Together/pdf/%s.pdf",OutFileName[i].Data()));

	}








}




