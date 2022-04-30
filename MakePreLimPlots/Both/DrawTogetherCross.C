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


	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.0308);

	gStyle->SetOptStat(0);
	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();
	c->SetBottomMargin(0.12);



	const int NTypes = 2;
	const int NDir = 2;

	TString infile;
	TString infile2;
	TString infile3;

	TString DataName[NDir]={"1.2 < y_{J/#psi} < 2.2","-2.2 < y_{J/#psi} < -1.2"};
	TString Together = "J/#psi #rightarrow #mu^{+} #mu^{-} subtracted";
	TLegend* leg[NTypes];

	float Min[NTypes] = {0,0};
	float Max[NTypes] = {18,18};
	TString OutFileName[NDir]={"FVTXN","FVTXS"};


	TF1 * func = new TF1("Func","x",0,19);
	func->SetLineStyle(2);

	for(int i = 0; i < NTypes; i++){


		leg[i] = new TLegend(0.12,0.66,0.45,0.80,NULL,"brNDC");
		leg[i]->SetBorderSize(0);
		leg[i]->SetTextSize(0.033);
		leg[i]->SetTextFont(42);
		leg[i]->SetFillStyle(0);
		leg[i]->SetLineWidth(3);






		infile = Form("/sphenix/user/zshi/PHENIX/JPsiAnaCodes/MakePreLimPlots/UnCorr/OutFiles/JPsiR_%d_%d.root",i,i);
		infile2 = Form("/sphenix/user/zshi/PHENIX/JPsiAnaCodes/MakePreLimPlots/Corr/OutFiles/JPsiR_%d_%d.root",i,i);
		infile3 = Form("/sphenix/user/zshi/PHENIX/JPsiAnaCodes/MakePreLimPlots/UnCorr/OutFiles/JPsiR_%d_%d.root",1-i,i);

		TFile * fin = new TFile(infile.Data());
		TFile * fin2 = new TFile(infile2.Data());
		TFile * fin3 = new TFile(infile3.Data());

		TH1D * RatioScaledCorr = (TH1D * ) fin->Get("RatioScaledCorr");



		TH1D * RatioScaledCorr2 = (TH1D * ) fin2->Get("RatioScaledCorr");


		TH1D * RatioScaledCorr3 = (TH1D * ) fin3->Get("RatioScaledCorr");



		TGraphAsymmErrors *Syst = (TGraphAsymmErrors *) fin->Get("Graph");
		TGraphAsymmErrors *Syst2 = (TGraphAsymmErrors *) fin2->Get("Graph");
		TGraphAsymmErrors *Syst3 = (TGraphAsymmErrors *) fin3->Get("Graph");

		RatioScaledCorr->SetMarkerColor(i+3);
		RatioScaledCorr->SetLineColor(i+3);


		if(i == 0){
		
			Syst->SetFillColorAlpha(kGreen-9,0.5);
			Syst->SetLineColor(kGreen-9);


			Syst3->SetFillColorAlpha(kBlue-9,0.5);
			Syst3->SetLineColor(kBlue-9);



		}


		if(i == 1){
		
			Syst->SetFillColorAlpha(kBlue-9,0.5);
			Syst->SetLineColor(kBlue-9);


			Syst3->SetFillColorAlpha(kGreen-9,0.5);
			Syst3->SetLineColor(kGreen-9);

		}


		RatioScaledCorr2->SetMarkerColor(kRed);
		RatioScaledCorr2->SetLineColor(kRed);



		Syst2->SetFillColorAlpha(kRed-9,0.5);
		Syst2->SetLineColor(kRed-9);



		RatioScaledCorr3->SetMarkerColor(4-i);
		RatioScaledCorr3->SetLineColor(4-i);



		RatioScaledCorr->SetMinimum(Min[i]);
		RatioScaledCorr->SetMaximum(Max[i]);				
		RatioScaledCorr->GetXaxis()->SetTitleOffset(0.90);

		RatioScaledCorr->Draw("ep");
		
		Syst->Draw("5SAME");



		leg[i]->AddEntry(RatioScaledCorr,DataName[i].Data(),"LP");
		leg[i]->Draw("SAME");

		lat->DrawLatex(0.125,0.55,"pp #sqrt{s} = 200 GeV");
		if(i == 0) lat->DrawLatex(0.125,0.49,"N^{N}_{ch}: 1.2 < #eta < 2.4");
		if(i == 1) lat->DrawLatex(0.125,0.49,"N^{S}_{ch}: -2.4 < #eta < -1.2");



		c->SaveAs(Form("Plots/CrossComp/png/%s_Part1.png",OutFileName[i].Data()));
		c->SaveAs(Form("Plots/CrossComp/pdf/%s_Part1.pdf",OutFileName[i].Data()));

		if(i == 0) RatioScaledCorr->GetXaxis()->SetTitle("N_{ch}^{N}/<N_{ch}^{N}>, #tilde{N}_{ch}^{N}/<#tilde{N}_{ch}^{N}>");
		if(i == 1) RatioScaledCorr->GetXaxis()->SetTitle("N_{ch}^{S}/<N_{ch}^{S}>, #tilde{N}_{ch}^{S}/<#tilde{N}_{ch}^{S}>");


		RatioScaledCorr->Draw("ep");
		Syst->Draw("5SAME");


		RatioScaledCorr2->Draw("epSAME");
		Syst2->Draw("5SAME");

		leg[i]->AddEntry(RatioScaledCorr2,DataName[i].Data(),"LP");
//		leg[i]->AddEntry("",Together.Data(),"");
	
		leg[i]->Draw("SAME");


		lat->DrawLatex(0.125,0.60,"pp #sqrt{s} = 200 GeV");
		if(i == 0) lat->DrawLatex(0.125,0.54,"N^{N}_{ch}: 1.2 < #eta < 2.4");
		if(i == 1) lat->DrawLatex(0.125,0.54,"N^{S}_{ch}: -2.4 < #eta < -1.2");

		if(i == 0) lat->DrawLatex(0.125,0.49,"#tilde{N}^{N}_{ch}: 1.2 < #eta < 2.4");
		if(i == 1) lat->DrawLatex(0.125,0.49,"#tilde{N}^{S}_{ch}: -2.4 < #eta < -1.2");

		lat->DrawLatex(0.125,0.44,"J/#psi #rightarrow #mu^{+} #mu^{-} subtracted");





		c->SaveAs(Form("Plots/CrossComp/png/%s_Part2.png",OutFileName[i].Data()));
		c->SaveAs(Form("Plots/CrossComp/pdf/%s_Part2.pdf",OutFileName[i].Data()));


			
		if(i == 0) RatioScaledCorr->GetXaxis()->SetTitle("N_{ch}^{N}/<N_{ch}^{N}>, #tilde{N}_{ch}^{N}/<#tilde{N}_{ch}^{N}>");
		if(i == 1) RatioScaledCorr->GetXaxis()->SetTitle("N_{ch}^{S}/<N_{ch}^{S}>, #tilde{N}_{ch}^{S}/<#tilde{N}_{ch}^{S}>");




		RatioScaledCorr->Draw("ep");
		Syst->Draw("5SAME");



		RatioScaledCorr2->Draw("epSAME");
		Syst2->Draw("5SAME");



		RatioScaledCorr3->Draw("epSAME");
		Syst3->Draw("5SAME");


		leg[i]->AddEntry(RatioScaledCorr3,DataName[1-i].Data(),"LP");
		leg[i]->Draw("SAME");


		lat->DrawLatex(0.125,0.60,"pp #sqrt{s} = 200 GeV");
		if(i == 0) lat->DrawLatex(0.125,0.54,"N^{N}_{ch}: 1.2 < #eta < 2.4");
		if(i == 1) lat->DrawLatex(0.125,0.54,"N^{S}_{ch}: -2.4 < #eta < -1.2");

		if(i == 0) lat->DrawLatex(0.125,0.49,"#tilde{N}^{N}_{ch}: 1.2 < #eta < 2.4");
		if(i == 1) lat->DrawLatex(0.125,0.49,"#tilde{N}^{S}_{ch}: -2.4 < #eta < -1.2");

		lat->DrawLatex(0.125,0.44,"J/#psi #rightarrow #mu^{+} #mu^{-} subtracted");


	
		c->SaveAs(Form("Plots/CrossComp/png/%s_Part3.png",OutFileName[i].Data()));
		c->SaveAs(Form("Plots/CrossComp/pdf/%s_Part3.pdf",OutFileName[i].Data()));

	}








}




