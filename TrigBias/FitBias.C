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



void FitBias(){
		gStyle->SetOptStat(0);


		TCanvas * c = new TCanvas("c","c",600,600);
		c->cd();



	const int NFiles = 3;

	TString FileName[NFiles] = {"Infiles/Run15pp_BBC_TrigEff_Multi_BBC_600-800kHz.root","Infiles/Run15pp_BBC_TrigEff_Multi_BBC_1000-1500kHz.root","Infiles/Run15pp_BBC_TrigEff_Multi_BBC_2000-2500kHz.root"};

	TString TitleName[NFiles] = {"BBC Rate: 600-800 kHz","BBC Rate: 1000-1500 kHz","BBC Rate: 2000-2500 kHz"};


	TH1D * h3FVTXS[NFiles];
	TH1D * h3FVTXN[NFiles];
	TH1D * h3SVX[NFiles];


	TF1 * FitNorth[NFiles];
	TF1 * FitSouth[NFiles];
	TF1 * FitSVX[NFiles];

	for(int i = 0; i < NFiles; i++){




		TFile * fin = new TFile(FileName[i].Data());
		fin->cd();


		h3FVTXS[i] = (TH1D *) fin->Get("h3FVTXS");
		h3FVTXN[i] = (TH1D *) fin->Get("h3FVTXN");
		h3SVX[i] = (TH1D *) fin->Get("h3SVX");

		h3FVTXS[i]->GetXaxis()->SetTitle("FVTX_South MB Event Multiplicity");
		h3FVTXS[i]->GetYaxis()->SetTitle("MB Trigger Bias");

		h3FVTXN[i]->GetXaxis()->SetTitle("FVTX_North MB Event Multiplicity");
		h3FVTXN[i]->GetYaxis()->SetTitle("MB Trigger Bias");


		h3SVX[i]->GetXaxis()->SetTitle("SVX MB Event Multiplicity");
		h3SVX[i]->GetYaxis()->SetTitle("MB Trigger Bias");


		FitNorth[i] = new TF1("FitNorth","[0] - [1] * TMath::Exp(-[2] *x)",1,12);
		FitSouth[i] = new TF1("FitSouth","[0] - [1] * TMath::Exp(-[2] *x)",1,12);
		FitSVX[i]= new TF1("FitSVX","[0] - [1] * TMath::Exp(-[2] *x)",1,12);
	

		h3FVTXS[i]->SetMarkerStyle(20);
		h3FVTXS[i]->SetMarkerSize(1);
		h3FVTXS[i]->SetMarkerColor(1);
		h3FVTXS[i]->SetLineColor(1);

		h3FVTXN[i]->SetMarkerStyle(20);
		h3FVTXN[i]->SetMarkerSize(1);
		h3FVTXN[i]->SetMarkerColor(1);
		h3FVTXN[i]->SetLineColor(1);

		h3SVX[i]->SetMarkerStyle(20);
		h3SVX[i]->SetMarkerSize(1);
		h3SVX[i]->SetMarkerColor(1);
		h3SVX[i]->SetLineColor(1);



		FitNorth[i]->SetParLimits(0,0.8,1.0);
		FitNorth[i]->SetParLimits(1,0,0.5);
		FitNorth[i]->SetParLimits(2,0,10);


		FitSouth[i]->SetParLimits(0,0.8,1.0);
		FitSouth[i]->SetParLimits(1,0,0.5);
		FitSouth[i]->SetParLimits(2,0,10);

		FitSVX[i]->SetParLimits(0,0.8,1.0);
		FitSVX[i]->SetParLimits(1,0,0.5);
		FitSVX[i]->SetParLimits(2,0,10);



		h3FVTXS[i]->Fit(FitSouth[i],"R");



		float p0S = FitSouth[i]->GetParameter(0);
		float p1S = FitSouth[i]->GetParameter(1);
		float p2S = FitSouth[i]->GetParameter(2);



		h3FVTXN[i]->Fit(FitNorth[i],"R");

		float p0N = FitNorth[i]->GetParameter(0);
		float p1N = FitNorth[i]->GetParameter(1);
		float p2N = FitNorth[i]->GetParameter(2);


		h3SVX[i]->Fit(FitSVX[i],"R");

		float p0SVX = FitSVX[i]->GetParameter(0);
		float p1SVX = FitSVX[i]->GetParameter(1);
		float p2SVX = FitSVX[i]->GetParameter(2);

		TString FuncNorthName = Form("%f - %f  * TMath::Exp(-%f  *x)",p0N,p1N,p2N);
		TString FuncSouthName = Form("%f - %f  * TMath::Exp(-%f  *x)",p0S,p1S,p2S);
		TString FuncSVXhName = Form("%f - %f  * TMath::Exp(-%f  *x)",p0SVX,p1SVX,p2SVX);

		cout << "FuncNorthName : " << FuncNorthName.Data() << endl;
		cout << "FuncSouthName : " << FuncSouthName.Data() << endl;
		cout << "FuncSVXhName : " << FuncSVXhName.Data() << endl;



		h3FVTXS[i]->Draw();
		c->SaveAs(Form("Plots/%d/h3FVTXS.png",i));


		h3FVTXN[i]->Draw();
		c->SaveAs(Form("Plots/%d/h3FVTXN.png",i));


		h3SVX[i]->Draw();
		c->SaveAs(Form("Plots/%d/h3SVX.png",i));

	
		h3FVTXN[i]->SetMarkerColor(kBlue);
		h3FVTXN[i]->SetLineColor(kBlue);
		
	
	
		h3FVTXS[i]->SetMarkerColor(kGreen);
		h3FVTXS[i]->SetLineColor(kGreen);
	
		h3SVX[i]->SetMarkerColor(6);
		h3SVX[i]->SetLineColor(6);

		FitNorth[i]->SetLineColor(kBlue);
		FitSouth[i]->SetLineColor(kGreen);
		FitSVX[i]->SetLineColor(6);

		h3FVTXS[i]->SetTitle(TitleName[i].Data());


		h3FVTXS[i]->Draw("ep");
		h3FVTXN[i]->Draw("epSAME");
		h3SVX[i]->Draw("epSAME");
		
		FitNorth[i]->Draw("SAME");
		FitSouth[i]->Draw("SAME");
		FitSVX[i]->Draw("SAME");
		


		TLegend *leg = new TLegend(0.42,0.25,0.69,0.45,NULL,"brNDC"); 
		leg->SetBorderSize(0);
		leg->SetTextSize(0.05);
		leg->SetTextFont(42);
		leg->SetFillStyle(0);

		leg->AddEntry(h3FVTXN[i],"FVTX North","pl");
		leg->AddEntry(h3FVTXS[i],"FVTX South","pl");
		leg->AddEntry(h3SVX[i],"SVX","pl");

		leg->Draw("SAME");

		
		c->SaveAs(Form("Plots/Summary/Plot_%d.png",i));

	}

	int Color[NFiles] = {2,3,6};
	TLegend *leg = new TLegend(0.18,0.25,0.69,0.45,NULL,"brNDC"); 

	for(int i = 0; i < NFiles; i++){
		h3FVTXN[i]->SetLineColor(Color[i]);
		h3FVTXN[i]->SetMarkerColor(Color[i]);		
		FitNorth[i]->SetLineColor(Color[i]);


		if(i == 0){
			h3FVTXN[i]->SetTitle("FVTX North");
			h3FVTXN[i]->Draw("ep");
		}
		if(i > 0){
	
			h3FVTXN[i]->Draw("epSAME");
		}
		FitNorth[i]->Draw("SAME");

		

		leg->SetBorderSize(0);
		leg->SetTextSize(0.05);
		leg->SetTextFont(42);
		leg->SetFillStyle(0);

		leg->AddEntry(h3FVTXN[i],TitleName[i].Data(),"pl");



	}
	    leg->Draw("SAME");

	c->SaveAs("Plots/Systems/North.png");



	for(int i = 0; i < NFiles; i++){
		h3FVTXS[i]->SetLineColor(Color[i]);
		h3FVTXS[i]->SetMarkerColor(Color[i]);		
		FitSouth[i]->SetLineColor(Color[i]);


		if(i == 0){
			h3FVTXS[i]->SetTitle("FVTX South");			
			h3FVTXS[i]->Draw("ep");
		}
		if(i > 0){
	
			h3FVTXS[i]->Draw("epSAME");
		}
		FitSouth[i]->Draw("SAME");


	}
	    leg->Draw("SAME");

	c->SaveAs("Plots/Systems/South.png");



	for(int i = 0; i < NFiles; i++){
		h3SVX[i]->SetLineColor(Color[i]);
		h3SVX[i]->SetMarkerColor(Color[i]);		
		FitSVX[i]->SetLineColor(Color[i]);


		if(i == 0){
			h3SVX[i]->SetTitle("SVX");						
			h3SVX[i]->Draw("ep");
		}
		if(i > 0){
	
			h3SVX[i]->Draw("epSAME");
		}
		FitSVX[i]->Draw("SAME");


	}
	    leg->Draw("SAME");

	c->SaveAs("Plots/Systems/SVX.png");

}

