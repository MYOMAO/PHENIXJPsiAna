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


	const int NFiles = 3;

	TString FileName[NFiles] = {"infiles/Run15pp_BBC_TrigEff_Multi_BBC_600-800kHz.root","infiles/Run15pp_BBC_TrigEff_Multi_BBC_1000-1500kHz.root","infiles/Run15pp_BBC_TrigEff_Multi_BBC_2000-2500kHz.root"};

	TString TitleName[NFiles] = {"BBC Rate: 600-800 kHz","BBC Rate: 1000-1500 kHz","BBC Rate: 2000-2500 kHz"};


	TH1D * h3FVTXS[NFiles];
	TH1D * h3FVTXN[NFiles];
	TH1D * h3SVX[NFiles];


	TF1 * FitNorth[NFiles];
	TF1 * FitSouth[NFiles];
	TF1 * FitSVX[NFiles];

	for(int i = 0; i < NFiles; i++){


		gStyle->SetOptStat(0);
		TCanvas * c = new TCanvas("c","c",600,600);
		c->cd();



		TFile * fin = new TFile(FileName[i].Data());
		fin->cd();


		h3FVTXS[i] = (TH1D *) fin->Get("h3FVTXS");
		h3FVTXN[i] = (TH1D *) fin->Get("h3FVTXN");
		h3SVX[i] = (TH1D *) fin->Get("h3SVX");

		h3FVTXS[i]->GetXaxis()->SetTitle("FVTX_South MB Event Multiplicity");
		h3FVTXS[i]->GetYaxis()->SetTitle("MB Trigger Bias");
		h3FVTXS[i]->GetYaxis()->SetTitleOffset(1.4);
		h3FVTXS[i]->GetXaxis()->CenterTitle();
		h3FVTXS[i]->GetYaxis()->CenterTitle();

		h3FVTXN[i]->GetXaxis()->SetTitle("FVTX_North MB Event Multiplicity");
		h3FVTXN[i]->GetYaxis()->SetTitle("MB Trigger Bias");
		h3FVTXN[i]->GetYaxis()->SetTitleOffset(1.4);
		h3FVTXN[i]->GetXaxis()->CenterTitle();
		h3FVTXN[i]->GetYaxis()->CenterTitle();


		h3SVX[i]->GetXaxis()->SetTitle("SVX MB Event Multiplicity");
		h3SVX[i]->GetYaxis()->SetTitle("MB Trigger Bias");
		h3SVX[i]->GetYaxis()->SetTitleOffset(1.4);
		h3SVX[i]->GetXaxis()->CenterTitle();
		h3SVX[i]->GetYaxis()->CenterTitle();


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
		c->SaveAs(Form("TrigBias/%d/h3FVTXS.png",i));


		h3FVTXN[i]->Draw();
		c->SaveAs(Form("TrigBias/%d/h3FVTXN.png",i));


		h3SVX[i]->Draw();
		c->SaveAs(Form("TrigBias/%d/h3SVX.png",i));


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


		c->SaveAs(Form("TrigBias/Summary/Plot_%d.png",i));

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

	c->SaveAs("TrigBias/Systems/North.png");



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

	c->SaveAs("TrigBias/Systems/South.png");



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

	c->SaveAs("TrigBias/Systems/SVX.png");

	//Evaluate Syst

	const int NBins = 10;
	double MultBin[NBins + 1] = {0,1,2,3,4,5,6,8,10,12,19};

	TString XName[NFiles] = {"FVTXN MB Event Multiplicity","FVTXS MB Event Multiplicity","SVX MB Event Multiplicity"};
	TString XNameHis[NFiles] = {"FVTXN.png","FVTXS.png","SVX.png"};
	TString Title[NFiles] = {"MB Trigger Bias Systematics: FVTXN","MB Trigger Bias Systematics: FVTXS","MB Trigger Bias Systematics: SVX"};

	TH1D * TrigBiasSystUp[NFiles];	
	TH1D * TrigBiasSystDown[NFiles];	


	float XCenter;
	float FVTXNVar[NFiles];
	float FVTXSVar[NFiles];
	float SVXVar[NFiles];

	float FVTXNSystUp;
	float FVTXNSystDown;
	float FVTXSSystUp;
	float FVTXSSystDown;
	float SVXSystUp;
	float SVXSystDown;


	for(int i = 0; i < NFiles; i++){

		TrigBiasSystUp[i] = new TH1D(Form("TrigBiasSystUp_%d",i),"",NBins,MultBin);
		TrigBiasSystUp[i]->GetXaxis()->SetTitle(XName[i].Data());
		TrigBiasSystUp[i]->GetYaxis()->SetTitle("MB Trigger Systematics (%)");
		TrigBiasSystUp[i]->GetYaxis()->SetTitleOffset(1.4);
		TrigBiasSystUp[i]->GetXaxis()->CenterTitle();
		TrigBiasSystUp[i]->GetYaxis()->CenterTitle();
		TrigBiasSystUp[i]->SetMarkerStyle(20);
		TrigBiasSystUp[i]->SetMarkerSize(1);
		TrigBiasSystUp[i]->SetMarkerColor(kRed);
		TrigBiasSystUp[i]->SetLineColor(kRed);
		TrigBiasSystUp[i]->SetTitle(Title[i].Data());

		TrigBiasSystDown[i] = new TH1D(Form("TrigBiasSystDown_%d",i),"",NBins,MultBin);
		TrigBiasSystDown[i]->GetXaxis()->SetTitle(XName[i].Data());
		TrigBiasSystDown[i]->GetYaxis()->SetTitle("MB Trigger Systematics (%)");
		TrigBiasSystDown[i]->GetYaxis()->SetTitleOffset(1.4);
		TrigBiasSystDown[i]->GetXaxis()->CenterTitle();
		TrigBiasSystDown[i]->GetYaxis()->CenterTitle();
		TrigBiasSystDown[i]->SetMarkerStyle(20);
		TrigBiasSystDown[i]->SetMarkerSize(1);
		TrigBiasSystDown[i]->SetMarkerColor(kBlue);
		TrigBiasSystDown[i]->SetLineColor(kBlue);


	}





	for(int q = 0; q < NBins; q++){

		XCenter = (MultBin[q] + MultBin[q+1])/2;

		for(int i = 0; i < NFiles; i++){

			FVTXNVar[i] = FitNorth[i]->Eval(XCenter);
			FVTXSVar[i] = FitSouth[i]->Eval(XCenter);
			SVXVar[i] = FitSVX[i]->Eval(XCenter);

		}

		FVTXNSystUp = TMath::Abs(FVTXNVar[2] - FVTXNVar[1])/FVTXNVar[1] * 100;
		FVTXSSystUp = TMath::Abs(FVTXSVar[2] - FVTXSVar[1])/FVTXSVar[1] * 100;
		SVXSystUp = TMath::Abs(SVXVar[2] - SVXVar[1])/SVXVar[1] * 100;

		FVTXNSystDown = -TMath::Abs(FVTXNVar[0] - FVTXNVar[1])/FVTXNVar[1] * 100;
		FVTXSSystDown = -TMath::Abs(FVTXSVar[0] - FVTXSVar[1])/FVTXSVar[1] * 100;
		SVXSystDown = -TMath::Abs(SVXVar[0] - SVXVar[1])/SVXVar[1] * 100;

		//cout << "FVTXNVar[2] = " << FVTXNVar[2] <<  "      FVTXNVar[1] = " << FVTXNVar[1]  << endl;


		TrigBiasSystUp[0]->SetBinContent(q+1,FVTXNSystUp);
		TrigBiasSystUp[1]->SetBinContent(q+1,FVTXNSystUp);
		TrigBiasSystUp[2]->SetBinContent(q+1,SVXSystUp);

		TrigBiasSystDown[0]->SetBinContent(q+1,FVTXNSystDown);
		TrigBiasSystDown[1]->SetBinContent(q+1,FVTXNSystDown);
		TrigBiasSystDown[2]->SetBinContent(q+1,SVXSystDown);

		TrigBiasSystUp[0]->SetBinError(q+1,0.0001);
		TrigBiasSystUp[1]->SetBinError(q+1,0.0001);
		TrigBiasSystUp[2]->SetBinError(q+1,0.0001);

		TrigBiasSystDown[0]->SetBinError(q+1,0.0001);
		TrigBiasSystDown[1]->SetBinError(q+1,0.0001);
		TrigBiasSystDown[2]->SetBinError(q+1,0.0001);


	}



	for(int i = 0; i < NFiles; i++){


		TrigBiasSystUp[i]->SetMaximum(10);
		TrigBiasSystUp[i]->SetMinimum(-10);
		
		TrigBiasSystUp[i]->Draw("ep");
		TrigBiasSystDown[i]->Draw("epSAME");

		TLegend *leg = new TLegend(0.42,0.15,0.69,0.30,NULL,"brNDC"); 
		leg->SetBorderSize(0);
		leg->SetTextSize(0.05);
		leg->SetTextFont(42);
		leg->SetFillStyle(0);

		leg->AddEntry(TrigBiasSystUp[i],"Upper Bound","LP");			
		leg->AddEntry(TrigBiasSystDown[i],"Lower Bound","LP");			

		leg->Draw("SAME");

		c->SaveAs(Form("SystVar/%s",XNameHis[i].Data()));
	}



}

