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



void FitBiasPreLim(){



	gStyle->SetOptStat(0);
	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();



	TFile * fin = new TFile("Infiles/Run15pp_BBC_TrigEff_Multi_Run421xxx.root");
	fin->cd();
	
	TH1D * h3FVTXS = (TH1D *) fin->Get("h3FVTXS");
	TH1D * h3FVTXN = (TH1D *) fin->Get("h3FVTXN");
	TH1D * h3SVX = (TH1D *) fin->Get("h3SVX");
	
	h3FVTXS->GetXaxis()->SetTitle("FVTX_South MB Event Multiplicity");
	h3FVTXS->GetYaxis()->SetTitle("MB Trigger Bias");

	h3FVTXN->GetXaxis()->SetTitle("FVTX_North MB Event Multiplicity");
	h3FVTXN->GetYaxis()->SetTitle("MB Trigger Bias");


	h3SVX->GetXaxis()->SetTitle("SVX MB Event Multiplicity");
	h3SVX->GetYaxis()->SetTitle("MB Trigger Bias");


	TF1 * FitNorth = new TF1("FitNorth","[0] - [1] * TMath::Exp(-[2] *x)",1,12);
	TF1 * FitSouth = new TF1("FitSouth","[0] - [1] * TMath::Exp(-[2] *x)",1,12);
	TF1 * FitSVX = new TF1("FitSVX","[0] - [1] * TMath::Exp(-[2] *x)",1,12);


	h3FVTXS->SetMarkerStyle(20);
	h3FVTXS->SetMarkerSize(1);
	h3FVTXS->SetMarkerColor(1);
	h3FVTXS->SetLineColor(1);

	h3FVTXN->SetMarkerStyle(20);
	h3FVTXN->SetMarkerSize(1);
	h3FVTXN->SetMarkerColor(1);
	h3FVTXN->SetLineColor(1);

	h3SVX->SetMarkerStyle(20);
	h3SVX->SetMarkerSize(1);
	h3SVX->SetMarkerColor(1);
	h3SVX->SetLineColor(1);



	FitNorth->SetParLimits(0,0.8,1.0);
	FitNorth->SetParLimits(1,0,0.5);
	FitNorth->SetParLimits(2,0,10);


	FitSouth->SetParLimits(0,0.8,1.0);
	FitSouth->SetParLimits(1,0,0.5);
	FitSouth->SetParLimits(2,0,10);

	FitSVX->SetParLimits(0,0.8,1.0);
	FitSVX->SetParLimits(1,0,0.5);
	FitSVX->SetParLimits(2,0,10);



	h3FVTXS->Fit(FitSouth,"R");



	float p0S = FitSouth->GetParameter(0);
	float p1S = FitSouth->GetParameter(1);
	float p2S = FitSouth->GetParameter(2);



	h3FVTXN->Fit(FitNorth,"R");
	
	float p0N = FitNorth->GetParameter(0);
	float p1N = FitNorth->GetParameter(1);
	float p2N = FitNorth->GetParameter(2);

	
	h3SVX->Fit(FitSVX,"R");
	
	float p0SVX = FitSVX->GetParameter(0);
	float p1SVX = FitSVX->GetParameter(1);
	float p2SVX = FitSVX->GetParameter(2);

	TString FuncNorthName = Form("%f - %f  * TMath::Exp(-%f  *x)",p0N,p1N,p2N);
	TString FuncSouthName = Form("%f - %f  * TMath::Exp(-%f  *x)",p0S,p1S,p2S);
	TString FuncSVXhName = Form("%f - %f  * TMath::Exp(-%f  *x)",p0SVX,p1SVX,p2SVX);

	cout << "FuncNorthName : " << FuncNorthName.Data() << endl;
	cout << "FuncSouthName : " << FuncSouthName.Data() << endl;
	cout << "FuncSVXhName : " << FuncSVXhName.Data() << endl;
	

	
	h3FVTXS->Draw();
	c->SaveAs("Plots/PreLim/h3FVTXS.png");

	
	h3FVTXN->Draw();
	c->SaveAs("Plots/PreLim/h3FVTXN.png");

	
	h3SVX->Draw();
	c->SaveAs("Plots/PreLim/h3SVX.png");



}

