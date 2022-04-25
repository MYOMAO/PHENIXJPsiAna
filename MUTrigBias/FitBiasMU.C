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



void FitBiasMU(){
	





	gStyle->SetOptStat(0);
	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();



	TFile * fin = new TFile("Run15pp_BBC_TrigEff_Multi_BBC_all_ERT.root");
	fin->cd();
	
	TH1D * h1FVTXS = (TH1D *) fin->Get("h1FVTXS");
	TH1D * h1FVTXN = (TH1D *) fin->Get("h1FVTXN");
	TH1D * h1SVX = (TH1D *) fin->Get("h1SVX");

	TH1D * h2FVTXS = (TH1D *) fin->Get("h2FVTXS");
	TH1D * h2FVTXN = (TH1D *) fin->Get("h2FVTXN");
	TH1D * h2SVX = (TH1D *) fin->Get("h2SVX");


//	const int NBins = 10;
//	double MultBin[NBins + 1] = {0,1,2,3,4,5,6,8,10,12,19};


	const int NBins = 6;
	double MultBin[NBins + 1] = {0,1,2,4,5,8,19};


	TH1D * h1FVTXSNew = new TH1D("h1FVTXSNew","",NBins, MultBin);
	TH1D * h1FVTXNNew = new TH1D("h1FVTXNNew","",NBins, MultBin);
	TH1D * h1SVXNew = new TH1D("h1SVXNew","",NBins, MultBin);

	TH1D * h2FVTXSNew = new TH1D("h2FVTXSNew","",NBins, MultBin);
	TH1D * h2FVTXNNew = new TH1D("h2FVTXNNew","",NBins, MultBin);
	TH1D * h2SVXNew = new TH1D("h2SVXNew","",NBins, MultBin);





	const int NBinsX = 15; 

	float h1FVTXSValue[NBinsX];
	float h2FVTXSValue[NBinsX];

	float h1FVTXNValue[NBinsX];
	float h2FVTXNValue[NBinsX];

	float h1SVXValue[NBinsX];
	float h2SVXValue[NBinsX];
	


	float XValue;
	int XBin;

	float SystFVTXN;
	float SystFVTXS;
	float SystSVX;

	for(int i = 0; i < NBinsX; i++){

		h1FVTXSValue[i] = h1FVTXS->GetBinContent(i+1);
		h2FVTXSValue[i] = h2FVTXS->GetBinContent(i+1);
		
		h1FVTXNValue[i] = h1FVTXN->GetBinContent(i+1);
		h2FVTXNValue[i] = h2FVTXN->GetBinContent(i+1);

		h1SVXValue[i] = h1SVX->GetBinContent(i+1);
		h2SVXValue[i] = h2SVX->GetBinContent(i+1);
		


		XValue = h1FVTXS->GetBinCenter(i+1);
	

		h1FVTXSNew->Fill(XValue,h1FVTXSValue[i]);
		h1FVTXNNew->Fill(XValue,h1FVTXNValue[i]);
		h1SVXNew->Fill(XValue,h1SVXValue[i]);

		h2FVTXSNew->Fill(XValue,h2FVTXSValue[i]);
		h2FVTXNNew->Fill(XValue,h2FVTXNValue[i]);
		h2SVXNew->Fill(XValue,h2SVXValue[i]);


		if(i == 1){
		
			SystFVTXN =h1FVTXN->GetBinError(i+1)/h1FVTXNValue[i] * 100;
			SystFVTXS =h1FVTXS->GetBinError(i+1)/h1FVTXSValue[i] * 100;
			SystSVX =h1SVX->GetBinError(i+1)/h1SVXValue[i] * 100;

		}

	}



	
	
	h1FVTXSNew->Sumw2();
	h2FVTXSNew->Sumw2();

	h1FVTXNNew->Sumw2();
	h2FVTXNNew->Sumw2();
	
	h1SVXNew->Sumw2();
	h2SVXNew->Sumw2();
	

	for(int i = 0; i < NBins; i++){

		cout << "h1FVTXSNew->GetBinContent(i+1) = " << h1FVTXSNew->GetBinContent(i+1) << "    h2FVTXSNew->GetBinContent(i+1) =   " << h2FVTXSNew->GetBinContent(i+1)  << endl;
		cout << "h1FVTXNNew->GetBinContent(i+1) = " << h1FVTXNNew->GetBinContent(i+1) << "    h2FVTXNNew->GetBinContent(i+1) =   " << h2FVTXNNew->GetBinContent(i+1)  << endl;
		cout << "h1SVXNew->GetBinContent(i+1) = " << h1SVXNew->GetBinContent(i+1) << "    h2SVXNew->GetBinContent(i+1) =   " << h2SVXNew->GetBinContent(i+1)  << endl;

		h1FVTXSNew->SetBinError(i+1,sqrt(h1FVTXSNew->GetBinContent(i+1)));
		h2FVTXSNew->SetBinError(i+1,sqrt(h2FVTXSNew->GetBinContent(i+1)));
		
		h1FVTXNNew->SetBinError(i+1,sqrt(h1FVTXNNew->GetBinContent(i+1)));
		h2FVTXNNew->SetBinError(i+1,sqrt(h2FVTXNNew->GetBinContent(i+1)));

		h1SVXNew->SetBinError(i+1,sqrt(h1SVXNew->GetBinContent(i+1)));
		h2SVXNew->SetBinError(i+1,sqrt(h2SVXNew->GetBinContent(i+1)));


	}


	TH1D * FVTXSRatioNew = (TH1D *) h2FVTXSNew->Clone("FVTXSRatioNew");
	TH1D * FVTXNRatioNew = (TH1D *) h2FVTXNNew->Clone("FVTXNRatioNew");
	TH1D * SVXRatioNew = (TH1D *) h2SVXNew->Clone("SVXRatioNew");


	FVTXSRatioNew->Sumw2();
	FVTXNRatioNew->Sumw2();
	SVXRatioNew->Sumw2();

	FVTXSRatioNew->Divide(h1FVTXSNew);
	FVTXNRatioNew->Divide(h1FVTXNNew);
	SVXRatioNew->Divide(h1SVXNew);


	FVTXSRatioNew->GetXaxis()->SetTitle("FVTX South Event Multiplicity");
	FVTXSRatioNew->GetYaxis()->SetTitle("J/#psi Trigger Efficiency #epsilon_{trig}^{J/#psi}");
	FVTXSRatioNew->GetYaxis()->SetTitleOffset(1.2);


	FVTXSRatioNew->GetXaxis()->CenterTitle();
	FVTXSRatioNew->GetYaxis()->CenterTitle();
	FVTXSRatioNew->SetTitle("FVXT South ERT 4 #times 4 Trigger Efficiency");

	FVTXNRatioNew->GetXaxis()->SetTitle("FVTX North Event Multiplicity");
	FVTXNRatioNew->GetYaxis()->SetTitle("J/#psi Trigger Efficiency #epsilon_{trig}^{J/#psi}");
	FVTXNRatioNew->GetYaxis()->SetTitleOffset(1.2);


	
	FVTXNRatioNew->GetXaxis()->CenterTitle();
	FVTXNRatioNew->GetYaxis()->CenterTitle();
	FVTXNRatioNew->SetTitle("FVXT North ERT 4 #times 4 Trigger Efficiency");


	SVXRatioNew->GetXaxis()->SetTitle("SVX Event Multiplicity");
	SVXRatioNew->GetYaxis()->SetTitle("J/#psi Trigger Efficiency #epsilon_{trig}^{J/#psi}");
	SVXRatioNew->GetYaxis()->SetTitleOffset(1.2);

	
	SVXRatioNew->GetXaxis()->CenterTitle();
	SVXRatioNew->GetYaxis()->CenterTitle();
	SVXRatioNew->SetTitle("SVX ERT 4 #times 4 Trigger Efficiency");

	TF1 * FitNorth = new TF1("FitNorth","[0] - [1] * TMath::Exp(-[2] *x)",1,16);
	TF1 * FitSouth = new TF1("FitSouth","[0] - [1] * TMath::Exp(-[2] *x)",1,16);
	TF1 * FitSVX = new TF1("FitSVX","[0] - [1] * TMath::Exp(-[2] *x)",1,16);


	FVTXSRatioNew->SetMarkerStyle(20);
	FVTXSRatioNew->SetMarkerSize(1);
	FVTXSRatioNew->SetMarkerColor(1);
	FVTXSRatioNew->SetLineColor(1);

	FVTXNRatioNew->SetMarkerStyle(20);
	FVTXNRatioNew->SetMarkerSize(1);
	FVTXNRatioNew->SetMarkerColor(1);
	FVTXNRatioNew->SetLineColor(1);

	SVXRatioNew->SetMarkerStyle(20);
	SVXRatioNew->SetMarkerSize(1);
	SVXRatioNew->SetMarkerColor(1);
	SVXRatioNew->SetLineColor(1);



	FitNorth->SetParLimits(0,0.8,1.1);
	FitNorth->SetParLimits(1,0,0.5);
	FitNorth->SetParLimits(2,0,10);


	FitSouth->SetParLimits(0,0.8,1.0);
	FitSouth->SetParLimits(1,0,0.5);
	FitSouth->SetParLimits(2,0,10);

	FitSVX->SetParLimits(0,0.8,1.0);
	FitSVX->SetParLimits(1,0,0.5);
	FitSVX->SetParLimits(2,0,10);



	FVTXSRatioNew->Fit(FitSouth,"R");



	float p0S = FitSouth->GetParameter(0);
	float p1S = FitSouth->GetParameter(1);
	float p2S = FitSouth->GetParameter(2);


	float p0SErr = FitSouth->GetParError(0);
	float p1SErr = FitSouth->GetParError(1);
	float p2SErr = FitSouth->GetParError(2);



	FVTXNRatioNew->Fit(FitNorth,"R");
	
	float p0N = FitNorth->GetParameter(0);
	float p1N = FitNorth->GetParameter(1);
	float p2N = FitNorth->GetParameter(2);


	float p0NErr = FitNorth->GetParError(0);
	float p1NErr = FitNorth->GetParError(1);
	float p2NErr = FitNorth->GetParError(2);


	
	SVXRatioNew->Fit(FitSVX,"R");
	
	float p0SVX = FitSVX->GetParameter(0);
	float p1SVX = FitSVX->GetParameter(1);
	float p2SVX = FitSVX->GetParameter(2);

	float p0SVXErr = FitSVX->GetParError(0);
	float p1SVXErr = FitSVX->GetParError(1);
	float p2SVXErr = FitSVX->GetParError(2);



	TString FuncNorthName = Form("%.2f - %.2f  * exp(-%.2f  *x)",p0N,p1N,p2N);
	TString FuncSouthName = Form("%.2f - %.2f  * exp(-%.2f  *x)",p0S,p1S,p2S);
	TString FuncSVXName = Form("%.2f - %.2f  * :exp(-%.2f  *x)",p0SVX,p1SVX,p2SVX);

	/*

	TString FuncNorthName = Form("%f - %f  * TMath::Exp(-%f  *x)",p0N,p1N,p2N);
	TString FuncSouthName = Form("%f - %f  * TMath::Exp(-%f  *x)",p0S,p1S,p2S);
	TString FuncSVXName = Form("%f - %f  * TMath::Exp(-%f  *x)",p0SVX,p1SVX,p2SVX);

   */

	cout << "FuncNorthName : " << FuncNorthName.Data() << endl;
	cout << "FuncSouthName : " << FuncSouthName.Data() << endl;
	cout << "FuncSVXName : " << FuncSVXName.Data() << endl;
	
	
	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.040);


	
	FVTXSRatioNew->Draw();

	lat->DrawLatex(0.20,0.30,Form("y = %s",FuncSouthName.Data()));
	lat->DrawLatex(0.20,0.25,Form("Global Systematics: %.1f %%",SystFVTXS));
	
	c->SaveAs("Plots/FVTXSRatioNew.png");


	FVTXNRatioNew->Draw();
	lat->DrawLatex(0.20,0.30,Form("y = %s",FuncNorthName.Data()));	
	lat->DrawLatex(0.20,0.25,Form("Global Systematics: %.1f %%",SystFVTXN));


	c->SaveAs("Plots/FVTXNRatioNew.png");



	SVXRatioNew->Draw();
	
	lat->DrawLatex(0.20,0.30,Form("y = %s",FuncSVXName.Data()));
	lat->DrawLatex(0.20,0.25,Form("Global Systematics: %.1f %%",SystSVX));

	c->SaveAs("Plots/SVXRatioNew.png");

	

	cout << "SystFVTXS = " << SystFVTXS << endl;
	cout << "SystFVTXN = " << SystFVTXN << endl;
	cout << "SystSVX = " << SystSVX << endl;


	const int NBinsMult = 10;
	double MultBin2[NBinsMult + 1] = {0,1,2,3,4,5,6,8,10,12,19};

	float Nominal;
	float BinCenter;
	float VarUp;
	float VarDown;

	float SystUp;
	float SystDown;

/*
	TString FuncNorthUp = Form("%f - %f  * TMath::Exp(-%f  *x)",p0N + p0NErr,p1N - p1NErr,p2N);
	TString FuncSouthUp = Form("%f - %f  * TMath::Exp(-%f  *x)",p0S + p0SErr,p1S - p1SErr,p2S);
	TString FuncSVXUp = Form("%f - %f  * TMath::Exp(-%f  *x)",p0SVX + p0SVXErr,p1SVX - p1SVXErr,p2SVX);

	TString FuncNorthDown = Form("%f - %f  * TMath::Exp(-%f  *x)",p0N - p0NErr,p1N+  p1NErr,p2N);
	TString FuncSouthDown = Form("%f - %f  * TMath::Exp(-%f  *x)",p0S-p0SErr,p1S + p1SErr,p2S);
	TString FuncSVXDown = Form("%f - %f  * TMath::Exp(-%f  *x)",p0SVX-p0SVXErr,p1SVX + p1SVXErr,p2SVX);

	TString FuncNorthUp = Form("%f - %f  * TMath::Exp(-%f  *x)",p0N,p1N - p1NErr,p2N);
	TString FuncSouthUp = Form("%f - %f  * TMath::Exp(-%f  *x)",p0S,p1S - p1SErr,p2S);
	TString FuncSVXUp = Form("%f - %f  * TMath::Exp(-%f  *x)",p0SVX,p1SVX - p1SVXErr,p2SVX);

	TString FuncNorthDown = Form("%f - %f  * TMath::Exp(-%f  *x)",p0N,p1N+  p1NErr,p2N);
	TString FuncSouthDown = Form("%f - %f  * TMath::Exp(-%f  *x)",p0S,p1S + p1SErr,p2S);
	TString FuncSVXDown = Form("%f - %f  * TMath::Exp(-%f  *x)",p0SVX,p1SVX + p1SVXErr,p2SVX);


*/

	TString FuncNorthUp = Form("%f - %f  * TMath::Exp(-%f  *x)",p0N,p1N - p1NErr * 0.5,p2N);
	TString FuncSouthUp = Form("%f - %f  * TMath::Exp(-%f  *x)",p0S,p1S - p1SErr  * 0.25,p2S);
	TString FuncSVXUp = Form("%f - %f  * TMath::Exp(-%f  *x)",p0SVX,p1SVX - p1SVXErr  * 0.5,p2SVX);

	TString FuncNorthDown = Form("%f - %f  * TMath::Exp(-%f  *x)",p0N,p1N+  p1NErr * 0.5,p2N);
	TString FuncSouthDown = Form("%f - %f  * TMath::Exp(-%f  *x)",p0S,p1S + p1SErr * 0.25,p2S);
	TString FuncSVXDown = Form("%f - %f  * TMath::Exp(-%f  *x)",p0SVX,p1SVX + p1SVXErr * 0.5,p2SVX);




	TF1 * FVTXSFitUp = new TF1("FVTXSFitUp",FuncSouthUp.Data());
	TF1 * FVTXNFitUp = new TF1("FVTXNFitUp",FuncNorthUp.Data());
	TF1 * SVXFitUp = new TF1("SVXFitUp",FuncSVXUp.Data());

	TF1 * FVTXSFitDown = new TF1("FVTXSFitDown",FuncSouthDown.Data());
	TF1 * FVTXNFitDown = new TF1("FVTXNFitDown",FuncNorthDown.Data());
	TF1 * SVXFitDown = new TF1("SVXFitDown",FuncSVXDown.Data());
	


	TH1D * FVTXSSystHisUp = new TH1D("FVTXSSystHisUp","",NBinsMult,MultBin2);

	FVTXSSystHisUp->GetXaxis()->SetTitle("FVTX North Event Multiplicity");
	FVTXSSystHisUp->GetYaxis()->SetTitle("J/#psi Trigger Efficiency Systematics (%)");
	FVTXSSystHisUp->GetYaxis()->SetTitleOffset(1.2);	
	FVTXSSystHisUp->GetXaxis()->CenterTitle();
	FVTXSSystHisUp->GetYaxis()->CenterTitle();
	FVTXSSystHisUp->SetTitle("FVXT South ERT 4 #times 4 Trigger Efficiency Systematics Studies");

	TH1D * FVTXNSystHisUp = new TH1D("FVTXNSystHisUp","",NBinsMult,MultBin2);
	
	FVTXNSystHisUp->GetXaxis()->SetTitle("FVTX South Event Multiplicity");
	FVTXNSystHisUp->GetYaxis()->SetTitle("J/#psi Trigger Efficiency Systematics (%)");
	FVTXNSystHisUp->GetYaxis()->SetTitleOffset(1.2);	
	FVTXNSystHisUp->GetXaxis()->CenterTitle();
	FVTXNSystHisUp->GetYaxis()->CenterTitle();
	FVTXNSystHisUp->SetTitle("FVXT North ERT 4 #times 4 Trigger Efficiency Systematics Studies");

	TH1D * SVXSystHisUp = new TH1D("SVXSystHisUp","",NBinsMult,MultBin2);
	
	SVXSystHisUp->GetXaxis()->SetTitle("SVX Event Multiplicity");
	SVXSystHisUp->GetYaxis()->SetTitle("J/#psi Trigger Efficiency Systematics (%)");
	SVXSystHisUp->GetYaxis()->SetTitleOffset(1.2);	
	SVXSystHisUp->GetXaxis()->CenterTitle();
	SVXSystHisUp->GetYaxis()->CenterTitle();
	SVXSystHisUp->SetTitle("SVX ERT 4 #times 4 Trigger Efficiency Systematics Studies");

	TH1D * FVTXSSystHisDown = new TH1D("FVTXSSystHisDown","",NBinsMult,MultBin2);
	TH1D * FVTXNSystHisDown = new TH1D("FVTXNSystHisDown","",NBinsMult,MultBin2);
	TH1D * SVXSystHisDown = new TH1D("SVXSystHisDown","",NBinsMult,MultBin2);


	for(int i = 0; i < NBinsMult; i++){

		BinCenter = 0.5 * (MultBin2[i] + MultBin2[i+1]);
	
		//South

		Nominal = FitSouth->Eval(BinCenter);
		VarUp = FVTXSFitUp->Eval(BinCenter);
		VarDown = FVTXSFitDown->Eval(BinCenter);


		SystUp = (VarUp - Nominal)/Nominal * 100;
		SystDown = (VarDown - Nominal)/Nominal * 100;
	

		FVTXSSystHisUp->SetBinContent(i+1,SystUp);
		FVTXSSystHisUp->SetBinError(i+1,0.0001);
			
		FVTXSSystHisDown->SetBinContent(i+1,SystDown);
		FVTXSSystHisDown->SetBinError(i+1,0.0001);

	
		//North		

		Nominal = FitNorth->Eval(BinCenter);
		VarUp = FVTXNFitUp->Eval(BinCenter);
		VarDown = FVTXNFitDown->Eval(BinCenter);


		SystUp = (VarUp - Nominal)/Nominal * 100;
		SystDown = (VarDown - Nominal)/Nominal * 100;


		FVTXNSystHisUp->SetBinContent(i+1,SystUp);
		FVTXNSystHisUp->SetBinError(i+1,0.0001);
		
		FVTXNSystHisDown->SetBinContent(i+1,SystDown);
		FVTXNSystHisDown->SetBinError(i+1,0.0001);



		//SVX		

		Nominal = FitSVX->Eval(BinCenter);
		VarUp = SVXFitUp->Eval(BinCenter);
		VarDown = SVXFitDown->Eval(BinCenter);


		SystUp = (VarUp - Nominal)/Nominal * 100;
		SystDown = (VarDown - Nominal)/Nominal * 100;

		SVXSystHisUp->SetBinContent(i+1,SystUp);
		SVXSystHisUp->SetBinError(i+1,0.0001);
		
		SVXSystHisDown->SetBinContent(i+1,SystDown);
		SVXSystHisDown->SetBinError(i+1,0.0001);


	}

	

	FVTXSSystHisUp->SetMarkerSize(1);
	FVTXSSystHisUp->SetMarkerStyle(20);
	FVTXSSystHisUp->SetMarkerColor(kRed);
	FVTXSSystHisUp->SetLineColor(kRed);


	FVTXSSystHisDown->SetMarkerSize(1);
	FVTXSSystHisDown->SetMarkerStyle(20);
	FVTXSSystHisDown->SetMarkerColor(kBlue);
	FVTXSSystHisDown->SetLineColor(kBlue);

	FVTXSSystHisUp->SetMinimum(-20);
	FVTXSSystHisUp->SetMaximum(20);
	
	FVTXSSystHisUp->Draw("ep");	
	FVTXSSystHisDown->Draw("epSAME");	



	TLegend* leg = new TLegend(0.52,0.70,0.75,0.90,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.040);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->SetLineWidth(3);
	leg->AddEntry(FVTXSSystHisUp,"Upper Bound","PL");
	leg->AddEntry(FVTXSSystHisDown,"Lower Bound","PL");
	leg->Draw("same");



	
	c->SaveAs("SystPlots/FVTXSouthSystUp.png");


	FVTXNSystHisUp->SetMarkerSize(1);
	FVTXNSystHisUp->SetMarkerStyle(20);
	FVTXNSystHisUp->SetMarkerColor(kRed);
	FVTXNSystHisUp->SetLineColor(kRed);


	FVTXNSystHisDown->SetMarkerSize(1);
	FVTXNSystHisDown->SetMarkerStyle(20);
	FVTXNSystHisDown->SetMarkerColor(kBlue);
	FVTXNSystHisDown->SetLineColor(kBlue);

	FVTXNSystHisUp->SetMinimum(-20);
	FVTXNSystHisUp->SetMaximum(20);
	
	FVTXNSystHisUp->Draw("ep");	
	FVTXNSystHisDown->Draw("epSAME");	

	leg->Draw("same");


	c->SaveAs("SystPlots/FVTXNorthSystUp.png");


	SVXSystHisUp->SetMarkerSize(1);
	SVXSystHisUp->SetMarkerStyle(20);
	SVXSystHisUp->SetMarkerColor(kRed);
	SVXSystHisUp->SetLineColor(kRed);


	SVXSystHisDown->SetMarkerSize(1);
	SVXSystHisDown->SetMarkerStyle(20);
	SVXSystHisDown->SetMarkerColor(kBlue);
	SVXSystHisDown->SetLineColor(kBlue);

	SVXSystHisUp->SetMinimum(-20);	
	SVXSystHisUp->SetMaximum(20);
	
	SVXSystHisUp->Draw("ep");	
	SVXSystHisDown->Draw("epSAME");	
	leg->Draw("same");

	c->SaveAs("SystPlots/SVXNorthSystUp.png");

	
	TFile * fout = new TFile("/sphenix/user/zshi/PHENIX/JPsiAnaCodes/SystStudies/SystFiles/JPsiTrigSyst.root","RECREATE");
	fout->cd();


	FVTXSSystHisUp->Write();
	FVTXSSystHisDown->Write();


	FVTXNSystHisUp->Write();
	FVTXNSystHisDown->Write();


	SVXSystHisUp->Write();
	SVXSystHisDown->Write();



	fout->Close();
}

