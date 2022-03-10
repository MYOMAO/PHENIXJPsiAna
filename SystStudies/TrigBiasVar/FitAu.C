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



void FitAu(){



	gStyle->SetOptStat(0);
	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();



	TFile * fin = new TFile("infiles/Run15pAu_BBC_TrigEff_Multi_SVX.root");
	fin->cd();

	TH1D * h3SVX = (TH1D *) fin->Get("h3SVX");
	


	h3SVX->GetXaxis()->SetTitle("SVX MB Event Multiplicity");
	h3SVX->GetYaxis()->SetTitle("MB Trigger Bias");
	h3SVX->SetTitle("BBCLL1_Eff 2015 pAu");

	TF1 * FitSVX = new TF1("FitSVX","[0] - [1] * TMath::Exp(-[2] *x)",1,12);


	h3SVX->SetMarkerStyle(20);
	h3SVX->SetMarkerSize(1);
	h3SVX->SetMarkerColor(1);
	h3SVX->SetLineColor(1);



	FitSVX->SetParLimits(0,0.8,1.0);
	FitSVX->SetParLimits(1,0,0.5);
	FitSVX->SetParLimits(2,0,10);



	h3SVX->Fit(FitSVX,"R");
	
	float p0SVX = FitSVX->GetParameter(0);
	float p1SVX = FitSVX->GetParameter(1);
	float p2SVX = FitSVX->GetParameter(2);
	TString FuncSVXhName = Form("%f - %f  * TMath::Exp(-%f  *x)",p0SVX,p1SVX,p2SVX);

	cout << "FuncSVXhName : " << FuncSVXhName.Data() << endl;
	

	h3SVX->Draw();
	c->SaveAs("TrigBias/h3SVX.png");



}

