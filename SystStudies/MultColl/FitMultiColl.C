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

void FitMultiColl(){


	gStyle->SetOptStat(0);
	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();

	TH1F *rFVTXS02_MB = new TH1F("rFVTXS02_MB","MB: Evt_Mult_FVTXS 2% double/single",20,0,20);
	rFVTXS02_MB->SetBinContent(1,0.9869716);
	rFVTXS02_MB->SetBinContent(2,0.993943);
	rFVTXS02_MB->SetBinContent(3,1.002147);
	rFVTXS02_MB->SetBinContent(4,1.010515);
	rFVTXS02_MB->SetBinContent(5,1.020512);
	rFVTXS02_MB->SetBinContent(6,1.032645);
	rFVTXS02_MB->SetBinContent(7,1.047229);
	rFVTXS02_MB->SetBinContent(8,1.064575);
	rFVTXS02_MB->SetBinContent(9,1.085827);
	rFVTXS02_MB->SetBinContent(10,1.112288);
	rFVTXS02_MB->SetBinContent(11,1.141278);
	rFVTXS02_MB->SetBinContent(12,1.181665);
	rFVTXS02_MB->SetBinContent(13,1.2263);
	rFVTXS02_MB->SetBinContent(14,1.290002);
	rFVTXS02_MB->SetBinContent(15,1.358501);
	rFVTXS02_MB->SetBinContent(16,1.447069);
	rFVTXS02_MB->SetBinContent(17,1.525226);
	rFVTXS02_MB->SetBinContent(18,1.56107);
	rFVTXS02_MB->SetBinContent(19,1.858142);
	rFVTXS02_MB->SetBinContent(20,2.152994);
	rFVTXS02_MB->SetBinContent(21,0.98);
	rFVTXS02_MB->SetBinError(1,0.0004337521);
	rFVTXS02_MB->SetBinError(2,0.0005053653);
	rFVTXS02_MB->SetBinError(3,0.0006391299);
	rFVTXS02_MB->SetBinError(4,0.0008124137);
	rFVTXS02_MB->SetBinError(5,0.001066722);
	rFVTXS02_MB->SetBinError(6,0.001440912);
	rFVTXS02_MB->SetBinError(7,0.001989862);
	rFVTXS02_MB->SetBinError(8,0.002796996);
	rFVTXS02_MB->SetBinError(9,0.004009258);
	rFVTXS02_MB->SetBinError(10,0.005862871);
	rFVTXS02_MB->SetBinError(11,0.008603301);
	rFVTXS02_MB->SetBinError(12,0.01304129);
	rFVTXS02_MB->SetBinError(13,0.01981832);
	rFVTXS02_MB->SetBinError(14,0.03122383);
	rFVTXS02_MB->SetBinError(15,0.04905561);
	rFVTXS02_MB->SetBinError(16,0.07895244);
	rFVTXS02_MB->SetBinError(17,0.1239293);
	rFVTXS02_MB->SetBinError(18,0.183404);
	rFVTXS02_MB->SetBinError(19,0.3656912);
	rFVTXS02_MB->SetBinError(20,0.6828727);
	rFVTXS02_MB->SetBinError(21,0.3578454);
	rFVTXS02_MB->SetMinimum(0.80);
	rFVTXS02_MB->SetMaximum(3.00);
	rFVTXS02_MB->SetEntries(955.4043);
	
	
	rFVTXS02_MB->Draw("ep");




	rFVTXS02_MB->SetMarkerStyle(20);
	rFVTXS02_MB->SetMarkerSize(1);
	rFVTXS02_MB->SetMarkerColor(1);
	rFVTXS02_MB->SetLineColor(1);


	TF1 * func = new TF1("func","[0] + [1] * x + [2] * x * x + [3] * x * x *x *x ",0,20);

	func->SetParLimits(0,0.9,1.1);
	func->SetParLimits(1,0,10);
	func->SetParLimits(2,0,10);
	func->SetParLimits(3,0,10);
	
	rFVTXS02_MB->Fit(func,"R");
	

	rFVTXS02_MB->GetXaxis()->SetTitle("Event Multiplicity");
	rFVTXS02_MB->GetYaxis()->SetTitle("Double Collision Correction Factor");

	rFVTXS02_MB->GetXaxis()->CenterTitle();
	rFVTXS02_MB->GetYaxis()->CenterTitle();
	rFVTXS02_MB->GetYaxis()->SetTitleOffset(1.4);

	rFVTXS02_MB->Draw("ep");
	func->Draw("SAME");


	c->SaveAs("FitMultiColl.png");


	
	

}
