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
using namespace std;

using std::cout;
using std::endl;
#endif


void CalPDFVar(){

	gStyle->SetOptStat(0);

	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();

	const int NBins = 10;
	double MultBin[NBins + 1] = {0,1,2,3,4,5,6,8,10,12,19};

	int NOpt = 2;
	int NMultOpt = 3;


	float SigErr;
	float BackErr;
	float PDFErr;

	float SigVarValue;
	float BackVarValue;
	float NominalValue;


	for(int i =0; i < NOpt; i++){
		for(int j =0; i < NMultOpt; i++){


			TString infile0 = Form("OutFiles/Nominal/FitResults_%d_%d.root",i,j);
			TFile * fin0 = new TFile(infile0.Data()); 	


			TString infile1 = Form("OutFiles/SingleCB/FitResults_%d_%d.root",i,j);
			TFile * fin1 = new TFile(infile1.Data()); 	




			TString infile2 = Form("OutFiles/Linear/FitResults_%d_%d.root",i,j);			
			TFile * fin2 = new TFile(infile2.Data()); 	


			TH1D * JPsiYield0 = (TH1D *) fin0->Get("JPsiYield");
			TH1D * JPsiYield1 = (TH1D *) fin1->Get("JPsiYield");
			TH1D * JPsiYield2 = (TH1D *) fin2->Get("JPsiYield");


			JPsiYield0->Draw("ep");
			JPsiYield0->SetLineColor(kBlack);
			JPsiYield0->SetMarkerColor(kBlack);

			JPsiYield1->Draw("epSAME");
			JPsiYield1->SetLineColor(kRed);
			JPsiYield1->SetMarkerColor(kRed);

			JPsiYield2->Draw("epSAME");
			JPsiYield2->SetLineColor(kBlue);
			JPsiYield2->SetMarkerColor(kBlue);




			TString outfile = Form("Results/Yield/YieldComp_%d_%d.png",i,j);




			TLegend* leg = new TLegend(0.17,0.60,0.50,0.80,NULL,"brNDC");
			leg->SetBorderSize(0);
			leg->SetTextSize(0.040);
			leg->SetTextFont(42);
			leg->SetFillStyle(0);
			leg->SetLineWidth(3);
			leg->AddEntry(JPsiYield0,"Nominal Results","PL");
			leg->AddEntry(JPsiYield1,"Signal: Single Crystal Ball","PL");
			leg->AddEntry(JPsiYield2,"Background: Linear","PL");
			leg->Draw("same");

			c->SaveAs(outfile.Data());

			TH1D * SigError = (TH1D *) fin0->Clone("SigError");
			SigError->GetYaxis()->SetTitle("Systematic Uncertainties (%)");

			TH1D * BackError = (TH1D *) fin0->Clone("BackError");
			BackError->GetYaxis()->SetTitle("Systematic Uncertainties (%)");

			TH1D * TotalError = (TH1D *) fin0->Clone("TotalError");
			TotalError->GetYaxis()->SetTitle("Systematic Uncertainties (%)");
			

			for(int q = 0; q < NBins; q++){


				NominalValue =  JPsiYield0->GetBinContent(q+1);
				SigVarValue =  JPsiYield1->GetBinContent(q+1);
				BackVarValue =  JPsiYield2->GetBinContent(q+1);
				
				SigErr = abs(SigVarValue - NominalValue)/NominalValue * 100;
				BackErr = abs(BackVarValue - NominalValue)/NominalValue * 100;
				PDFErr = sqrt(SigErr * SigErr + BackErr * BackErr);

				SigError->SetBinContent(q+1,SigErr);
				BackError->SetBinContent(q+1,BackErr);
				TotalError->SetBinContent(q+1,PDFErr);
				
			}

			SigError->Draw("ep");
			SigError->SetLineColor(kBlue);
			SigError->SetMarkerColor(kBlue);

			BackError->Draw("epSAME");
			BackError->SetLineColor(kBlue);
			BackError->SetMarkerColor(kBlue);




			TLegend* leg2 = new TLegend(0.17,0.60,0.50,0.80,NULL,"brNDC");
			leg2->SetBorderSize(0);
			leg2->SetTextSize(0.040);
			leg2->SetTextFont(42);
			leg2->SetFillStyle(0);
			leg2->SetLineWidth(3);
			leg2->AddEntry(JPsiYield1,"Signal: Single Crystal Ball","PL");
			leg2->AddEntry(JPsiYield2,"Background: Linear","PL");
			leg2->Draw("same");


			TString outfile2 = Form("Results/SigBakErr/SigBackComp_%d_%d.png",i,j);


			c->SaveAs(outfile2.Data());


			TotalError->Draw("ep");
			TString outfile3 = Form("Results/PDFVar/SigBackComp_%d_%d.png",i,j);

			c->SaveAs(outfile3.Data());
			
		}
	}

	

	


}
