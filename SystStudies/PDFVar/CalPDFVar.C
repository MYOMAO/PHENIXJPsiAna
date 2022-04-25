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
	int MultBin[NBins + 1] = {0,1,2,3,4,5,6,8,10,12,19};


	const int NOpt = 2;
	const int NMultOpt = 3;


	float SigErr;
	float BackErr;
	float PDFErr;

	float SigVarValue;
	float BackVarValue;
	float NominalValue;


	float PDFSyst[NOpt][NMultOpt][NBins];

	TString TitileNameDir[NOpt] = {"North","South"};
	TString TitileNameMult[NMultOpt] = {"FVTXN","FVTXS","SVX"};

	std::ofstream OutText("SystText.txt");
	cout << "Pass 1" << endl;


	TH1D * TotalError[NOpt][NMultOpt];

	for(int i =0; i < NOpt; i++){
		for(int j =0; j < NMultOpt; j++){


			TString infile0 = Form("OutFiles/Nominal/FitResults_%d_%d.root",i,j);
			TFile * fin0 = new TFile(infile0.Data()); 	


			TString infile1 = Form("OutFiles/SingleCB/FitResults_%d_%d.root",i,j);
			TFile * fin1 = new TFile(infile1.Data()); 	




			TString infile2 = Form("OutFiles/Linear/FitResults_%d_%d.root",i,j);			
			TFile * fin2 = new TFile(infile2.Data()); 	


			cout << "Pass 2" << endl;

			TH1D * JPsiYield0 = (TH1D *) fin0->Get("JPsiYield");
			TH1D * JPsiYield1 = (TH1D *) fin1->Get("JPsiYield");
			TH1D * JPsiYield2 = (TH1D *) fin2->Get("JPsiYield");
	
	
			TString TitleName = Form("%s J/#psi with %s Multiplicity",TitileNameDir[i].Data(),TitileNameMult[j].Data());

			JPsiYield0->GetXaxis()->SetTitle(TitleName.Data());
			JPsiYield0->GetYaxis()->SetTitle("J/#psi Signal Raw Yield");
			JPsiYield0->GetXaxis()->CenterTitle();
			JPsiYield0->GetYaxis()->CenterTitle();
	

			JPsiYield1->GetXaxis()->SetTitle(TitleName.Data());
			JPsiYield1->GetYaxis()->SetTitle("J/#psi Signal Raw Yield");
			JPsiYield1->GetXaxis()->CenterTitle();
			JPsiYield1->GetYaxis()->CenterTitle();

			JPsiYield2->GetXaxis()->SetTitle(TitleName.Data());
			JPsiYield2->GetYaxis()->SetTitle("J/#psi Signal Raw Yield");
			JPsiYield2->GetXaxis()->CenterTitle();
			JPsiYield2->GetYaxis()->CenterTitle();

			JPsiYield0->SetLineColor(kBlack);
			JPsiYield0->SetMarkerColor(kBlack);
			JPsiYield0->SetMarkerStyle(20);
			JPsiYield0->SetMarkerSize(1);
			JPsiYield0->Draw("ep");

			JPsiYield1->SetLineColor(kRed);
			JPsiYield1->SetMarkerColor(kRed);
			JPsiYield1->SetMarkerStyle(20);
			JPsiYield1->SetMarkerSize(1);
			JPsiYield1->Draw("epSAME");


			JPsiYield2->SetLineColor(kBlue);
			JPsiYield2->SetMarkerColor(kBlue);
			JPsiYield2->SetMarkerStyle(20);
			JPsiYield2->SetMarkerSize(1);
			JPsiYield2->Draw("epSAME");



			cout << "Pass 3" << endl;


			TString outfile = Form("Results/Yield/YieldComp_%d_%d.png",i,j);




			TLegend* leg = new TLegend(0.40,0.60,0.80,0.80,NULL,"brNDC");
			leg->SetBorderSize(0);
			leg->SetTextSize(0.040);
			leg->SetTextFont(42);
			leg->SetFillStyle(0);
			leg->SetLineWidth(3);
			leg->AddEntry(JPsiYield0,"Nominal Results","PL");
			leg->AddEntry(JPsiYield1,"Signal: Johnson SU Function","PL");
			leg->AddEntry(JPsiYield2,"Background: Linear","PL");
			leg->Draw("same");

			cout << "Pass 4" << endl;


			c->SaveAs(outfile.Data());

			TH1D * SigError = (TH1D *) JPsiYield0->Clone("SigError");
			SigError->GetYaxis()->SetTitle("Systematic Uncertainties (%)");

			TH1D * BackError = (TH1D *) JPsiYield0->Clone("BackError");
			BackError->GetYaxis()->SetTitle("Systematic Uncertainties (%)");

			TotalError[i][j] = (TH1D *) JPsiYield0->Clone(Form("TotalError_%d_%d",i,j));
			TotalError[i][j]->GetYaxis()->SetTitle("Systematic Uncertainties (%)");


			cout << "Pass 5" << endl;


			for(int q = 0; q < NBins; q++){


				NominalValue =  JPsiYield0->GetBinContent(q+1);
				SigVarValue =  JPsiYield1->GetBinContent(q+1);
				BackVarValue =  JPsiYield2->GetBinContent(q+1);

				SigErr = abs(SigVarValue - NominalValue)/NominalValue * 100;
				BackErr = abs(BackVarValue - NominalValue)/NominalValue * 100;
				PDFErr = sqrt(SigErr * SigErr + BackErr * BackErr);

				SigError->SetBinContent(q+1,SigErr);
				BackError->SetBinContent(q+1,BackErr);
				TotalError[i][j]->SetBinContent(q+1,PDFErr);


				SigError->SetBinError(q+1,SigErr/20);
				BackError->SetBinError(q+1,BackErr/20);
				TotalError[i][j]->SetBinError(q+1,PDFErr/20);


			}

			for(int q = 0; q < NBins; q++){


				PDFSyst[i][j][q] = TotalError[i][j]->GetBinContent(q+1);


			}

			cout << "Pass 6" << endl;


			SigError->Draw("ep");
			SigError->SetLineColor(kRed);
			SigError->SetMarkerColor(kRed);

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

			cout << "Pass 7" << endl;


			TString outfile2 = Form("Results/SigBakErr/SigBackComp_%d_%d.png",i,j);


			c->SaveAs(outfile2.Data());

			for(int q = 0; q < NBins; q++){


				PDFSyst[i][j][q] = TotalError[i][j]->GetBinContent(q+1);


			}


	
			TotalError[i][j]->Draw("ep");
			TString outfile3 = Form("Results/PDFVar/PDFVar_%d_%d.png",i,j);





			c->SaveAs(outfile3.Data());
		
				
		}
	}


	OutText.precision(3);


	for(int q = 0; q < NBins; q++){
		OutText <<  Form("%d - %d",MultBin[q],MultBin[q+1]) << "  & " <<	PDFSyst[0][0][q]   << "\\%  & " << PDFSyst[0][1][q] << "\\%  & " << PDFSyst[0][2][q] << "\\%  & " <<	PDFSyst[1][0][q] << "\\%  & " << PDFSyst[1][1][q]<< "\\%  & " << PDFSyst[1][2][q]<< "% \\\\ "  << endl;


	}



	TFile * fout = new TFile("/sphenix/user/zshi/PHENIX/JPsiAnaCodes/SystStudies/SystFiles/FitPDFSyst.root","RECREATE");
	fout->cd();


	for(int i = 0; i < NOpt; i++){

		for(int j = 0; j < NMultOpt; j++){

			TotalError[i][j]->Write();

		}

	}



	fout->Close();
	

	



}
