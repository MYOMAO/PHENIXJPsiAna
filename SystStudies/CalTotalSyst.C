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




void CalTotalSyst(){


	gStyle->SetOptStat(0);

	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();

	const int NOpt = 2;
	const int NOptMult = 3;

	const int NBins = 10;
	double MultBin[NBins + 1] = {0,1,2,3,4,5,6,8,10,12,19};




	TFile * FitSyst = new TFile("SystFiles/FitPDFSyst.root");
	FitSyst->cd();


	TFile * MBTrigSyst = new TFile("SystFiles/MBTrigSyst.root");
	MBTrigSyst->cd();

	TFile * JPsiTrigSyst = new TFile("SystFiles/JPsiTrigSyst.root");
	MBTrigSyst->cd();


	TFile * MultCollSyst = new TFile("SystFiles/MultCollSyst.root");
	MBTrigSyst->cd();


	TH1D* FinalErrorUp[NOpt][NOptMult]; 
	TH1D* FinalErrorDown[NOpt][NOptMult]; 

	TString JPsiType[NOpt] = {"FVTXNorth","FVTXSouth"};

	TString FVTXName[NOptMult] = {"FVTXN","FVTXS","SVX"};


	float SystValueFit;
	float SystValueMultCorr;
	float SystValueReco = 5;
	float SystValueMBTrigUp;
	float SystValueMBTrigDown;
	float SystValueJPsiTrigUp;
	float SystValueJPsiTrigDown;
	float SystValueMBTrigSym;

	float SystValueTotalUp;
	float SystValueTotalDown;

	for(int i = 0; i < NOpt; i++){

		for(int j = 0; j < NOptMult; j++){

			FinalErrorDown[i][j] = new TH1D(Form("FinalErrorDown_%d_%d",i,j),"",NBins,MultBin);
			FinalErrorDown[i][j]->GetYaxis()->SetTitle("Total Systematic Error #sigma_{Total} (%)");
			FinalErrorDown[i][j]->GetXaxis()->SetTitle(Form("%s Event Multiplicity",FVTXName[j].Data()));
			FinalErrorDown[i][j]->SetTitle(Form("%s Multiplicity Systematics (%)",FVTXName[j].Data()));
			FinalErrorDown[i][j]->GetXaxis()->CenterTitle();
			FinalErrorDown[i][j]->GetYaxis()->CenterTitle();
			FinalErrorDown[i][j]->GetYaxis()->SetTitleOffset(1.3);


			FinalErrorDown[i][j]->SetLineColor(kBlue);
			FinalErrorDown[i][j]->SetLineWidth(2);



			FinalErrorUp[i][j] = new TH1D(Form("FinalErrorUp_%d_%d",i,j),"",NBins,MultBin);
			FinalErrorUp[i][j]->GetYaxis()->SetTitle("Total Systematic Error #sigma_{Total} (%)");
			FinalErrorUp[i][j]->GetXaxis()->SetTitle(Form("%s Event Multiplicity",FVTXName[j].Data()));
			FinalErrorUp[i][j]->SetTitle(Form("%s Multiplicity Systematics (%)",FVTXName[j].Data()));
			FinalErrorUp[i][j]->GetXaxis()->CenterTitle();
			FinalErrorUp[i][j]->GetYaxis()->CenterTitle();
			FinalErrorUp[i][j]->GetYaxis()->SetTitleOffset(1.3);


			FinalErrorUp[i][j]->SetLineColor(kRed);
			FinalErrorUp[i][j]->SetLineWidth(2);

		}
	}


	TH1D * MultCollSystHis = (TH1D *) MultCollSyst->Get("MultCollSyst"); 
	MultCollSystHis->SetBinContent(1,0);
	/*
	   for(int i = 0; i < NOpt; i++){

	   for(int j = 0; j < NOptMult; j++){

	   TH1D * FitSystError = (TH1D * ) FitSyst->Get(Form("TotalError_%d_%d",i,j));
	   TH1D * MBSystUp = (TH1D *) MBTrigSyst->Get(Form("TrigBiasSystUp_%d",j));
	   TH1D * MBSystDown = (TH1D *) MBTrigSyst->Get(Form("TrigBiasSystDown_%d",j));
	   TH1D * JPsiSystUp = (TH1D *) JPsiTrigSyst->Get(Form("%sSystHisUp",FVTXName[j].Data()));
	   TH1D * JPsiSystDown = (TH1D *) JPsiTrigSyst->Get(Form("%sSystHisDown",FVTXName[j].Data()));

	   for(int q = 0; q < NBins; q++){

	   SystValueFit = FitSystError->GetBinContent(q+1);
	   SystValueMultCorr =	MultCollSystHis->GetBinContent(q+1);
	   SystValueMBTrigUp = MBSystUp->GetBinContent(q+1);
	   SystValueMBTrigDown = MBSystDown->GetBinContent(q+1);
	   SystValueJPsiTrigUp = JPsiSystUp->GetBinContent(q+1);
	   SystValueJPsiTrigDown = JPsiSystDown->GetBinContent(q+1);


	   SystValueTotalUp = TMath::Sqrt(SystValueFit * SystValueFit + SystValueMultCorr * SystValueMultCorr + SystValueReco * SystValueReco + SystValueMBTrigUp * SystValueMBTrigUp + SystValueJPsiTrigUp * SystValueJPsiTrigUp);

	   SystValueTotalDown = -TMath::Sqrt(SystValueFit * SystValueFit + SystValueMultCorr * SystValueMultCorr + SystValueReco * SystValueReco + SystValueMBTrigDown * SystValueMBTrigDown + SystValueJPsiTrigDown * SystValueJPsiTrigDown);


	   FinalErrorDown[i][j]->SetBinContent(q+1,SystValueTotalDown);
	   FinalErrorUp[i][j]->SetBinContent(q+1,SystValueTotalUp);




	   }




	   }
	   */


	float PrintMB[NBins][NOptMult];
	float PrintJPsiUp[NBins][NOptMult];
	float PrintJPsiDown[NBins][NOptMult];	
	float PrintColl[NBins][NOptMult];
	float PrintTotalUp[NBins][NOptMult];
	float PrintTotalDown[NBins][NOptMult];

	for(int i = 0; i < NOpt; i++){

		for(int j = 0; j < NOptMult; j++){

			TH1D * FitSystError = (TH1D * ) FitSyst->Get(Form("TotalError_%d_%d",i,j));
			TH1D * MBSystUp = (TH1D *) MBTrigSyst->Get(Form("TrigBiasSystUp_%d",j));
			TH1D * MBSystDown = (TH1D *) MBTrigSyst->Get(Form("TrigBiasSystDown_%d",j));
			TH1D * JPsiSystUp = (TH1D *) JPsiTrigSyst->Get(Form("%sSystHisUp",FVTXName[j].Data()));
			TH1D * JPsiSystDown = (TH1D *) JPsiTrigSyst->Get(Form("%sSystHisDown",FVTXName[j].Data()));	
			TH1D * MBSystSym = (TH1D *) MBTrigSyst->Get(Form("TrigBiasSystSym_%d",j));

			for(int q = 0; q < NBins; q++){

				SystValueFit = FitSystError->GetBinContent(q+1);
				SystValueMultCorr =	MultCollSystHis->GetBinContent(q+1);
				SystValueMBTrigUp = MBSystUp->GetBinContent(q+1);
				SystValueMBTrigDown = MBSystDown->GetBinContent(q+1);
				SystValueJPsiTrigUp = JPsiSystUp->GetBinContent(q+1);
				SystValueJPsiTrigDown = JPsiSystDown->GetBinContent(q+1);		
				SystValueMBTrigSym = MBSystSym->GetBinContent(q+1);



				SystValueTotalUp = TMath::Sqrt(SystValueMultCorr * SystValueMultCorr + SystValueReco * SystValueReco + SystValueMBTrigSym * SystValueMBTrigSym + SystValueJPsiTrigUp * SystValueJPsiTrigUp);
				SystValueTotalDown = -TMath::Sqrt(SystValueMultCorr * SystValueMultCorr + SystValueReco * SystValueReco + SystValueMBTrigSym * SystValueMBTrigSym + SystValueJPsiTrigDown * SystValueJPsiTrigDown);


				FinalErrorDown[i][j]->SetBinContent(q+1,SystValueTotalDown);
				FinalErrorUp[i][j]->SetBinContent(q+1,SystValueTotalUp);
	



				PrintMB[q][j] = SystValueMBTrigSym;
				PrintJPsiUp[q][j] = SystValueJPsiTrigUp;
				PrintJPsiDown[q][j] = -SystValueJPsiTrigDown;
				PrintColl[q][j] = SystValueMultCorr;
				PrintTotalUp[q][j] = SystValueTotalUp;
				PrintTotalDown[q][j] = -SystValueTotalDown;


			}

		}
	}
	TFile * fout = new TFile("SystTotalOut.root","RECREATE");	
	fout->cd();


	TLegend* leg = new TLegend(0.13,0.75,0.40,0.90,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.040);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->SetLineWidth(3);
	leg->AddEntry(FinalErrorUp[0][0],"Upper Bound","PL");	
	leg->AddEntry(FinalErrorDown[0][0],"Lower Bound","PL");


	for(int i = 0; i < NOpt; i++){

		for(int j = 0; j < NOptMult; j++){

			FinalErrorDown[i][j]->SetMinimum(-30);
			FinalErrorDown[i][j]->SetMaximum(30);

			FinalErrorDown[i][j]->Write();
			FinalErrorUp[i][j]->Write();


			FinalErrorDown[i][j]->Draw("hist");
			FinalErrorUp[i][j]->Draw("histSAME");

			leg->Draw("same");


			c->SaveAs(Form("FinalPlots/%s%s.png",JPsiType[i].Data(),FVTXName[j].Data()));


		}

	}

	fout->Close();


	std::cout.precision(3);

	for(int j = 0; j < NOptMult; j++){

		cout << "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;

		cout << "$\\sigma{\\epsilon_{trig}^{MB}}$ " << "  & " << PrintMB[0][j] << "\\%  & " << PrintMB[1][j]  << "\\%  & " << PrintMB[2][j] << "\\%  & " << PrintMB[3][j] << "\\%  & " << PrintMB[4][j]  << "\\%  & "  << PrintMB[5][j] << "\\%  & " << PrintMB[6][j]  << "\\%  & "  << PrintMB[7][j] << "\\%  & " << PrintMB[8][j]  << "\\%  & " << PrintMB[9][j] <<  "  //" <<endl;
		cout << "$\\sigma{\\epsilon_{trig}^{J/\psi}}$ Up ($+$) " << " & " << PrintJPsiUp[0][j] << "% & " << PrintJPsiUp[1][j]  << "\\%  & " << PrintJPsiUp[2][j] << "\\%  & " << PrintJPsiUp[3][j] << "\\%  & " << PrintJPsiUp[4][j]  << "\\%  & "  << PrintJPsiUp[5][j] << "\\%  & " << PrintJPsiUp[6][j]  << "\\%  & "  << PrintJPsiUp[7][j] << "\\%  & " << PrintJPsiUp[8][j]  << "\\% & " << PrintJPsiUp[9][j]<<  "  //" <<endl;
		cout <<"$\\sigma{\\epsilon_{trig}^{J/\psi}}$ Down ($-$) " << " & " << PrintJPsiDown[0][j] << "\\%  & " << PrintJPsiDown[1][j]  << "\\%  & " << PrintJPsiDown[2][j] << "\\%  & " << PrintJPsiDown[3][j] << "\\%  & " << PrintJPsiDown[4][j]  << "\\%  & "  << PrintJPsiDown[5][j] << "\\%  & " << PrintJPsiDown[6][j]  << "\\%  & "  << PrintJPsiDown[7][j] << "\\%  & " << PrintJPsiDown[8][j]  << "\\%  & " << PrintJPsiDown[9][j] <<  "  //"<<endl;
		cout << "$\\sigma{\\epsilon_{reco}^{J/\psi}}$ & 5\\% & 5\\%  &  5\\%  &  5\\%  &  5\\%  &  5\\%  &  5\\%  &  5\\%  &  5\\%  &  5\\%  \\" << endl;
		cout << "$\\sigma{f_{coll}}$  " << " & " << PrintColl[0][j] << "\\%  & " << PrintColl[1][j]  << "\\%  & " << PrintColl[2][j] << "\\%  & " << PrintColl[3][j] << "\\%  & " << PrintColl[4][j]  << "\\%  & "  << PrintColl[5][j] << "\\%  & " << PrintColl[6][j]  << "\\%  & "  << PrintColl[7][j] << "\\%  & " << PrintColl[8][j]  << "\\%  & " << PrintColl[9][j] <<  "  //"<<endl;
		cout << "$\\sigma_{Total}$ Up ($+$) " << " & " << PrintTotalUp[0][j] << "\\%  & " << PrintTotalUp[1][j]  << "\\%  & " << PrintTotalUp[2][j] << "\\%  & " << PrintTotalUp[3][j] << "\\%  & " << PrintTotalUp[4][j]  << "\\%  & "  << PrintTotalUp[5][j] << "\\%  & " << PrintTotalUp[6][j]  << "\\%  & "  << PrintTotalUp[7][j] << "\\%  & " << PrintTotalUp[8][j]  << "\\%  & " << PrintTotalUp[9][j]<<  "  //" <<endl;
		cout << "$\\sigma_{Total}$ Down ($-$) "  << " & " << PrintTotalDown[0][j] << "\\%  & " << PrintTotalDown[1][j]  << "\\%  & " << PrintTotalDown[2][j] << "\\%  & " << PrintTotalDown[3][j] << "\\%  & " << PrintTotalDown[4][j]  << "\\%  & "  << PrintTotalDown[5][j] << "\\%  & " << PrintTotalDown[6][j]  << "\\%  & "  << PrintTotalDown[7][j] << "\\%  & " << PrintTotalDown[8][j]  << "\\%  & " << PrintTotalDown[9][j]<<  "  //" <<endl;


		cout << "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;

		
	}






}
