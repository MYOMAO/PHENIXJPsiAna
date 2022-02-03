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
#include "CMS_lumi.C"

using namespace std;

using std::cout;
using std::endl;
#endif

void DrawFinal(){
	cout << "Pass 0" << endl;

	gStyle->SetOptStat(0);

	const int NCases = 3;

	TString Name[NCases] = {"FVTX_N + FVTX_S","FVTX_N","FVTX_S"};
	cout << "Pass 1" << endl;

	TString HisName[NCases] = {"South + North","North","South"};
	cout << "Pass 2" << endl;

	TString OutName[NCases] = {"North_AND_South","North","South"};

	cout << "Pass 3" << endl;

	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();
	c->SetLeftMargin(0.12);

	TLine * l = new TLine(0,3.0969,19,3.0969);

	l->SetLineWidth(2);
	l->SetLineStyle(2);
	l->SetLineColor(kRed);

	TLine * Unity = new TLine(0,1.0,19,1.0);

	Unity->SetLineWidth(2);
	Unity->SetLineStyle(2);
	Unity->SetLineColor(kRed);


	float Min[NCases] = {3.08,3.08,3.00};


	

	for(int i = 0; i < NCases; i++){

		CMS_lumi(c,19011,0);
		c->Update();


		TFile * fin = new TFile(Form("OutFiles/FitResults_%d.root",i));
		fin->cd();

		TLatex* tex_y;
		tex_y = new TLatex(0.61,0.73,HisName[i].Data());


		tex_y->SetNDC();
		tex_y->SetTextFont(42);
		tex_y->SetTextSize(0.045);
		tex_y->SetLineWidth(2);

		cout << "Pass 1" << endl;

		TH1D * JPsiYield = (TH1D *) fin->Get("JPsiYield");
		JPsiYield->GetXaxis()->SetTitle(Form("Event Multiplicity: %s",Name[i].Data()));
		JPsiYield->GetYaxis()->SetTitle("J/#psi Signal Raw Yield From Fits");
		JPsiYield->GetXaxis()->CenterTitle();
		JPsiYield->GetYaxis()->CenterTitle();
		JPsiYield->GetYaxis()->SetTitleOffset(1.6);

		JPsiYield->SetMarkerStyle(20);
		JPsiYield->SetMarkerSize(1);
		JPsiYield->SetMarkerColor(1);
		JPsiYield->SetLineColor(1);
		cout << "Pass 2" << endl;


		JPsiYield->Draw("ep");
		tex_y->Draw();


		c->SaveAs(Form("FinalResults/Parameters/Sig%s.png",OutName[i].Data()));

		TH1D * JPsiMean = (TH1D *) fin->Get("JPsiMean");
		JPsiMean->GetXaxis()->SetTitle(Form("Event Multiplicity: %s",Name[i].Data()));
		JPsiMean->GetYaxis()->SetTitle("J/#psi Peak Mean From Fits (GeV/c^{2})");
		JPsiMean->GetXaxis()->CenterTitle();
		JPsiMean->GetYaxis()->CenterTitle();
		JPsiMean->GetYaxis()->SetTitleOffset(1.6);

		JPsiMean->SetMarkerStyle(20);
		JPsiMean->SetMarkerSize(1);
		JPsiMean->SetMarkerColor(1);
		JPsiMean->SetLineColor(1);

		JPsiMean->SetMinimum(Min[i]);



		JPsiMean->Draw("ep");

		tex_y->Draw();
		l->Draw("SAME");

		c->SaveAs(Form("FinalResults/Parameters/Mean%s.png",OutName[i].Data()));


	



	}

	TFile * finData = new TFile("merged_anaT2_run15pp_MU_v14.root");

	finData->cd();

	TTree * anaT2 = (TTree *) finData->Get("anaT2");

	const int NBins = 5;
	int MultBin[NBins + 1] = {0,2,5,8,12,19};

	TString Var[NCases] = {"Evt_Mult_FVTXN + Evt_Mult_FVTXS","Evt_Mult_FVTXN","Evt_Mult_FVTXS"};	
	TString EtaCut[NCases] = {"dimuon_rapidity > -999 ","dimuon_rapidity > 0","dimuon_rapidity < 0"};	
	
	int NEventStat[NCases][NBins];
	int NEventTotal[NCases];

	for(int i = 0; i < NCases; i++){


				
		TH1D * MultHisTotal = new TH1D("MultHisTotal","",200,0,20);
		anaT2->Project("MultHisTotal",Var[i].Data(),EtaCut[i].Data());
	
		NEventTotal[i] = MultHisTotal->Integral();

		
	
		for(int j = 0; j < NBins; j++){
		
			TString Cut = Form("%s && %s >= %d && %s < %d",EtaCut[i].Data(),Var[i].Data(),MultBin[j],Var[i].Data(),MultBin[j+1]);
			
			TH1D * MultHis = new TH1D("MultHis","",200,0,20);
			anaT2->Project("MultHis",Var[i].Data(),Cut.Data());
			
			NEventStat[i][j] = MultHis->Integral();


			cout << "NEventStat[i][j]  = " << NEventStat[i][j]  << endl;
		}

		
	
		cout << "NEventTotal[i] = " << NEventTotal[i] << endl;

	}


	for(int i = 0; i < NCases; i++){


		TFile * fin = new TFile(Form("OutFiles/FitResults_%d.root",i));
		fin->cd();

		TLatex* tex_y;
		tex_y = new TLatex(0.61,0.73,HisName[i].Data());


		tex_y->SetNDC();
		tex_y->SetTextFont(42);
		tex_y->SetTextSize(0.045);
		tex_y->SetLineWidth(2);

		TFile * finInc = new TFile(Form("OutFiles/FitResults_Inc_%d.root",i));
		
		TH1D * JPsiYieldInc = (TH1D *) finInc->Get("JPsiYield");

		float IncYield = JPsiYieldInc->GetBinContent(1)/NEventTotal[i];
		float IncYieldError = JPsiYieldInc->GetBinError(1)/NEventTotal[i];
		
	

		TH1D * JPsiYield = (TH1D *) fin->Get("JPsiYield");

		TH1D * JPsiYieldRatio = (TH1D * ) JPsiYield->Clone("JPsiYieldRatio");
		JPsiYieldRatio->GetXaxis()->SetTitle(Form("Event Multiplicity: %s",Name[i].Data()));
		JPsiYieldRatio->GetYaxis()->SetTitle("(J/#psi/N)/(J/#psi_{MB}/N_{MB})");
		JPsiYieldRatio->GetXaxis()->CenterTitle();
		JPsiYieldRatio->GetYaxis()->CenterTitle();
		JPsiYieldRatio->GetYaxis()->SetTitleOffset(1.6);

		JPsiYieldRatio->SetMarkerStyle(20);
		JPsiYieldRatio->SetMarkerSize(1);
		JPsiYieldRatio->SetMarkerColor(1);
		JPsiYieldRatio->SetLineColor(1);

	
		for(int j = 0; j < NBins; j++){

			float Yield = JPsiYield->GetBinContent(j+1)/NEventStat[i][j];
			float YieldError = JPsiYield->GetBinError(j+1)/NEventStat[i][j];

			float Ratio = Yield/IncYield;
			float RatioError = Ratio * TMath::Sqrt(IncYieldError/IncYield * IncYieldError/IncYield + YieldError/Yield * YieldError/Yield);
	

			JPsiYieldRatio->SetBinContent(j+1,Ratio);
			JPsiYieldRatio->SetBinError(j+1,RatioError);


		}

		JPsiYieldRatio->SetMinimum(0.5);
		JPsiYieldRatio->SetMaximum(1.2);

		JPsiYieldRatio->Draw("ep");
		tex_y->Draw("SAME");
		Unity->Draw("SAME");


		c->SaveAs(Form("FinalResults/Ratio/JPsiRatio_%s.png",OutName[i].Data()));

	}





}
