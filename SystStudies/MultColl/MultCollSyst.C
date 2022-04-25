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



void MultCollSyst(){




	gStyle->SetOptStat(0);
	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();



	TFile * fin = new TFile("run15pp_anaT0MB_anaT2MU_PDFs_preliminary.root");
	fin->cd();

	TH1D * MultCollData = (TH1D *) fin->Get("rFVTXS_MB_high1000");


	TH1D *rFVTXS02_MB = new TH1D("rFVTXS02_MB","MB: Evt_Mult_FVTXS 2% double/single",20,0,20);
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




	const int NBins = 10;
	double MultBin[NBins + 1] = {0,1,2,3,4,5,6,8,10,12,19};



	TH1D * MultCollModelRebin =  new TH1D("MultCollModelRebin","",NBins,MultBin);
	TH1D * MultCollDataRebin =  new TH1D("MultCollDataRebin","",NBins,MultBin);

	int NTotalBins = rFVTXS02_MB->GetNbinsX();

	float Value;
	float ValueErr;

	float TotalValue;
	float TotalValueErr;

	float TotalWeight;

	float FinalValue;
	float FinalValueErr;

	for(int j = 0; j < NBins; j++){
		TotalValue = 0;	
		TotalValueErr = 0;
		TotalWeight = 0;
		for(int i = MultBin[j]; i < MultBin[j+1]; i++){


			Value = rFVTXS02_MB->GetBinContent(i+1);
			ValueErr = rFVTXS02_MB->GetBinError(i+1);
			TotalValue = TotalValue + Value/(ValueErr * ValueErr);
			TotalValueErr = TotalValueErr + ValueErr/(ValueErr * ValueErr);
			TotalWeight = TotalWeight + 1/(ValueErr * ValueErr);
		

		}

		FinalValue = TotalValue/TotalWeight;
		FinalValueErr = TotalValueErr/TotalWeight/((MultBin[j+1]-MultBin[j])* (MultBin[j+1]-MultBin[j]) );
		MultCollModelRebin->SetBinContent(j+1,FinalValue);
		MultCollModelRebin->SetBinError(j+1,FinalValueErr);

	}


	for(int j = 0; j < NBins; j++){
		TotalValue = 0;	
		TotalValue = 0;
		TotalWeight = 0;
		for(int i = MultBin[j]; i < MultBin[j+1]; i++){


			Value = MultCollData->GetBinContent(i+1);
			ValueErr = MultCollData->GetBinError(i+1);
			TotalValue = TotalValue + Value/(ValueErr * ValueErr);
			TotalValueErr = TotalValueErr + ValueErr/(ValueErr * ValueErr);
			TotalWeight = TotalWeight + 1/(ValueErr * ValueErr);
	
		}
	
		FinalValue = TotalValue/TotalWeight;
		FinalValueErr = TotalValueErr/TotalWeight/((MultBin[j+1]-MultBin[j])* (MultBin[j+1]-MultBin[j]));
		MultCollDataRebin->SetBinContent(j+1,FinalValue);
		MultCollDataRebin->SetBinError(j+1,FinalValueErr);

	}





	//MultCollModelRebin->Scale(1/TotalWeight);





	MultCollDataRebin->SetMarkerStyle(20);
	MultCollDataRebin->SetMarkerColor(3);
	MultCollDataRebin->SetMarkerSize(1);
	MultCollDataRebin->SetLineColor(3);


	MultCollModelRebin->SetMarkerStyle(20);
	MultCollModelRebin->SetMarkerColor(2);
	MultCollModelRebin->SetMarkerSize(1);
	MultCollModelRebin->SetLineColor(2);
	

	MultCollDataRebin->SetMinimum(0.8);
	MultCollDataRebin->SetMaximum(1.3);

	MultCollDataRebin->SetTitle("Double Collision Data-Model Comparison");
	MultCollDataRebin->GetXaxis()->SetTitle("Event Multiplicity");
	MultCollDataRebin->GetYaxis()->SetTitle("Multiple Collision Correction Factor");
	MultCollDataRebin->GetXaxis()->CenterTitle();
	MultCollDataRebin->GetYaxis()->CenterTitle();
	MultCollDataRebin->GetYaxis()->SetTitleOffset(1.4);

	MultCollDataRebin->Draw("ep");
	MultCollModelRebin->Draw("epSAME");



	TLegend* leg = new TLegend(0.12,0.65,0.40,0.85,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.040);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->SetLineWidth(3);
	leg->AddEntry(MultCollModelRebin,"Model","PL");
	leg->AddEntry(MultCollDataRebin,"PHENIX Data: 1000 - 1500 kHz","PL");
	leg->Draw("same");


	c->SaveAs("MultCollComp.png");



	TH1D *  MultCollSyst = (TH1D *) MultCollModelRebin->Clone("MultCollSyst");

	MultCollSyst->Divide(MultCollDataRebin);

	float SystValue;

	for(int i = 0; i < NBins; i++){

		Value = MultCollSyst->GetBinContent(i+1);
		SystValue = abs(Value - 1) * 100;

		MultCollSyst->SetBinContent(i+1,SystValue);		
		MultCollSyst->SetBinError(i+1,0.0001);
	}




	MultCollSyst->SetMarkerStyle(20);
	MultCollSyst->SetMarkerColor(1);
	MultCollSyst->SetMarkerSize(1);
	MultCollSyst->SetLineColor(1);

	MultCollSyst->SetTitle("Double Collision Systematic Error");
	MultCollSyst->GetXaxis()->SetTitle("Event Multiplicity");
	MultCollSyst->GetYaxis()->SetTitle("Systematic Error (%)");

	MultCollSyst->GetXaxis()->CenterTitle();
	MultCollSyst->GetYaxis()->CenterTitle();
	MultCollSyst->GetYaxis()->SetTitleOffset(1.4);

	MultCollSyst->Draw("ep");


	c->SaveAs("MultCollSyst.png");

	

	TFile * fout = new TFile("/sphenix/user/zshi/PHENIX/JPsiAnaCodes/SystStudies/SystFiles/MultCollSyst.root","RECREATE");
	fout->cd();
	MultCollSyst->Write();
	fout->Close();


}
