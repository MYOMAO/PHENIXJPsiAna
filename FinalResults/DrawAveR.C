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



void DrawAveR(int Opt, int MultOpt){




	TString FileName;

	if(Opt == 0) FileName = "North";
	if(Opt == 1) FileName = "South";



	TString TitleName;

	if(Opt == 0) TitleName = "FVTX North Muons";
	if(Opt == 1) TitleName = "FVTX South Muons";



	TString XNameBefore;

	if(MultOpt == 0) XNameBefore = "N_{ch}^{FVTXN}";
	if(MultOpt == 1) XNameBefore = "N_{ch}^{FVTXS}";
	if(MultOpt == 2) XNameBefore = "N_{ch}^{SVX}";


	TString XName;

	if(MultOpt == 0) XName = "N_{ch}^{FVTXN}/<N_{ch}^{FVTXN}>";
	if(MultOpt == 1) XName = "N_{ch}^{FVTXS}/<N_{ch}^{FVTXS}>";
	if(MultOpt == 2) XName = "N_{ch}^{SVX}/<N_{ch}^{SVX}>";




	TString MBHisName;
	if(Opt == 0) MBHisName = "Evt_Mult_FVTXN_His";
	if(Opt == 1) MBHisName = "Evt_Mult_FVTXS_His";
	if(Opt == 2) MBHisName = "Evt_Mult_SVX_His";


	gStyle->SetOptStat(0);
	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();

	TFile * finMU = new TFile(Form("Infiles/FitResults_%d_%d.root",Opt, MultOpt));
	finMU->cd();
	TH1D * JPsiYield = (TH1D *) finMU->Get("JPsiYield"); 

	JPsiYield->GetXaxis()->SetTitle(XNameBefore.Data());
	JPsiYield->GetYaxis()->SetTitle("J/#psi Signal Raw Yield");
	JPsiYield->GetXaxis()->CenterTitle();
	JPsiYield->GetYaxis()->CenterTitle();
	JPsiYield->GetYaxis()->SetTitleOffset(1.4);

	JPsiYield->Draw("ep");
	JPsiYield->SetMarkerStyle(20);
	JPsiYield->SetMarkerSize(1);
	JPsiYield->SetMarkerColor(1);
	JPsiYield->SetLineColor(1);

	c->SaveAs(Form("FinalPlots/RawYields/JPsiYield_%d_%d.png",Opt,MultOpt));

	TFile * finMB = new TFile(Form("Infiles/%sMBEvtMult.root",FileName.Data()));
	finMB->cd();
	TH1D * Evt_Mult_FVTXN_His = (TH1D *) finMB->Get(MBHisName.Data()); 

	float TotalJPsi = JPsiYield->Integral();
	float TotalMB = Evt_Mult_FVTXN_His->Integral();
	float AveJPsi = TotalJPsi/TotalMB;
	float AverageMB = Evt_Mult_FVTXN_His->GetMean();





	Evt_Mult_FVTXN_His->GetXaxis()->SetTitle(XNameBefore.Data());
	Evt_Mult_FVTXN_His->GetYaxis()->SetTitle("Number of Events");
	Evt_Mult_FVTXN_His->GetXaxis()->CenterTitle();
	Evt_Mult_FVTXN_His->GetYaxis()->CenterTitle();
	Evt_Mult_FVTXN_His->GetYaxis()->SetTitleOffset(1.4);

	Evt_Mult_FVTXN_His->SetMarkerStyle(20);
	Evt_Mult_FVTXN_His->SetMarkerSize(1);
	Evt_Mult_FVTXN_His->SetMarkerColor(1);
	Evt_Mult_FVTXN_His->SetLineColor(1);


	TCanvas * c2 = new TCanvas("c2","c2",600,600);	
	c2->cd();
	c2->SetLogy();



	Evt_Mult_FVTXN_His->Draw("ep");


	c2->SaveAs(Form("FinalPlots/MBEvent/MBHisName_%d_%d.png",Opt,MultOpt));



	TH1D * Ratio = (TH1D *) JPsiYield->Clone("Ratio");

	Ratio->Sumw2();
	Evt_Mult_FVTXN_His->Sumw2();

	Ratio->Divide(Evt_Mult_FVTXN_His);

	Ratio->Scale(AverageMB/AveJPsi);



	Ratio->SetTitle(TitleName.Data());
	Ratio->GetYaxis()->SetTitle("<R>");
	Ratio->GetYaxis()->SetTitleOffset(1.4);
	Ratio->GetXaxis()->SetTitle(XName.Data());
	Ratio->GetYaxis()->CenterTitle();
	Ratio->GetXaxis()->CenterTitle();


	Ratio->SetMarkerStyle(20);
	Ratio->SetMarkerSize(1);
	Ratio->SetMarkerColor(1);
	Ratio->SetLineColor(1);

	/*

	   const int NBins = 5;
	   double MultBinScale[NBins + 1] = {0,2.0/AverageMB,5.0/AverageMB,8.0/AverageMB,12.0/AverageMB,19.0/AverageMB};
	   double MultBin[NBins + 1] = {0,2.0,5.0,8.0,12.0,19.0};
	   */
	const int NBins = 10;
	double MultBinScale[NBins + 1] = {0,1/AverageMB,2/AverageMB,3/AverageMB,4/AverageMB,5/AverageMB,6/AverageMB,8/AverageMB,10/AverageMB,12/AverageMB,19/AverageMB};
	double MultBin[NBins + 1] = {0,1,2,3,4,5,6,8,10,12,19};


	float Value;
	float ValueErr;



	TH1D * RatioScaled = new TH1D("RatioScaled","",NBins,MultBinScale);
	RatioScaled->SetTitle(TitleName.Data());
	RatioScaled->GetYaxis()->SetTitle("<R>");
	RatioScaled->GetYaxis()->SetTitleOffset(1.4);
	RatioScaled->GetXaxis()->SetTitle(XName.Data());
	RatioScaled->GetYaxis()->CenterTitle();
	RatioScaled->GetXaxis()->CenterTitle();


	RatioScaled->SetMarkerStyle(20);
	RatioScaled->SetMarkerSize(1);
	RatioScaled->SetMarkerColor(1);
	RatioScaled->SetLineColor(1);
	RatioScaled->SetMinimum(0);

	for(int i = 0; i < NBins; i++){


		Value =  Ratio->GetBinContent(i+1);
		ValueErr =  Ratio->GetBinError(i+1);


		RatioScaled->SetBinContent(i+1,Value);
		RatioScaled->SetBinError(i+1,ValueErr);

	}
	c->cd();

	RatioScaled->Draw("ep");

	TF1 * func = new TF1("Func","x",0,19);
	func->SetLineStyle(2);
	func->Draw("SAME");

	c->SaveAs(Form("FinalPlots/UnCorr/JPsiRatio_%d_%d.png",Opt,MultOpt));


	TString FuncForm; 


	if(MultOpt == 0) FuncForm = "0.906048 - 0.433618  * TMath::Exp(-0.283628  *x)";
	if(MultOpt == 1) FuncForm = "0.866227 - 0.386111  * TMath::Exp(-0.286963  *x)";
	if(MultOpt == 2) FuncForm = "0.985577 - 0.412983  * TMath::Exp(-0.196578  *x)";

	TF1 * TrigBias = new TF1("TrigBias",FuncForm.Data(),0,19);


	TH1D * RatioScaledCorr = new TH1D("RatioScaledCorr","",NBins,MultBinScale);
	RatioScaledCorr->SetTitle(TitleName.Data());
	RatioScaledCorr->GetYaxis()->SetTitle("<R> With Trig Bias Corrected");
	RatioScaledCorr->GetYaxis()->SetTitleOffset(1.2);
	RatioScaledCorr->GetXaxis()->SetTitle(XName.Data());
	RatioScaledCorr->GetYaxis()->CenterTitle();
	RatioScaledCorr->GetXaxis()->CenterTitle();


	RatioScaledCorr->SetMarkerStyle(20);
	RatioScaledCorr->SetMarkerSize(1);
	RatioScaledCorr->SetMarkerColor(1);
	RatioScaledCorr->SetLineColor(1);

	RatioScaledCorr->SetMinimum(0);

	float R;
	float RErr;

	float CorrFactor;
	float MultFactor;

	float RFinal;
	float RErrFinal;

	float BinCenter;



	TF1 * MultiFunc = new TF1("MultiFunc","9.83665e-01 + 6.44703e-03 * x + 3.02238e-04 * x * x + 4.83052e-06 * x * x * x * x",0,19);

	
	//Box for Systematics

	double MultBinCenter[NBins] = {0.5/AverageMB,1.5/AverageMB,2.5/AverageMB,3.5/AverageMB,4.5/AverageMB,5.5/AverageMB,7/AverageMB,9/AverageMB,11/AverageMB,15.5/AverageMB};
	double MultBinHigh[NBins] = {0.5/AverageMB,0.5/AverageMB,0.5/AverageMB,0.5/AverageMB,0.5/AverageMB,0.5/AverageMB,1/AverageMB,1/AverageMB,1/AverageMB,3.5/AverageMB};
	double MultBinLow[NBins] = {0.5/AverageMB,0.5/AverageMB,0.5/AverageMB,0.5/AverageMB,0.5/AverageMB,0.5/AverageMB,1/AverageMB,1/AverageMB,1/AverageMB,3.5/AverageMB};




	double YCenter[NBins];
	double YUp[NBins];
	double YDown[NBins];


	for(int i = 0; i < NBins; i++){


		R = RatioScaled->GetBinContent(i+1);
		RErr = RatioScaled->GetBinError(i+1);
		BinCenter = (MultBin[i] + MultBin[i+1]) * 0.5;

		CorrFactor = TrigBias->Eval(BinCenter);
		MultFactor = MultiFunc->Eval(BinCenter);

		//		cout << "BinCenter = " << BinCenter << "   MultFactor = " << MultFactor << endl;
		//		cout << "CorrFactor = " << CorrFactor << endl;

		RFinal = R * CorrFactor * MultFactor;
		RErrFinal = RErr * CorrFactor * MultFactor;


		RatioScaledCorr->SetBinContent(i+1,RFinal);
		RatioScaledCorr->SetBinError(i+1,RErrFinal);
	
		YCenter[i] = RFinal;
		YUp[i] = RFinal * 0.342;
		YDown[i] = RFinal * 0.342;

	}


	TGraphAsymmErrors *Syst = new TGraphAsymmErrors(NBins, MultBinCenter, YCenter,
			MultBinLow, MultBinHigh,
			YDown,YUp);

	Syst->SetFillColorAlpha(kBlue-9,0.5);
	Syst->SetLineColor(kBlue-9);


	c->cd();



	RatioScaledCorr->Draw("ep");



	c->SaveAs(Form("FinalPlots/Corr/JPsiRatio_%d_%d.png",Opt,MultOpt));



	TF1 * FitFunc = new TF1("FitFunc","[0] + [1] *x",0.6,7);
	FitFunc->SetParLimits(0,-2,2);
	FitFunc->SetParLimits(1,0.2,2.5);
	FitFunc->SetLineColor(kGreen);
	FitFunc->SetLineStyle(2);
	FitFunc->SetLineWidth(2);


	RatioScaledCorr->Draw("ep");
	RatioScaledCorr->Fit(FitFunc,"R");

	
	RatioScaledCorr->SetMaximum(RFinal * 1.5);

	if(RFinal < 10) RatioScaledCorr->SetMaximum(20);
	if(RFinal > 20) RatioScaledCorr->SetMaximum(15);


	RatioScaledCorr->Draw("ep");
	FitFunc->Draw("SAME");
	func->Draw("SAME");

	float p0 = FitFunc->GetParameter(0);
	float p1 = FitFunc->GetParameter(1);
	float XInter = -p0/p1;

	float p0Err = FitFunc->GetParError(0);
	float p1Err = FitFunc->GetParError(1);
	float XInterErr = (p0/p1) * sqrt(p0Err/p0 * p0Err/p0 + p1Err/p1 * p1Err/p1);


	TLegend* leg = new TLegend(0.17,0.60,0.50,0.80,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.040);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->SetLineWidth(3);
	leg->AddEntry(RatioScaledCorr,"PHENIX Data","PL");
	leg->AddEntry(FitFunc,Form("Fit to Data: y = %.2f + %.2f x",p0,p1),"L");	
	leg->AddEntry(func,Form("Reference: y = x"),"L");	
	leg->Draw("same");



	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.035);

	lat->DrawLatex(0.18,0.55,Form("x-intercept = %.2f",XInter));	

	Syst->Draw("5SAME");

	c->SaveAs(Form("FinalPlots/Final/JPsiRatio_%d_%d.png",Opt,MultOpt));

	TString OutFileName = Form("OutFiles/JPsiR_%d_%d.root",Opt,MultOpt);

	TFile * fout = new TFile(OutFileName.Data(),"RECREATE");
	fout->cd();


	RatioScaledCorr->Write();
	FitFunc->Write();
	Syst->Write();

	fout->Close();

	cout << "---------------------------------------------------------------------------------------------" << endl;
	cout << "m = " << p0 << " \\pm " << p0Err << endl;
	cout << "n = " << p1 << " \\pm " << p1Err << endl;
	cout << "XInter = " << XInter << " \\pm " << XInterErr << endl;
	
	cout << "---------------------------------------------------------------------------------------------" << endl;


}
