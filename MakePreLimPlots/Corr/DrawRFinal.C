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



void DrawRFinal(int Opt, int MultOpt){



	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.034);
	

	TString FileName;

	if(Opt == 0) FileName = "North";
	if(Opt == 1) FileName = "South";

	const int NDir = 2;
	TString DataName[NDir]={"1.2 < y_{J/#psi} < 2.2","-2.2 < y_{J/#psi} < -1.2"};
	TString Together = "(J/#psi #rightarrow #mu^{+} #mu^{-} subtracted)";

	
	TString TitleName;

	if(Opt == 0) TitleName = "FVTX North Muons";
	if(Opt == 1) TitleName = "FVTX South Muons";


	TString FVTXName;

	if(Opt == 0) FVTXName = "FVTX North";
	if(Opt == 1) FVTXName = "FVTX South";



	TString XNameBefore;

	if(MultOpt == 0) XNameBefore = "N_{ch}^{FVTXN}";
	if(MultOpt == 1) XNameBefore = "N_{ch}^{FVTXS}";
	if(MultOpt == 2) XNameBefore = "N_{ch}^{SVX}";


	TString XName;
/*
	if(MultOpt == 0) XName = "N_{ch}^{FVTXN}/<N_{ch}^{FVTXN}>";
	if(MultOpt == 1) XName = "N_{ch}^{FVTXS}/<N_{ch}^{FVTXS}>";
	if(MultOpt == 2) XName = "N_{ch}^{SVX}/<N_{ch}^{SVX}>";
*/


	if(MultOpt == 0) XName = "#tilde{N}_{ch}^{N}/<#tilde{N}_{ch}^{N}>";
	if(MultOpt == 1) XName = "#tilde{N}_{ch}^{S}/<#tilde{N}_{ch}^{S}>";
	if(MultOpt == 2) XName = "#tilde{N}_{ch}^{M}/<#tilde{N}_{ch}^{M}>";


	TString MBHisName;
	if(Opt == 0) MBHisName = "Evt_Mult_FVTXN_His";
	if(Opt == 1) MBHisName = "Evt_Mult_FVTXS_His";
	if(Opt == 2) MBHisName = "Evt_Mult_SVX_His";


	gStyle->SetOptStat(0);
	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();
	c->SetBottomMargin(0.12);

	TFile * finMU = new TFile(Form("/sphenix/user/zshi/PHENIX/JPsiAnaCodes/FinalResults/Infiles/FitResults_%d_%d.root",Opt, MultOpt));
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

	//c->SaveAs(Form("FinalPlots/RawYields/JPsiYield_%d_%d.png",Opt,MultOpt));

	TFile * finMB = new TFile(Form("/sphenix/user/zshi/PHENIX/JPsiAnaCodes/FinalResults/Infiles/%sMBEvtMult.root",FileName.Data()));
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

	lat->DrawLatex(0.58,0.65,Form("%s",FVTXName.Data()));	

//	c2->SaveAs(Form("FinalPlots/MBEvent/MBHisName_%d_%d.png",Opt,MultOpt));



	TH1D * Ratio = (TH1D *) JPsiYield->Clone("Ratio");

	Ratio->Sumw2();
	Evt_Mult_FVTXN_His->Sumw2();

	Ratio->Divide(Evt_Mult_FVTXN_His);

	Ratio->Scale(1/AveJPsi);



//	Ratio->SetTitle(TitleName.Data());
	Ratio->GetYaxis()->SetTitle("<R>");
	Ratio->GetYaxis()->SetTitleOffset(1.4);
	Ratio->GetXaxis()->SetTitle(XName.Data());
	Ratio->GetYaxis()->CenterTitle();
	Ratio->GetXaxis()->CenterTitle();


	Ratio->SetMarkerStyle(20);
	Ratio->SetMarkerSize(1);
	Ratio->SetMarkerColor(1);
	Ratio->SetLineColor(1);

	const int NBins = 10;
	double MultBinScale[NBins + 1] = {0,1/AverageMB,2/AverageMB,3/AverageMB,4/AverageMB,5/AverageMB,6/AverageMB,8/AverageMB,10/AverageMB,12/AverageMB,19/AverageMB};
	double MultBin[NBins + 1] = {0,1,2,3,4,5,6,8,10,12,19};

	double MultBinShift[NBins + 1];


	float Shift = 1.4;

	for(int i = 0; i < NBins+ 1; i++){

		MultBinShift[i] = (MultBin[i] - Shift)/AverageMB;

	}
	
	

	float Value;
	float ValueErr;



	TString FuncForm; 


	if(MultOpt == 0) FuncForm = "0.906048 - 0.433618  * TMath::Exp(-0.283628  *x)";
	if(MultOpt == 1) FuncForm = "0.866227 - 0.386111  * TMath::Exp(-0.286963  *x)";
	if(MultOpt == 2) FuncForm = "0.985577 - 0.412983  * TMath::Exp(-0.196578  *x)";



	TF1 * func = new TF1("Func","x",0,19);
	func->SetLineStyle(2);
	func->Draw("SAME");

	TF1 * TrigBias = new TF1("TrigBias",FuncForm.Data(),0,19);




	//Add J/Psi Trig Bias Here//


	TString MUFuncForm; 


	if(MultOpt == 0) MUFuncForm = "1.099998 - 0.332713  * TMath::Exp(-0.069429  *x)";
	if(MultOpt == 1) MUFuncForm = "0.976232 - 0.222794  * TMath::Exp(-0.170423  *x)";
	if(MultOpt == 2) MUFuncForm = "0.864865 - 0.158782  * TMath::Exp(-0.434453  *x)";


	TF1 * MUTrigBias = new TF1("MUTrigBias",MUFuncForm.Data(),0,19);



	TH1D * RatioScaledCorr = new TH1D("RatioScaledCorr","",NBins,MultBinScale);
//	RatioScaledCorr->SetTitle(TitleName.Data());
//	RatioScaledCorr->GetYaxis()->SetTitle("<R> With Trig Bias Corrected");
	RatioScaledCorr->GetYaxis()->SetTitle("R (J/#psi)");	
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
	float MUCorrFactor;



	float RFinal;
	float RErrFinal;

	float BinCenter;



	TF1 * MultiFunc = new TF1("MultiFunc","9.83665e-01 + 6.44703e-03 * x + 3.02238e-04 * x * x + 4.83052e-06 * x * x * x * x",0,19);

	
	//Box for Systematics

	double MultBinCenter[NBins] = {0.5/AverageMB,1.5/AverageMB,2.5/AverageMB,3.5/AverageMB,4.5/AverageMB,5.5/AverageMB,7/AverageMB,9/AverageMB,11/AverageMB,15.5/AverageMB};
	double MultBinHigh[NBins] = {0.5/AverageMB,0.5/AverageMB,0.5/AverageMB,0.5/AverageMB,0.5/AverageMB,0.5/AverageMB,1/AverageMB,1/AverageMB,1/AverageMB,3.5/AverageMB};
	double MultBinLow[NBins] = {0.5/AverageMB,0.5/AverageMB,0.5/AverageMB,0.5/AverageMB,0.5/AverageMB,0.5/AverageMB,1/AverageMB,1/AverageMB,1/AverageMB,3.5/AverageMB};

	double MultBinCenterShifted[NBins];

	for(int i = 0; i < NBins; i++){

		MultBinCenterShifted[i] = MultBinCenter[i] - 1.4/AverageMB;

	}


	double YCenter[NBins];
	double YUp[NBins];
	double YDown[NBins];

	float AveJPsiEff = 0.79;
	float AveMBEff = 0.55;

	float scalefactor = AveJPsiEff/AveMBEff;

//	float scalefactor = 0.92/0.72;


	//Import Systematics//

	TFile * fin = new TFile("/sphenix/user/zshi/PHENIX/JPsiAnaCodes/SystStudies/SystTotalOut.root");
	fin->cd();

	TH1D * SystDownHis = (TH1D *)  fin->Get(Form("FinalErrorDown_%d_%d",Opt,MultOpt));
	TH1D * SystUpHis = (TH1D *)  fin->Get(Form("FinalErrorUp_%d_%d",Opt,MultOpt));
	
	float SystUpValue[NBins];
	float SystDownValue[NBins];

	for(int i = 0; i < NBins; i++){

		SystUpValue[i] = SystUpHis->GetBinContent(i+1)/100;
		SystDownValue[i] = -1 *SystDownHis->GetBinContent(i+1)/100;
		
		cout << "i = " << i << "   SystUpValue[i] = " << SystUpValue[i] << "    SystDownValue[i]  = " << SystDownValue[i]  << endl;
	}




	for(int i = 0; i < NBins; i++){


		R = Ratio->GetBinContent(i+1);
		RErr = Ratio->GetBinError(i+1);
		BinCenter = (MultBin[i] + MultBin[i+1]) * 0.5;

		CorrFactor = TrigBias->Eval(BinCenter);
		MultFactor = MultiFunc->Eval(BinCenter);
		MUCorrFactor = MUTrigBias->Eval(BinCenter);  //New Corr Factor



		RFinal = R * CorrFactor * MultFactor/MUCorrFactor * scalefactor;
		RErrFinal = RErr * CorrFactor * MultFactor/MUCorrFactor * scalefactor;


		RatioScaledCorr->SetBinContent(i+1,RFinal);
		RatioScaledCorr->SetBinError(i+1,RErrFinal);

		YCenter[i] = RFinal;
		YUp[i] = RFinal * SystUpValue[i];
		YDown[i] = RFinal * SystDownValue[i];



	}


	TGraphAsymmErrors *Syst = new TGraphAsymmErrors(NBins, MultBinCenter, YCenter,
			MultBinLow, MultBinHigh,
			YDown,YUp);

	
	if(Opt == 0){
				Syst->SetFillColorAlpha(kGreen-9,0.5);
				Syst->SetLineColor(kGreen-9);
	
	}
	
	if(Opt == 1){
				Syst->SetFillColorAlpha(kBlue-9,0.5);
				Syst->SetLineColor(kBlue-9);
	
	}

	c->cd();



	RatioScaledCorr->SetMarkerColor(Opt+3);
	RatioScaledCorr->SetLineColor(Opt+3);

	RatioScaledCorr->GetYaxis()->SetTitleSize(0.05);
	RatioScaledCorr->GetXaxis()->SetTitleSize(0.05);
	RatioScaledCorr->GetXaxis()->SetTitleOffset(0.90);
	RatioScaledCorr->GetYaxis()->SetTitleOffset(0.9);
	

	RatioScaledCorr->SetMaximum(20);
	RatioScaledCorr->SetMinimum(0);
	
	

	RatioScaledCorr->Draw("ep");
	

	Syst->Draw("5SAME");

	TLegend* leg = new TLegend(0.12,0.65,0.50,0.80,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.034);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->SetLineWidth(3);
//	if(Opt != MultOpt) leg->AddEntry(RatioScaledCorr,DataName[Opt].Data(),"PL");
//	if(Opt == MultOpt) leg->AddEntry(RatioScaledCorr,Form("%s %s",DataName[Opt].Data(),Together.Data()),"PL");
	leg->AddEntry(RatioScaledCorr,DataName[Opt].Data(),"PL");

	leg->AddEntry(Syst,"Systematics","f");	
	leg->Draw("same");


	lat->DrawLatex(0.15,0.61,"pp #sqrt{s} = 200 GeV");
	if(MultOpt == 0) lat->DrawLatex(0.15,0.55,"#tilde{N}_{ch}^{N}: 1.2 < #eta < 2.4, J/#psi #rightarrow #mu^{+} #mu^{-} subtracted");
	if(MultOpt == 1) lat->DrawLatex(0.15,0.55,"#tilde{N}_{ch}^{S}: -2.4 < #eta < -1.2, J/#psi #rightarrow #mu^{+} #mu^{-} subtracted");
	if(MultOpt == 2) lat->DrawLatex(0.15,0.55,"#tilde{N}_{ch}^{M}: -1.0 < #eta < 1.0, J/#psi #rightarrow #mu^{+} #mu^{-} subtracted");


	c->SaveAs(Form("Plots/Single/png/JPsiRatio_%d_%d.png",Opt,MultOpt));
	c->SaveAs(Form("Plots/Single/pdf/JPsiRatio_%d_%d.pdf",Opt,MultOpt));





	TString OutFileName = Form("OutFiles/JPsiR_%d_%d.root",Opt,MultOpt);

	TFile * fout = new TFile(OutFileName.Data(),"RECREATE");
	fout->cd();


	RatioScaledCorr->Write();
	Syst->Write();
	fout->Close();

}
