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



void DrawReference(){



	gStyle->SetOptStat(0);
	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();



	const int NTypes = 3;
	const int NDir = 2;

	TString infile;

	TString OutFileName[NTypes]={"FVTXN","FVTXS","SVX"};
	TString DataName[NDir]={"FVTX North J/#psi","FVTX South J/#psi"};

	TLegend* leg[NTypes];

	float Min[NTypes] = {0,0,0};
	float Max[NTypes] = {18,18,18};


	TF1 * func = new TF1("Func","x",0,19);
	func->SetLineStyle(2);

	//STAR Reference//

	const int NPointSTAR = 4;
	float XSTAR[NPointSTAR] = {0.4553571428571437,1.321428571428572,2.1875,3.0982142857142865};	
	float XSTARErr[NPointSTAR] = {0,0,0,0};	
	float XSTARSystErr[NPointSTAR] = {0.1,0.1,0.1,0.1};	

	float YSTAR[NPointSTAR] = {0.4210526315789451,1.6842105263157876,3.0526315789473664,5.210526315789473};
	float YSTARStatLow[NPointSTAR] = {0.4210526315789451,0.9999999999999964,2.6842105263157876,4.263157894736839};
	float YSTARSystLow[NPointSTAR] = {0.4210526315789451,1.3157894736842088,2.7894736842105274,3.3684210526315788};

	float StatErrorSTAR[NPointSTAR];
	float SystErrorSTAR[NPointSTAR];

	for(int i = 0; i < NPointSTAR; i++){

		StatErrorSTAR[i] = YSTAR[i] - YSTARStatLow[i];
		SystErrorSTAR[i] = YSTAR[i] - YSTARSystLow[i];

	}

	/*
	//ALICE OLD Reference//

	const int NPointALICE = 5;
	float XALICE[NPointALICE] = {0.4553571428571437,1.1785714285714293,1.7767857142857149,2.6250000000000004,4};	
	float XALICEErr[NPointALICE] = {0,0,0,0,0};	
	float YALICE[NPointALICE] = {0.4210526315789451,1.0526315789473664,2.315789473684209,3.1052631578947327,7.94736842105263};
	float YALICEStatLow[NPointALICE] = {0.4210526315789451,1.0526315789473664,2.315789473684209,3.1052631578947327,6.473684210526315};
	float YALICESystLow[NPointALICE] = {0.4210526315789451,0.8947368421052637,2.0526315789473664,2.7894736842105274,7.7368421052631575};
	float XALICESystError[NPointALICE] = {0.1,0.1,0.1,0.1,0.1};	

	float StatErrorALICE[NPointALICE];
	float SystErrorALICE[NPointALICE];

	for(int i = 0; i < NPointALICE; i++){

	StatErrorALICE[i] = YALICE[i] - YALICEStatLow[i];
	SystErrorALICE[i] = YALICE[i] - YALICESystLow[i];

	}
	*/
	//ALICE New Reference//


	const int NPointALICE = 10;
	float XALICE[NPointALICE];
	float XALICEErr[NPointALICE] = {0,0,0,0,0};	
	float YALICE[NPointALICE] = {0.081,	0.387,1.230,2.498,4.201,6.208, 7.974, 10.071, 13.274, 14.968};
	float YALICECorr[NPointALICE] = {0.355,	0.640,1.003,1.186,1.411,1.642, 1.740, 1.871, 2.153, 2.100};


	//	float YALICEStatLow[NPointALICE] = {0.4210526315789451,1.0526315789473664,2.315789473684209,3.1052631578947327,6.473684210526315};
	//	float YALICESystLow[NPointALICE] = {0.4210526315789451,0.8947368421052637,2.0526315789473664,2.7894736842105274,7.7368421052631575};
	float XALICESystError[NPointALICE] = {0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1};	

	float StatErrorALICE[NPointALICE] = {0.011,0.029,0.054,0.120,0.182,0.229,0.314,0.480,0.836,1.420};
	float SystErrorALICE[NPointALICE] = {0.004,0.019,0.035,0.080,0.143,0.169,0.216,0.506,0.794,1.056};

	for(int i = 0; i < NPointALICE; i++){
		XALICE[i] = YALICE[i]/YALICECorr[i];

		//		StatErrorALICE[i] = YALICE[i] - YALICEStatLow[i];
		//		SystErrorALICE[i] = YALICE[i] - YALICESystLow[i];

	}



	const int NPointALICE2 = 9;
	float XALICE2[NPointALICE2] = {0.37777171317459235,1.0115985535921403,2.025243910759364,2.557412840281095,3.229037320051852,3.913283755202292,4.6734666029883325,5.674421777990039,7.270450979054377};
	float XALICEErr2[NPointALICE2] = {0,0,0,0,0};	
	float YALICE2[NPointALICE2] = {0.32673807736917837,	0.8711878283414087 ,2.1186463805690146,2.7695299174455883,3.5294398580882884,4.343248959541519, 5.31916490414137, 6.566487002797299, 8.89547656409907};



	float YALICEStatLow2[NPointALICE2] = {0.32673807736917837,	0.8711878283414087 ,2.1186463805690146,2.7695299174455883,3.5294398580882884,4.343248959541519, 5.31916490414137, 6.566487002797299, 9.43310363648769};
	float YALICESystLow2[NPointALICE2] = {0.32673807736917837,	0.8711878283414087 ,2.1186463805690146,2.7695299174455883,3.5294398580882884,4.612062495735827, 5.587978440335679, 6.8889267926588005, 9.48686634372655};

	float XALICESystError2[NPointALICE2] = {0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1};	

	float StatErrorALICE2[NPointALICE2];
	float SystErrorALICE2[NPointALICE2];

	for(int i = 0; i < NPointALICE2; i++){

		StatErrorALICE2[i] = YALICEStatLow2[i] - YALICE2[i];
		SystErrorALICE2[i] = YALICESystLow2[i] - YALICE2[i];

	}



	TGraphErrors * STAR = new TGraphErrors(NPointSTAR,XSTAR,YSTAR,XSTARErr,StatErrorSTAR);
	STAR->SetLineColor(kOrange-3);
	STAR->SetMarkerColor(kOrange-3);
	STAR->SetMarkerSize(1);
	STAR->SetMarkerStyle(21);

	TGraphErrors * STARSyst = new TGraphErrors(NPointSTAR,XSTAR,YSTAR,XSTARSystErr,SystErrorSTAR);
	STARSyst->SetFillColorAlpha(kOrange-3,0.5);
	STARSyst->SetLineColor(kOrange-3);



	TGraphErrors * ALICEMidRap = new TGraphErrors(NPointALICE,XALICE,YALICE,XALICEErr,StatErrorALICE);
	ALICEMidRap->SetLineColor(kMagenta-3);
	ALICEMidRap->SetMarkerColor(kMagenta-3);
	ALICEMidRap->SetMarkerSize(1);
	ALICEMidRap->SetMarkerStyle(33);

	TGraphErrors * ALICEMidRapSyst = new TGraphErrors(NPointALICE,XALICE,YALICE,XALICESystError,SystErrorALICE);
	ALICEMidRapSyst->SetFillColorAlpha(kMagenta-3,0.5);
	ALICEMidRapSyst->SetLineColor(kMagenta-3);





	TGraphErrors * ALICEPeri = new TGraphErrors(NPointALICE2,XALICE2,YALICE2,XALICEErr2,StatErrorALICE2);
	ALICEPeri->SetLineColor(kCyan+2);
	ALICEPeri->SetMarkerColor(kCyan+2);
	ALICEPeri->SetMarkerSize(1);
	ALICEPeri->SetMarkerStyle(34);

	TGraphErrors * ALICEPeriSyst = new TGraphErrors(NPointALICE2,XALICE2,YALICE2,XALICESystError2,SystErrorALICE2);
	ALICEPeriSyst->SetFillColorAlpha(kCyan+2,0.5);
	ALICEPeriSyst->SetLineColor(kCyan+2);




	for(int i = 0; i < NTypes; i ++){


		leg[i] = new TLegend(0.12,0.65,0.50,0.80,NULL,"brNDC");
		leg[i]->SetBorderSize(0);
		leg[i]->SetTextSize(0.035);
		leg[i]->SetTextFont(42);
		leg[i]->SetFillStyle(0);
		leg[i]->SetLineWidth(3);




		infile = Form("OutFiles/JPsiR_%d_%d.root",i,i);

		TFile * fin = new TFile(infile.Data());

		TH1D * RatioScaledCorr = (TH1D * ) fin->Get("RatioScaledCorr");
		TF1 * FitFunc = (TF1 * ) fin->Get("FitFunc");


		TGraphAsymmErrors *Syst = (TGraphAsymmErrors *) fin->Get("Graph");

		RatioScaledCorr->SetMarkerColor(j+3);
		RatioScaledCorr->SetLineColor(j+3);
		FitFunc->SetLineColor(j+3);

		if(i == 0){
			Syst->SetFillColorAlpha(kGreen-9,0.5);
			Syst->SetLineColor(kGreen-9);
		}
		if(i == 1){
			Syst->SetFillColorAlpha(kBlue-9,0.5);
			Syst->SetLineColor(kBlue-9);
		}


		RatioScaledCorr->SetMinimum(Min[i]);
		RatioScaledCorr->SetMaximum(Max[i]);				
		RatioScaledCorr->SetTitle("");
		RatioScaledCorr->Draw("ep");




		leg[i]->AddEntry(RatioScaledCorr,DataName[j].Data(),"LP");
		Syst->Draw("5SAME");


		STAR->Draw("epSAME");
		STARSyst->Draw("5SAME");

		ALICEMidRap->Draw("epSAME");
		ALICEMidRapSyst->Draw("5SAME");

		ALICEPeri->Draw("epSAME");
		ALICEPeriSyst->Draw("5SAME");

		leg[i]->AddEntry(STAR,"STAR: #sqrt{s} = 200 GeV, |y| < 1","LP");
		leg[i]->AddEntry(ALICEMidRap,"ALICE: #sqrt{s} = 13 TeV, |y| < 0.9","LP");
		leg[i]->AddEntry(ALICEPeri,"ALICE: #sqrt{s} = 13 TeV, 2.5 < y < 4","LP");


		leg[i]->Draw("SAME");

		c->SaveAs(Form("Plots/Reference/png/%s.png",OutFileName[i].Data()));
		c->SaveAs(Form("Plots/Reference/pdf/%s.pdf",OutFileName[i].Data()));

	}








}




