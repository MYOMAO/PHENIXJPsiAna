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




//#include "his.h"
using namespace std;
using std::cout;
using std::endl;

void PrintOutInc(){

	const int NBins = 1;
	int MultiBin[NBins + 1] = {0,19};

	const int NDir = 2;
	const int NOpt = 3;


	TString infile;

	TString Slash = "\\";

	TString MuonDir[NDir] = {"FVTXN","FVTXS"};
	TString MultiName[NOpt] = {"FVTXN","FVTXS","SVX"};


	for(int i = 0; i < NDir; i++){



		infile = Form("OutFiles/Inc/FitResults_%d.root",i);

		TFile * fin = new TFile(infile.Data());

		TH1D * JPsiYield = (TH1D *) fin->Get("JPsiYield");

		TH1D * JPsiMean = (TH1D *) fin->Get("JPsiMean");

		TH1D * JPsiWidth = (TH1D *) fin->Get("JPsiWidth");
		std::ofstream OutFile(Form("OutFiles/Text/%s.txt",MuonDir[i].Data()));


		for(int q = 0; q < NBins; q++){


			float Yield = JPsiYield->GetBinContent(q+1);
			float YieldErr = JPsiYield->GetBinError(q+1);

			float Mean = JPsiMean->GetBinContent(q+1);
			float MeanErr = JPsiMean->GetBinError(q+1);

			float Width = JPsiWidth->GetBinContent(q+1);
			float WidthErr = JPsiWidth->GetBinError(q+1);


			OutFile << Form("%d - %d & %.1f $\\pm$ %.1f & %.3f $\\pm$ %.3f & %.3f $\\pm$ %.3f  \\\\",MultiBin[q], MultiBin[q+1], Yield, YieldErr, Mean, MeanErr, Width, WidthErr)  << endl;

		}

	}

}



