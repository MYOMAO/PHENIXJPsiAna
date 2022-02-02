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
#include "Fit/Fitter.h"





//#include "his.h"
using namespace std;

using std::cout;
using std::endl;


void PreFilterPlot(){

	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);


	const int NPtBins = 4;
	const int NCentBins = 3;

	double PtMinCut[NPtBins+1] = {7,10,15,20,10};
	double PtMaxCut[NPtBins+1] = {10,15,20,50,50};

	double CentMinCut[NCentBins+1] = {0,0,0,60};
	double CentMaxCut[NCentBins+1] = {180,180,60,180};



	TString PreCut = "((hiBin < 181) && Btrk1Pt > 1.0 && Btrk2Pt > 1.0 && Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2.6 && Bpt > 5 && abs(Btrk1Eta-0.0) < 2.4  && abs(Btrk2Eta-0.0) < 2.4 && (TMath::Abs(By)<2.4&&TMath::Abs(Bmumumass-3.096916)<0.15&&((abs(Bmu1eta)<1.2&&Bmu1pt>3.5)||(abs(Bmu1eta)>1.2&&abs(Bmu1eta)<2.1&&Bmu1pt>(5.47-1.89*abs(Bmu1eta)))||(abs(Bmu1eta)>2.1&&abs(Bmu1eta)<2.4&&Bmu1pt>1.5))&&((abs(Bmu2eta)<1.2&&Bmu2pt>3.5)||(abs(Bmu2eta)>1.2&&abs(Bmu2eta)<2.1&&Bmu2pt>(5.47-1.89*abs(Bmu2eta)))||(abs(Bmu2eta)>2.1&&abs(Bmu2eta)<2.4&&Bmu2pt>1.5))&&Bmu1TMOneStationTight&&Bmu2TMOneStationTight&&Bmu1InPixelLayer>0&&(Bmu1InPixelLayer+Bmu1InStripLayer)>5&&Bmu2InPixelLayer>0&&(Bmu2InPixelLayer+Bmu2InStripLayer)>5&&Bmu1dxyPV<0.3&&Bmu2dxyPV<0.3&&Bmu1dzPV<20&&Bmu2dzPV<20&&Bmu1isTrackerMuon&&Bmu2isTrackerMuon&&Bmu1isGlobalMuon&&Bmu2isGlobalMuon&&Btrk1highPurity&&Btrk2highPurity&&abs(Btrk1Eta)<2.4&&abs(Btrk2Eta)<2.4&&Btrk1Pt>1.&&Btrk2Pt>1.&&abs(Btktkmass-1.019455)<0.015) && (abs(PVz)<15&&pclusterCompatibilityFilter&&pprimaryVertexFilter) && (Btrk1PixelHit + Btrk1StripHit > 10) && (Btrk2PixelHit + Btrk2StripHit > 10) && (Btrk1PtErr/Btrk1Pt < 0.1)&& (Btrk2PtErr/Btrk2Pt < 0.1) && Btrk1Chi2ndf/(Btrk1nStripLayer+Btrk1nPixelLayer) < 0.18 && Btrk2Chi2ndf/(Btrk2nStripLayer+Btrk2nPixelLayer) < 0.18 && (HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1) && (abs(PVz)<15&&pclusterCompatibilityFilter&&pprimaryVertexFilter && phfCoincFilter2Th4))";
	

	TString Cut;
	TFile * fin = new TFile("/export/d00/scratch/zzshi/CMSSW_7_5_8_patch3/Merge/2018Ana/Samples/MUONJSONBs/ntphi_20190808_Bfinder_20190513_HIDoubleMuon__PsiPeri__HIRun2018A_04Apr2019_v1_HF_and_MuonJSON.root");


	TTree * nt = (TTree *) fin->Get("Bfinder/ntphi");

	nt->AddFriend("hltanalysis/HltTree");
	nt->AddFriend("hiEvtAnalyzer/HiTree");
	nt->AddFriend("skimanalysis/HltTree");



	TH1D * HisMass = new TH1D("HisMass","",50,5.0,6.0);
	HisMass->GetXaxis()->SetTitle("B_{s} Mass (GeV/c^{2})");
	HisMass->GetYaxis()->SetTitle("Counts");

	HisMass->SetMarkerStyle(20);
	HisMass->SetMarkerSize(1);


	TCanvas *c = new TCanvas("c","c",600,600);
	c->cd();

	for(int i = 0; i < NPtBins; i++){

		Cut = Form("%s &&(Bpt > %f && Bpt < %f)",PreCut.Data(),PtMinCut[i],PtMaxCut[i]);

		nt->Project("HisMass","Bmass",Cut.Data());
		
		HisMass->Draw("ep");

		c->SaveAs(Form("PreCutPlots/Mass_Pt_%d.png",i));

	}

	for(int i = 0; i < NCentBins; i++){

		Cut = Form("%s &&(hiBin > %f && hiBin < %f) && (Bpt > 10 && Bpt < 50)",PreCut.Data(),CentMinCut[i],CentMaxCut[i]);

		nt->Project("HisMass","Bmass",Cut.Data());
		
		HisMass->Draw("ep");

		c->SaveAs(Form("PreCutPlots/Mass_Cent_%d.png",i));


	}



}
