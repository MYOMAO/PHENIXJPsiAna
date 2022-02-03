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

void Skim(){


	TFile * fin = new TFile("merged_anaT2_run15pp_MU_v14.root");

	TTree * anaT2 = (TTree *) fin->Get("anaT2");

	TFile *fout = new TFile("Skimmed.root","RECREATE");
	fout->cd();
	TTree * anaT2_new = anaT2->CloneTree(0);
	
	int NEvents = 200;

	for(int i = 0; i < NEvents; i++){

		anaT2->GetEntry(i);
		anaT2_new->Fill();

	
	}

	anaT2_new->Write();
	fout->Close();
}
