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

void CheckRunList(int Opt){



	TString infile;

	if(Opt == 0) infile = "CommonFiles/NorthMU.root";
	if(Opt == 1) infile = "CommonFiles/SouthMU.root";


	TFile * fin = new TFile(infile.Data());
	fin->cd();

	TTree * anaT2 = (TTree *) fin->Get("anaT2");

	int RunNumber;

	anaT2->SetBranchAddress("RunNumber",&RunNumber);

	int NEvents = anaT2->GetEntries();


	std::vector<int> DiffRuns;


	for(int i = 0; i < NEvents; i++){

		int Distinct = 1;

		anaT2->GetEntry(i);

		int Size = DiffRuns.size();

		if(Size == 0) DiffRuns.push_back(RunNumber);

		if(Size > 0){
			for(int j = 0; j < Size; j++){

				if(RunNumber == DiffRuns[j]) Distinct = 0;

			}
		}
	
		if(Distinct == 1) DiffRuns.push_back(RunNumber);

	}


	int FinalSize = DiffRuns.size();


	

	for(int i = 0; i < FinalSize; i++){


		cout  << DiffRuns[i] << endl;

	}




	cout << "FinalSize = " << FinalSize << endl;
}
