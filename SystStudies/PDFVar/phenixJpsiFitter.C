#include "CMS_lumi.C"
#include "phenixJpsiFitter.h"

float MassRange = 0.15;
float JPsiMass = 3.096916;

void phenixJpsiFitter(int Opt, int MultOpt, int SigVar, int BackVar) {

	TString FolderName;

	if(SigVar == 0 && BackVar == 0){

		FolderName = "Nominal";

	}


	if(SigVar == 1){

		FolderName = "SingleCB";

	}

	if(BackVar == 1){

		FolderName = "Linear";

	}




	TString inputdata;
	if (Opt == 0)
		inputdata = "CommonFiles/NorthMU.root";
	if (Opt == 1)
		inputdata = "CommonFiles/SouthMU.root";

	gStyle->SetTextSize(0.05);
	gStyle->SetTextFont(42);
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadTopMargin(0.15);
	gStyle->SetPadBottomMargin(0.50);
	gStyle->SetTitleX(.0f);

	TFile *inf = new TFile(inputdata.Data());
	TTree *skimtree_new = (TTree *)inf->Get("anaT2");
	TH1F *h;
	h = new TH1F(Form("h%d", _count), "", nbinsmasshisto, minhisto, maxhisto);

	RooRealVar *mass =
		new RooRealVar("dimuon_mass", "dimuon_mass", minhisto, maxhisto);
	RooRealVar *y = new RooRealVar("dimuon_rapidity", "dimuon_rapidity", -3, 3);
	RooRealVar *Evt_Mult_FVTXN =
		new RooRealVar("Evt_Mult_FVTXN", "Evt_Mult_FVTXN", 0, 20);
	RooRealVar *Evt_Mult_FVTXS =
		new RooRealVar("Evt_Mult_FVTXS", "Evt_Mult_FVTXS", 0, 20);
	RooRealVar *Evt_Mult_SVX =
		new RooRealVar("Evt_Mult_SVX", "Evt_Mult_SVX", 0, 20);
	RooRealVar *MuID_N2D_Scale_Down =
		new RooRealVar("MuID_N2D_Scale_Down", "MuID_N2D_Scale_Down", -10, 10);
	RooRealVar *MuID_S2D_Scale_Down =
		new RooRealVar("MuID_S2D_Scale_Down", "MuID_S2D_Scale_Down", -10, 10);

	RooPlot *frame = new RooPlot();
	RooHist *datahist = new RooHist();
	RooCurve *modelcurve = new RooCurve();
	RooDataSet *ds = new RooDataSet();
	RooDataHist *dh = new RooDataHist("dh", "", *mass, Import(*h));


	TFile* outf = new TFile(Form("OutFiles/%s/FitResults_%d_%d.root",FolderName.Data(),Opt,MultOpt),"recreate");

	outf->cd();
	ds = new RooDataSet(
			"ds0", "", skimtree_new,
			RooArgSet(*mass, *y, *Evt_Mult_FVTXN, *Evt_Mult_FVTXS, *Evt_Mult_SVX));

	const int NBins = 10;
	int MultBin[NBins + 1] = {0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 19};
	float f_MultBin[NBins + 1] = {0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 19};


	TH1F *JPsiYield = new TH1F("JPsiYield", "", NBins, f_MultBin);
	TH1F *JPsiMean = new TH1F("JPsiMean", "", NBins, f_MultBin);
	TH1F *JPsiWidth = new TH1F("JPsiWidth", "", NBins, f_MultBin);


	for (int i = 0; i < NBins; i++) {
		TString Cut;
		w_val = new RooWorkspace("w_val");

		if (Opt == 0 && MultOpt == 0)
			Cut = Form(
					"dimuon_rapidity > 0 && Evt_Mult_FVTXN >= %d && Evt_Mult_FVTXN < %d",
					MultBin[i], MultBin[i + 1]);
		if (Opt == 0 && MultOpt == 1)
			Cut = Form(
					"dimuon_rapidity > 0 && Evt_Mult_FVTXS >= %d && Evt_Mult_FVTXS < %d",
					MultBin[i], MultBin[i + 1]);
		if (Opt == 0 && MultOpt == 2)
			Cut =
				Form("dimuon_rapidity > 0 && Evt_Mult_SVX >= %d && Evt_Mult_SVX < %d",
						MultBin[i], MultBin[i + 1]);
		if (Opt == 1 && MultOpt == 0)
			Cut = Form(
					"dimuon_rapidity < 0 && Evt_Mult_FVTXN >= %d && Evt_Mult_FVTXN < %d",
					MultBin[i], MultBin[i + 1]);
		if (Opt == 1 && MultOpt == 1)
			Cut = Form(
					"dimuon_rapidity < 0 && Evt_Mult_FVTXS >= %d && Evt_Mult_FVTXS < %d",
					MultBin[i], MultBin[i + 1]);
		if (Opt == 1 && MultOpt == 2)
			Cut =
				Form("dimuon_rapidity < 0 && Evt_Mult_SVX >= %d && Evt_Mult_SVX < %d",
						MultBin[i], MultBin[i + 1]);


		TString PreScale;
		if (Opt == 0)
			PreScale = "MuID_N2D_Scale_Down";
		if (Opt == 1)
			PreScale = "MuID_S2D_Scale_Down";

		RooDataSet *ds_cut;
		ds_cut = new RooDataSet(
				"ds", "", ds,
				RooArgSet(*mass, *y, *Evt_Mult_FVTXN, *Evt_Mult_FVTXS, *Evt_Mult_SVX),
				Cut.Data());

		TCanvas *c = new TCanvas(Form("c%d", _count), "", 600, 600);
		outputw->import(*ds);

		outputw->import(*dh);

		RooFitResult *f = fit("", "", c, ds_cut, dh, mass, frame,SigVar, BackVar, Opt);

		bool runToyStudy = false;
		if (runToyStudy) {
			validate_fit(w_val, Opt, MultOpt);
			modelcurve = frame->getCurve(Form("model%d", _count));

			RooRealVar *fitYield = static_cast<RooRealVar *>(f->floatParsFinal().at(
						f->floatParsFinal().index(Form("nsig%d", _count))));
		}

		CMS_lumi(c, 19011, 0);
		c->Update();

		TString Name;
		//		if(Opt == 0)  Name = "South + North"; //Never combine north and
		//south FVTX muons
		if (Opt == 0)
			Name = "North";
		if (Opt == 1)
			Name = "South";
		TLatex *tex_y;
		tex_y = new TLatex(0.21, 0.73, Name);
		tex_y->SetNDC();
		tex_y->SetTextFont(42);
		tex_y->SetTextSize(0.045);
		tex_y->SetLineWidth(2);
		tex_y->Draw();
		TString OutName;
		//	if(Opt == 0)  OutName = "North_AND_South"; //Never combine north and
		// south FVTX muons
		if (Opt == 0)
			OutName = "North";
		if (Opt == 1)
			OutName = "South";
		TString MultOutName;
		if (MultOpt == 0)
			MultOutName = "FVTXN";
		if (MultOpt == 1)
			MultOutName = "FVTXS";
		if (MultOpt == 2)
			MultOutName = "SVX";
		TString MultRange;
		MultRange = Form("%d <= N^{%s}_{TRK} <= %d", MultBin[i], MultOutName.Data(),
				MultBin[i + 1] - 1);

		TLatex *tex_Mult;
		tex_Mult = new TLatex(0.21, 0.68, MultRange.Data());
		tex_Mult->SetNDC();
		tex_Mult->SetTextFont(42);
		tex_Mult->SetTextSize(0.045);
		tex_Mult->SetLineWidth(2);
		tex_Mult->Draw();
		c->SaveAs(Form("Plots/FitResults/JPsiFit_%s_%s_%d.png", OutName.Data(),
					MultOutName.Data(), _count));
		c->SaveAs(Form("Plots/FitResults/JPsiFit_%s_%s_%d.pdf", OutName.Data(),
					MultOutName.Data(), _count));

		_count++;
		JPsiYield->SetBinContent(i + 1, yield);
		JPsiYield->SetBinError(i + 1, yieldErr);
		JPsiMean->SetBinContent(i + 1, FitMean);
		JPsiMean->SetBinError(i + 1, FitMeanErr);
		JPsiWidth->SetBinContent(i + 1, FitWidth1);
		JPsiWidth->SetBinError(i + 1, FitWidth1Err);

		cout << "yield OUTSIDE = " << yield << endl;
	}

	JPsiYield->Write();
	JPsiMean->Write();
	JPsiWidth->Write();
	// outf->Write();
	outf->Close();
}
