#include "roofitB.h"
#include "CMS_lumi.C"
//int _nBins = nBins;
//double *_ptBins = ptBins;



float MassRange = 0.15;
float JPsiMass = 3.096916;





void roofitB(int Opt,int MultOpt)
{


	cout << "Pass 1" << endl;
	//	TString inputdata = "merged_anaT2_run15pp_MU_v14.root";

	TString inputdata;

	if(Opt == 0) inputdata = "CommonFiles/NorthMU.root";
	if(Opt == 1) inputdata = "CommonFiles/SouthMU.root";

	//return;
	gStyle->SetTextSize(0.05);
	gStyle->SetTextFont(42);
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadTopMargin(0.15);
	gStyle->SetPadBottomMargin(0.50);

	gStyle->SetTitleX(.0f);

	TFile* inf = new TFile(inputdata.Data());
	TTree* skimtree_new = (TTree*)inf->Get("anaT2");

	TH1D* h;
	h = new TH1D(Form("h%d",_count),"",nbinsmasshisto,minhisto,maxhisto);

	cout << "Pass 2" << endl;

	RooRealVar* mass = new RooRealVar("dimuon_mass","dimuon_mass",minhisto,maxhisto);
	RooRealVar* y = new RooRealVar("dimuon_rapidity","dimuon_rapidity",-3,3);
	RooRealVar* Evt_Mult_FVTXN = new RooRealVar("Evt_Mult_FVTXN","Evt_Mult_FVTXN",0,20);
	RooRealVar* Evt_Mult_FVTXS = new RooRealVar("Evt_Mult_FVTXS","Evt_Mult_FVTXS",0,20);
	RooRealVar* Evt_Mult_SVX = new RooRealVar("Evt_Mult_SVX","Evt_Mult_SVX",0,20);

	RooRealVar* MuID_N2D_Scale_Down = new RooRealVar("MuID_N2D_Scale_Down","MuID_N2D_Scale_Down",-10,10);
	RooRealVar* MuID_S2D_Scale_Down = new RooRealVar("MuID_S2D_Scale_Down","MuID_S2D_Scale_Down",-10,10);


	RooPlot* frame = new RooPlot();
	RooHist* datahist = new RooHist();
	RooCurve* modelcurve = new RooCurve();
	RooDataSet* ds = new RooDataSet();


	RooDataHist* dh = new RooDataHist();   
	dh = new RooDataHist("dh","",*mass,Import(*h));

	TFile* outf = new TFile(Form("OutFiles/FitResults_%d_%d.root",Opt,MultOpt),"recreate");
	outf->cd();

	cout << "Pass 3" << endl;


	ds = new RooDataSet("ds0","",skimtree_new, RooArgSet(*mass, *y, *Evt_Mult_FVTXN, *Evt_Mult_FVTXS, *Evt_Mult_SVX));


	const int NBins = 5;
	int MultBin[NBins + 1] = {0,2,5,8,12,19};
	double MultBinHis[NBins + 1] = {0,2,5,8,12,19};


	TH1D * JPsiYield = new TH1D("JPsiYield","",NBins,MultBinHis);
	TH1D * JPsiMean = new TH1D("JPsiMean","",NBins,MultBinHis);
	TH1D * JPsiWidth1 = new TH1D("JPsiWidth1","",NBins,MultBinHis);
	TH1D * JPsiWidth2 = new TH1D("JPsiWidth2","",NBins,MultBinHis);




	for(int i=0;i<NBins;i++)
	{
		TString Cut;

		//		if(Opt==0) Cut = Form("dimuon_rapidity > -999 && Evt_Mult_FVTXN + Evt_Mult_FVTXS >= %d && Evt_Mult_FVTXN + Evt_Mult_FVTXS < %d",MultBin[i],MultBin[i+1]);
		if(Opt==0 && MultOpt == 0) Cut = Form("dimuon_rapidity > 0  && Evt_Mult_FVTXN >= %d && Evt_Mult_FVTXN < %d",MultBin[i],MultBin[i+1]);
		if(Opt==0 && MultOpt == 1) Cut = Form("dimuon_rapidity > 0  && Evt_Mult_FVTXS >= %d && Evt_Mult_FVTXS < %d",MultBin[i],MultBin[i+1]);
		if(Opt==0 && MultOpt == 2) Cut = Form("dimuon_rapidity > 0  && Evt_Mult_SVX >= %d && Evt_Mult_SVX < %d",MultBin[i],MultBin[i+1]);

		if(Opt==1 && MultOpt == 0) Cut = Form("dimuon_rapidity < 0  && Evt_Mult_FVTXN >= %d && Evt_Mult_FVTXN < %d",MultBin[i],MultBin[i+1]);		
		if(Opt==1 && MultOpt == 1) Cut = Form("dimuon_rapidity < 0  && Evt_Mult_FVTXS >= %d && Evt_Mult_FVTXS < %d",MultBin[i],MultBin[i+1]);
		if(Opt==1 && MultOpt == 2) Cut = Form("dimuon_rapidity < 0  && Evt_Mult_SVX >= %d && Evt_Mult_SVX < %d",MultBin[i],MultBin[i+1]);



		TString PreScale;

		if(Opt==0) PreScale = "MuID_N2D_Scale_Down";
		if(Opt==1) PreScale = "MuID_S2D_Scale_Down";
		RooDataSet* ds_cut;


		//ds_cut = new RooDataSet("ds","",ds, RooArgSet(*mass, *y, *Evt_Mult_FVTXN, *Evt_Mult_FVTXS), Cut.Data(),PreScale.Data()); 
		ds_cut = new RooDataSet("ds","",ds, RooArgSet(*mass, *y, *Evt_Mult_FVTXN, *Evt_Mult_FVTXS, *Evt_Mult_SVX), Cut.Data()); 

		cout << "Pass 3.2" << endl;


		//		TCanvas* c= new TCanvas("c","",600,600);

		TCanvas* c= new TCanvas(Form("c%d",_count),"",600,600);


		outputw->import(*ds);
		cout << "Pass 3.3" << endl;

		outputw->import(*dh);

		cout << "Pass 4" << endl;

		RooFitResult* f = fit("", "", c,  ds_cut, dh, mass, frame);

		std::cout << "Now Finall We Validate Our Fits" << std::endl;
		//	validate_fit(w_val);
		modelcurve = frame->getCurve(Form("model%d",_count));

		cout << "Got Curve" << endl;

		RooRealVar* fitYield = static_cast<RooRealVar*>(f->floatParsFinal().at(f->floatParsFinal().index(Form("nsig%d",_count))));
		double yield = fitYield->getVal();


		cout << "Pass Here 8" << endl;

		CMS_lumi(c,19011,0);
		c->Update();

		/*
		   TLatex* tex_pt;
		   TLatex* tex_hibin;

		   TLatex* tex_y;
		   tex_pt->SetNDC();
		   tex_pt->SetTextFont(42);
		   tex_pt->SetTextSize(0.045);
		   tex_pt->SetLineWidth(2);

		   tex_pt->Draw();


		   tex_hibin->SetNDC();
		   tex_hibin->SetTextFont(42);
		   tex_hibin->SetTextSize(0.045);
		   tex_hibin->SetLineWidth(2);
		   tex_hibin->Draw();

		   std::cout<<"CHEGUEI AQUI"<<std::endl;


		   tex_y->SetNDC();
		   tex_y->SetTextFont(42);
		   tex_y->SetTextSize(0.045);
		   tex_y->SetLineWidth(2);

		   tex_y->Draw();

		   std::cout<<"CHEGUEI AQUI"<<std::endl;
		   */

		TString Name;
		//		if(Opt == 0)   Name = "South + North";  //Never combine north and south FVTX muons
		if(Opt == 0)   Name = "North";
		if(Opt == 1)   Name = "South";


		TLatex* tex_y;
		tex_y = new TLatex(0.21,0.73,Name);


		tex_y->SetNDC();
		tex_y->SetTextFont(42);
		tex_y->SetTextSize(0.045);
		tex_y->SetLineWidth(2);

		tex_y->Draw();

		TString OutName;
		//	if(Opt == 0)   OutName = "North_AND_South"; //Never combine north and south FVTX muons
		if(Opt == 0)   OutName = "North";
		if(Opt == 1)   OutName = "South";


		TString MultOutName;
		if(MultOpt == 0)   MultOutName = "FVTXN";
		if(MultOpt == 1)   MultOutName = "FVTXS";
		if(MultOpt == 2)   MultOutName = "SVX";

		TString MultRange;

		MultRange = Form("%d <= N^{%s}_{TRK} <= %d",MultBin[i],MultOutName.Data(),MultBin[i+1]-1);


		TLatex* tex_Mult;
		tex_Mult = new TLatex(0.21,0.68,MultRange.Data());


		tex_Mult->SetNDC();
		tex_Mult->SetTextFont(42);
		tex_Mult->SetTextSize(0.045);
		tex_Mult->SetLineWidth(2);

		tex_Mult->Draw();


		c->SaveAs(Form("Plots/FitResults/JPsiFit_%s_%s_%d.png",OutName.Data(),MultOutName.Data(),_count));
		_count++;

		JPsiYield->SetBinContent(i+1,yield);
		JPsiYield->SetBinError(i+1,yieldErr);

		JPsiMean->SetBinContent(i+1,FitMean);
		JPsiMean->SetBinError(i+1,FitMeanErr);


		JPsiWidth1->SetBinContent(i+1,FitWidth1);
		JPsiWidth1->SetBinError(i+1,FitWidth1Err);


		JPsiWidth2->SetBinContent(i+1,FitWidth2);
		JPsiWidth2->SetBinError(i+1,FitWidth2Err);



	}

	JPsiYield->Write();
	JPsiMean->Write();
	JPsiWidth1->Write();	
	JPsiWidth2->Write();

	//outf->Write();
	outf->Close();

}
