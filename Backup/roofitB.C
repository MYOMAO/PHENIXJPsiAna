#include "roofitB.h"
#include "CMS_lumi.C"
int _nBins = nBins;
double *_ptBins = ptBins;

int syst=0;
int bp_bs = 0u;

TTree* makeTTree(TTree* intree, TString treeTitle) 
{
	TTree* outtree = new TTree(treeTitle.Data(),treeTitle.Data());
	float mass;
	outtree->Branch("mass",&mass,"mass/F") ;
	int Bsize;
	float BmassNew[20000];
	intree->SetBranchAddress("Bsize",&Bsize);
	intree->SetBranchAddress("BmassNew",BmassNew);
	int nentries = intree->GetEntries();
	for(int n=0; n<nentries; n++){
		intree->GetEntry(n);
		for(int b=0; b<Bsize; b++){
			mass = BmassNew[b];
			outtree->Fill();
		}
	}
	return outtree ;
}

//void roofitB(int usePbPb = 0, int fitOnSaved = 0, TString inputmc = "", TString _varExp = "", TString trgselection = "",  TString cut = "", TString cutmcgen = "", int isMC = 0, Double_t luminosity = 1., int doweight = 1, TString collsyst = "", TString outputfile = "", TString outplotf = "", TString npfit = "", int doDataCor = 0, Float_t centmin = 0., Float_t centmax = 100.)

void roofitB(int doubly = 0, TString tree = "ntphi", int full = 1, int usePbPb = 0, int fitOnSaved = 0, TString inputdata = "", TString inputmc = "", TString _varExp = "", TString trgselection = "",  TString cut = "", TString cutmcgen = "", int isMC = 0, Double_t luminosity = 1., int doweight = 1, TString collsyst = "", TString outputfile = "", TString outplotf = "", TString npfit = "", int doDataCor = 0, Float_t centmin = 0., Float_t centmax = 100.)
{
	std::cout<<"DEBUG "<<std::endl;
	double MyBackground;
	//	return;
	double	yieldRec;
	collisionsystem=collsyst;

	TString varExp = _varExp;
	if(_varExp == "Bpt750"){
		_nBins = nBins750;
		_ptBins = ptBins750;
		varExp = "BptNew";
	}
	/*	if(_varExp == "abs(By)"){
		_nBins = nBinsY;
		_ptBins = ptBinsY;
		}*/
	if(collisionsystem=="ppInc"||collisionsystem=="PbPbInc"){
		_nBins = nBinsInc;
		_ptBins = ptBinsInc;
	}

	hiBinMin = centmin*2;
	hiBinMax = centmax*2;
	centMin = centmin;
	centMax = centmax;
	std::cout<<"DEBUG 2 "<<std::endl;


	if (!(usePbPb==1||usePbPb==0)) std::cout<<"ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!, you are using a non valid isPbPb option"<<std::endl;
	bool isPbPb=(bool)(usePbPb);

	std::cout<<"DEBUG 3"<<std::endl;
	if(!isPbPb)
	{
		std::cout<<"DEBUG 4 !PbPb"<<std::endl;
		seldata = Form("%s",trgselection.Data());
		std::cout<<"DEBUG 5"<<std::endl;
		selmc = Form("%s",cut.Data());
		std::cout<<"DEBUG 6"<<std::endl;
		selmcgen = Form("%s",cutmcgen.Data());
		std::cout<<"DEBUG 7"<<std::endl;
	}
	else
	{
		std::cout<<"DEBUG 4 else"<<std::endl;
		seldata = Form("%s&&%s&&BsizeNew>0",trgselection.Data(),cut.Data());
		std::cout<<"DEBUG 5"<<std::endl;
		selmc = Form("%s&&BsizeNew>0",cut.Data());
		std::cout<<"DEBUG 6"<<std::endl;
		selmcgen = Form("%s&&BsizeNew>0",cutmcgen.Data());
		std::cout<<"DEBUG 7"<<std::endl;
		std::cout<<"cut = "<<cut<<std::endl;
		std::cout<<"seldata= "<<seldata<<std::endl;
		std::cout<<"selmc= "<<selmc<<std::endl;
	}
	//return;
	gStyle->SetTextSize(0.05);
	gStyle->SetTextFont(42);
	gStyle->SetPadRightMargin(cRightMargin);
	gStyle->SetPadLeftMargin(cLeftMargin);
	gStyle->SetPadTopMargin(cTopMargin);
	//gStyle->SetPadBottomMargin(cBottomMargin);
	cout << "cBottomMargin = " << cBottomMargin << endl;
	gStyle->SetPadBottomMargin(0.50);

	gStyle->SetTitleX(.0f);

	TFile* inf = new TFile(inputdata.Data());
	TTree* skimtree_new = (TTree*)inf->Get(tree);
	//TTree* skimtree_new =skimtree->CloneTree();
	//  inf->Close();
	TFile* infMC = new TFile(inputmc.Data());
	std::cout<<"MC file: "<<inputmc.Data()<<std::endl;
	//  infMC->cd();
	std::cout<<"cd"<<std::endl;
	TTree* skimtreeMC_new = (TTree*)infMC->Get(tree);
	//TTree* skimtreeMC_new =skimtreeMC->CloneTree();
	//  fout->cd();
	std::cout<<"depois do get"<<std::endl;
	std::cout<<"clonou"<<std::endl;
	//  inf->Close();
	//  fout->Close();

	TH1D* h;
	TH1D* hMC;
	TH1D* hpull;

	//TTree* ntGen = new TTree();
	/*TTree* skimtree= new TTree();
	  TTree* skimtreeMC = new TTree();
	  TTree* skimtree_new= new TTree();
	  TTree* skimtreeMC_new = new TTree();*/

	//if(fitOnSaved == 0){

	std::cout<<"Creating trees "<<std::endl;

	/*skimtree->SetObject("ntphi_data", "ntphi_data");  
	  skimtree_new->SetObject("ntphi_data", "ntphi_data");  */
	//nt->AddFriend("hltanalysis/HltTree");
	//nt->AddFriend("hiEvtAnalyzer/HiTree");
	//nt->AddFriend("skimanalysis/HltTree");
	/*
	   nt->AddFriend("BDT_pt_15_20");
	   nt->AddFriend("BDT_pt_7_15");
	   nt->AddFriend("BDT_pt_5_7");
	   nt->AddFriend("BDT_pt_20_50");	
	   */
	//nt->AddFriend("BDT");	

	/*
	   ntGen = (TTree*)infMC->Get("ntGen");
	   ntGen->AddFriend("ntHlt");
	   ntGen->AddFriend("ntHi");
	   cout << "Pass 2" << endl;
	   ntMC = (TTree*)infMC->Get("ntphi");
	   ntMC->AddFriend("ntHlt");
	   ntMC->AddFriend("ntHi");
	   ntMC->AddFriend("ntSkim");
	   ntMC->AddFriend("mvaTree");
	   ntMC->AddFriend("ntGen");
	   */
	/*	ntGen = (TTree*)infMC->Get("Bfinder/ntGen");
		ntGen->AddFriend("hltanalysis/HltTree");
		ntGen->AddFriend("hiEvtAnalyzer/HiTree");*/
	/*  infMC->cd();
		std::cout<<"cd"<<std::endl;
		skimtreeMC = (TTree*)infMC->Get("ntphi");
		std::cout<<"depois do get"<<std::endl;
		skimtreeMC_new =skimtreeMC->CloneTree();
		std::cout<<"clonou"<<std::endl;
		TFile* fout = new TFile("/afs/cern.ch/work/j/jusilva/CMSSW_7_5_8_patch5/src/UserCode/Bmeson2018Ana/Bs/CrossSection/test_trees/test_mc_new.root", "RECREATE");
		fout->cd();
		skimtreeMC->Write();
		fout->Close();*/
	/*		ntMC->AddFriend("hltanalysis/HltTree");
			ntMC->AddFriend("hiEvtAnalyzer/HiTree");
			ntMC->AddFriend("skimanalysis/HltTree");*/
	/*
	   ntMC->AddFriend("BDT_pt_15_20");
	   ntMC->AddFriend("BDT_pt_7_15");
	   ntMC->AddFriend("BDT_pt_5_7");
	   ntMC->AddFriend("BDT_pt_20_50");	
	   */
	//		ntMC->AddFriend("BDT");
	//	ntMC->AddFriend("Bfinder/ntGen");	
	//}
	std::cout<<"cd'ed"<<std::endl;
	int nBins_check = 0;
	double *ptBins_check;

	if(_varExp=="BptNew"){
		if(full==0) {
			nBins_check = nBins;
			ptBins_check = ptBins;
			if(tree=="ntKp" && bp_bs==0){
				nBins_check = nBins_bp;
				ptBins_check = ptBins_bp;

			}
		}
		else if(full==1){
			nBins_check = nBins_full;
			ptBins_check = ptBins_full;

		}

	}
	if(_varExp=="HiBin"){
		if(full==0) {
			nBins_check = nBinshi;
			ptBins_check = Binshi;

		}
		else if(full==1){
			nBins_check = nBins_full;
			ptBins_check = hiBins_full;

		}

	}
	else if(_varExp=="abs(By)"){
		if(full==0) {
			nBins_check = nBinsY;
			ptBins_check = ptBinsY;

		}
		else if(full==1){
			nBins_check = nBins_full;
			ptBins_check = hiBins_full;

		}

	}


	TH1D* hPt = new TH1D("hPt","",nBins_check,ptBins_check);


	TH1D* hPtMC = new TH1D("hPtMC","",_nBins,_ptBins);
	//TH1D* hPtGen = new TH1D("hPtGen","",_nBins,_ptBins);
	TH1D* hMean = new TH1D("hMean","",nBins_check,ptBins_check);                       
	TH1D* hSigmaGaus1 = new TH1D("hSigmaGaus1","",_nBins,_ptBins); 
	TH1D* hSigmaGaus2 = new TH1D("hSigmaGaus2","",_nBins,_ptBins); 
	std::cout<<"Histograms"<<std::endl;
	//RooWorkspace* w = new RooWorkspace("w");
	RooRealVar* mass = new RooRealVar("BmassNew","BmassNew",minhisto,maxhisto);
	RooRealVar* pt = new RooRealVar("BptNew","BptNew",5,50);
	RooRealVar* y = new RooRealVar("ByNew","ByNew",-2.4, 2.4);
//	RooRealVar* HiBin = new RooRealVar("hiBinNew","hiBinNew",0.,90.*2);
	RooRealVar* BDT_pt_5_10 = new RooRealVar("BDT_pt_5_10_New","BDT_pt_5_10_New",-1.0, 1.0);
	RooRealVar* BDT_pt_10_15 = new RooRealVar("BDT_pt_10_15_New","BDT_pt_10_15_New",-1.0, 1.0);
	RooRealVar* BDT_pt_15_20 = new RooRealVar("BDT_pt_15_20_New","BDT_pt_15_20_New",-1.0, 1.0);
	RooRealVar* BDT_pt_20_50 = new RooRealVar("BDT_pt_20_50_New","BDT_pt_20_50_New",-1.0, 1.0);


	RooRealVar* w = new RooRealVar("pthatweightNew","pthatweightNew",0.,12.) ;
	RooDataSet* ds = new RooDataSet();
	RooDataSet* dsMC = new RooDataSet();   
	std::cout<<"Created dataset"<<std::endl;
	RooDataHist* dh = new RooDataHist();   
	RooDataHist* dhMC = new RooDataHist();   
	std::cout<<"Created roodatahists"<<std::endl;
	RooPlot* frame = new RooPlot();
	std::cout<<"Created d"<<std::endl;
	RooHist* datahist = new RooHist();
	RooCurve* modelcurve = new RooCurve();
	std::cout<<"Created roocurve"<<std::endl;

	//weightgen = weightgen_pp;
	weightmc  = weightmc_pp;
	if(usePbPb){
		weightgen = weightgen_PbPb;
		//weightmc = weightmc_PbPb;
		if(doweight == 0 ) weightmc = "1"; 
		if(doweight == 1) weightmc = "pthatweightNew";
		std::cout<<"weights defined"<<std::endl;	

	}

	TString _prefix = "";
	TString _isMC = "data";
	if(isMC) _isMC = "mcAsData";
	TString _isPbPb = "pp";
	if(usePbPb) _isPbPb = "PbPb";
	TString _postfix = "";
	if(weightdata!="1") _postfix = "_EFFCOR";


	std::cout<<"Passed AddFriends"<<std::endl;
	TString outputf;
	outputf = Form("%s",outputfile.Data());
	std::cout<<"Formed outputf "<<outputf<<std::endl;
	TFile* outf = new TFile(outputf.Data(),"recreate");
	outf->cd();

	std::cout<<"varExp= "<<_varExp<<std::endl;

	ds = new RooDataSet(Form("ds%d",_count),"",skimtree_new, RooArgSet(*mass, *pt, *y));
	dsMC = new RooDataSet(Form("dsMC%d",_count),"",skimtreeMC_new,RooArgSet(*mass, *pt, *w, *y));
	//    TString ptbinning;

	std::vector<std::string> background = {"1st", "2nd", "3rd"};
	std::vector<std::vector<double>> background_syst;
	// std::vector<std::string> signal = {"3gauss", "fixed", "scal", "perr", "merr", "scal+", "scal-"};
	//std::vector<std::string> signal = {"3gauss", "fixed","scal+", "scal-"};
	//	 std::vector<std::string> signal = {"3gauss"};
	std::vector<std::string> signal = {"3gauss", "fixed", "scal+", "scal-"};
	//std::vector<std::string> signal = {"3gauss", "1gauss","fixed", "scal", "perr", "merr"};
	// std::vector<std::string> signal = {"3gauss", "1gauss", "fixed_mean", "scaling", "-_error", "+_error"};
	std::vector<std::vector<double>> signal_syst;
	std::vector<std::vector<double>> general_syst;
	std::vector<double> nominal_yields;
	std::vector<std::vector<double>> back_syst_rel_values;
	std::vector<std::vector<double>> sig_syst_rel_values;



	for(int i=0;i<nBins_check;i++)
	{
		_count++;
		std::cout<<"I'm in the for"<<std::endl;
		TCanvas* c= new TCanvas(Form("c%d",_count),"",600,600);
		TCanvas* cMC= new TCanvas(Form("cMC%d",_count),"",600,600);
		//		c->cd();
		//	if(fitOnSaved == 0){
		drawOpt = 1;
		// tree with selecitons applied, create RooDataSet
		/*	if(isMC==1) {
			skimtree = (TTree*)nt->CopyTree(Form("%s*(%s&&%s>%f&&%s<%f)*(1/%s)&&BgenNew==23333",weightmc.Data(),seldata.Data(),varExp.Data(),_ptBins[i],varExp.Data(),_ptBins[i+1],weightdata.Data()), "skimtree");
			}
			else{
			skimtree = (TTree*)nt->CopyTree(   Form("(%s&&%s>%f&&%s<%f)*(1/%s)",                seldata.Data(),varExp.Data(),_ptBins[i],varExp.Data(),_ptBins[i+1],weightdata.Data()), "skimtree");		
			std::cout<<"copy tree for data"<<std::endl;
			}
			fout->cd();*/ 
		//skimtree->Write(); 
		//fout->Close();
		//      skimtreeMC = (TTree*)ntMC->CopyTree(Form("%s*(%s&&%s>%f&&%s<%f)",weightmc.Data(),Form("%s&&BgenNew==23333",selmc.Data()),varExp.Data(),_ptBins[i],varExp.Data(),_ptBins[i+1]), "skimtreeMC");
		//			ds = new RooDataSet(Form("ds%d",_count),"",skimtree_new, RooArgSet(*mass, *pt));
		std::cout<<_varExp<<" bin = "<<ptBins_check[i]<<" "<<ptBins_check[i+1]<<std::endl;
		//      return;
		RooDataSet* ds_cut;
		std::cout<< "Pass 1" << std::endl;
		//if(doubly==0) ds_cut = new RooDataSet(Form("ds_cut%d",_count),"",ds, RooArgSet(*mass, *pt, *y, *BDT_pt_5_10, *BDT_pt_10_15, *BDT_pt_15_20,*BDT_pt_20_50), Form("%s>=%f&&%s<=%f&&BmassNew>%f&&BmassNew<%f&&(BptNew > 5 &&  BptNew < 10 && BDT_pt_5_10_New > 0.10)&&(BptNew > 10 &&  BptNew < 15 && BDT_pt_10_15_New > 0.10)&&(BptNew > 15 &&  BptNew < 20 && BDT_pt_15_20_New > 0.10)&&(BptNew > 20 &&  BptNew < 50 && BDT_pt_20_50_New > 0.10)",varExp.Data(),ptBins_check[i],varExp.Data(),ptBins_check[i+1],minhisto, maxhisto));
		if(doubly==0) ds_cut = new RooDataSet(Form("ds_cut%d",_count),"",ds, RooArgSet(*mass, *pt, *y, *BDT_pt_5_10, *BDT_pt_10_15, *BDT_pt_15_20,*BDT_pt_20_50), Form("%s>=%f&&%s<=%f&&BmassNew>%f&&BmassNew<%f",varExp.Data(),ptBins_check[i],varExp.Data(),ptBins_check[i+1],minhisto, maxhisto)); 

		if(doubly==1)ds_cut = new RooDataSet(Form("ds_cut%d",_count),"",ds, RooArgSet(*mass, *pt, *y,  *BDT_pt_5_10, *BDT_pt_10_15, *BDT_pt_15_20,*BDT_pt_20_50), Form("%s>=%f&&%s<=%f&&BmassNew>%f&&BmassNew<%f",varExp.Data(),ptBins_check[i],varExp.Data(),ptBins_check[i+1],minhisto, maxhisto)); 
		if(doubly==2)ds_cut = new RooDataSet(Form("ds_cut%d",_count),"",ds, RooArgSet(*mass, *pt, *y,  *BDT_pt_5_10, *BDT_pt_10_15, *BDT_pt_15_20,*BDT_pt_20_50), Form("%s>=%f&&%s<=%f&&BmassNew>%f&&BmassNew<%f",varExp.Data(),ptBins_check[i],varExp.Data(),ptBins_check[i+1],minhisto, maxhisto)); 
		std::cout<< "Pass 2" << std::endl;
		//  RooRealVar* w = (RooRealVar*) data->addColumn(wFunc) ;
		//      RooDataSet wdata(data->GetName(),data->GetTitle(),data,*data->get(),0,w->GetName()) ;
		//dsMC = new RooDataSet(Form("dsMC%d",_count),"",skimtreeMC_new,RooArgSet(*mass, *pt));
		RooDataSet* dsMC_cut;
		if(doubly==0) dsMC_cut = new RooDataSet(Form("dsMC_cut%d",_count),"",dsMC, RooArgSet(*mass, *pt, *w, *y), Form("%s>=%f&&%s<=%f&&BmassNew>%f&&BmassNew<%f",varExp.Data(),ptBins_check[i],varExp.Data(),ptBins_check[i+1],minhisto, maxhisto), "pthatweightNew"); 
		if(doubly==1) dsMC_cut = new RooDataSet(Form("dsMC_cut%d",_count),"",dsMC, RooArgSet(*mass, *pt, *w, *y), Form("%s>=%f&&%s<=%f&&BmassNew>%f&&BmassNew<%f",varExp.Data(),ptBins_check[i],varExp.Data(),ptBins_check[i+1],minhisto, maxhisto), "pthatweightNew"); 
		if(doubly==2) dsMC_cut = new RooDataSet(Form("dsMC_cut%d",_count),"",dsMC, RooArgSet(*mass, *pt, *w, *y), Form("%s>=%f&&%s<=%f&&BmassNew>%f&&BmassNew<%f",varExp.Data(),ptBins_check[i],varExp.Data(),ptBins_check[i+1],minhisto, maxhisto), "pthatweightNew"); 
		// return;
		std::cout<<"Really created datasets"<<std::endl;
		// create RooDataHist
		h = new TH1D(Form("h%d",_count),"",nbinsmasshisto,minhisto,maxhisto);
		hMC = new TH1D(Form("hMC%d",_count),"",nbinsmasshisto,minhisto,maxhisto);
		std::cout<<"TH1D "<<std::endl;
		if(isMC==1) skimtree_new->Project(Form("h%d",_count),"BmassNew",Form("%s*(%s&&%s>%f&&%s<%f && BgenNew == 23333)*(1/%s)",weightmc.Data(),seldata.Data(),varExp.Data(),_ptBins[i],varExp.Data(),_ptBins[i+1],weightdata.Data()));
		else        skimtree_new->Project(Form("h%d",_count),"BmassNew",   Form("(%s&&%s>%f&&%s<%f)*(1/%s)",                seldata.Data(),varExp.Data(),_ptBins[i],varExp.Data(),_ptBins[i+1],weightdata.Data()));
		std::cout<<"isMC? "<<isMC<<std::endl;
//		std::string project_str=Form("(%s&&%s>%f&&%s<%f)*(1/%s)",                seldata.Data(),varExp.Data(),_ptBins[i],varExp.Data(),_ptBins[i+1],weightdata.Data());
//		std::cout<<"string= "<<project_str<<std::endl;
//		std::string project_str_mc= Form("%s*(%s&&%s>%f&&%s<%f && BgenNew == 23333)*(1/%s)",weightmc.Data(),seldata.Data(),varExp.Data(),_ptBins[i],varExp.Data(),_ptBins[i+1],weightdata.Data());
//		std::cout<<"string isMC= "<<project_str_mc<<std::endl;
		skimtreeMC_new->Project(Form("hMC%d",_count),"BmassNew",Form("%s*(%s&&%s>%f&&%s<%f)",weightmc.Data(),Form("%s&&BgenNew==23333",selmc.Data()),varExp.Data(),_ptBins[i],varExp.Data(),_ptBins[i+1]));
		std::cout<<"skim_tree "<<std::endl;
		dh = new RooDataHist(Form("dh%d",_count),"",*mass,Import(*h));
		dhMC = new RooDataHist(Form("dhMC%d",_count),"",*mass,Import(*hMC));
		std::cout<<"RooDataHist "<<std::endl;
		h->SetAxisRange(0,h->GetMaximum()*1.4,"Y");
		outputw->import(*ds);
		outputw->import(*dsMC);
		outputw->import(*dh);
		outputw->import(*dhMC);
		std::cout<<"import "<<std::endl;
		//}
		/*if(fitOnSaved == 1){
		  drawOpt = 1;
		  inputw = (RooWorkspace*)inf->Get("w");
		  ds = (RooDataSet*)inputw->data(Form("ds%d",_count));
		  dsMC = (RooDataSet*)inputw->data(Form("dsMC%d",_count));
		  dh = (RooDataHist*)inputw->data(Form("dh%d",_count));
		  dhMC = (RooDataHist*)inputw->data(Form("dhMC%d",_count));
		  }*/


		//RooFitResult* f = fit(c, cMC, ds, dsMC, dh, dhMC, mass, frame, _ptBins[i], _ptBins[i+1], isMC, isPbPb, centmin, centmax, npfit);
		RooFitResult* f = fit("", "", tree, c, cMC, ds_cut, dsMC_cut, dh, dhMC, mass, frame, ptBins_check[i], ptBins_check[i+1], isMC, isPbPb, centmin, centmax, npfit);
		//RooFitResult* f = fit("sigonly", "", tree, c, cMC, ds_cut, dsMC_cut, dh, dhMC, mass, frame, ptBins_check[i], ptBins_check[i+1], isMC, isPbPb, centmin, centmax, npfit);

		//  return;

		std::cout << "Now Finall We Validate Our Fits" << std::endl;
	//	validate_fit(w_val, tree, _varExp, full,centmin, centmax, ptBins_check[i], ptBins_check[i+1]);
		//std::cout << "Now Finall We Scan Our Significance" << std::endl;

		//scan_significance(w_val, tree, _varExp, full,centmin, centmax, ptBins_check[i], ptBins_check[i+1]);

		/*
		for(int q= 0; q < 100; q++){
			std::cout << "BRO NOW IS WORKING ON  " << q << std::endl; 
			validate_fit(w_val, tree, _varExp, full,q);
		}
		*/
	
	
		//datahist = frame->getHist("ds");
		//TGraphAsymmErrors* datagraph = static_cast<TGraphAsymmErrors*>(datahist);
		modelcurve = frame->getCurve(Form("model%d",_count));

		RooRealVar* fitYield = static_cast<RooRealVar*>(f->floatParsFinal().at(f->floatParsFinal().index(Form("nsig%d",_count))));
		double yield = fitYield->getVal();
		yieldRec = fitYield->getVal();


		RooRealVar* BackGround = static_cast<RooRealVar*>(f->floatParsFinal().at(f->floatParsFinal().index(Form("nbkg%d",_count))));
		MyBackground = BackGround->getVal();

		nominal_yields.push_back(yield);
		double yieldErr = fitYield->getError();
		printf("yield: %f, yieldErr: %f\n", yield, yieldErr);
		yieldErr = yieldErr*_ErrCor;
		if(fitOnSaved == 0){
			TH1D* htest = new TH1D(Form("htest%d",_count),"",nbinsmasshisto,minhisto,maxhisto);
			TString sideband = "(abs(BmassNew-5.367)>0.2&&abs(BmassNew-5.367)<0.3";
			//    	skimtree_new->Project(Form("htest%d",_count),"BmassNew",Form("%s&&%s&&%s>%f&&%s<%f)*(1/%s)",sideband.Data(),seldata.Data(),varExp.Data(),_ptBins[i],varExp.Data(),_ptBins[i+1],weightdata.Data()));
			std::cout<<"yield bkg sideband: "<<htest->GetEntries()<<std::endl;
		}
		std::cout<<"pt bin = "<<ptBins_check[i]<<" "<<ptBins_check[i+1]<<std::endl;
		std::cout<<"count= "<<i<<std::endl;
		std::cout<<"pt bin = "<<ptBins_check[i]<<" "<<ptBins_check[i+1]<<std::endl;
		std::cout<<"count= "<<i<<std::endl;
		std::cout<<"pt bin = "<<ptBins_check[i]<<" "<<ptBins_check[i+1]<<std::endl;
		std::cout<<"count= "<<i<<std::endl;

		if(_varExp!="HiBin"){
			hPt->SetBinContent(i+1,yield/(ptBins_check[i+1]-ptBins_check[i]));
			hPt->SetBinError(i+1,yieldErr/(ptBins_check[i+1]-ptBins_check[i]));
		}




		else{
			hPt->SetBinContent(i+1,yield/(_ptBins[1] - _ptBins[0]));
			hPt->SetBinError(i+1,yieldErr/(_ptBins[1] - _ptBins[0]));
		}




		if(f->floatParsFinal().index(Form("nsig%d",_count)) != -1){
			RooRealVar* fitMean = static_cast<RooRealVar*>(f->floatParsFinal().at(f->floatParsFinal().index(Form("mean%d",_count))));
			hMean->SetBinContent(i+1,fitMean->getVal());
			hMean->SetBinError(i+1,fitMean->getError());  
		}



		TLatex* tex_pt;
	//	TLatex* tex_hibin;
		TLatex* tex_y;
		//tex = new TLatex(0.55,0.85,Form("%.0f < p_{T} < %.0f GeV/c",_ptBins[i],_ptBins[i+1]));
		//if(varExp=="abs(By)") tex = new TLatex(0.55,0.85,Form("%.1f < y < %.1f",_ptBins[i],_ptBins[i+1]));
		if(varExp=="BptNew"){
			tex_pt = new TLatex(0.21,0.75,Form("%d < p_{T} < %d GeV/c",(int)ptBins_check[i],(int)ptBins_check[i+1]));
		//	if(centmin==0&&centmax==90)tex_hibin = new TLatex(0.21,0.62,"Cent. 0-90%");
		//	if(centmin==0&&centmax==30)tex_hibin = new TLatex(0.21,0.62,"Cent. 0-30%");
		//	if(centmin==30&&centmax==90)tex_hibin = new TLatex(0.21,0.62,"Cent. 30-90%");
			tex_y = new TLatex(0.21,0.69,"|y| < 2.4"); 
			//	    tex = new TLatex(0.21,0.72,Form("%.0f < p_{T} < %.0f GeV/c",_ptBins[i],_ptBins[i+1]));

		}

		std::cout<<"CHEGUEI AQUI"<<std::endl;

		if(varExp=="HiBin"){
			tex_pt = new TLatex(0.21,0.75,"5 < p_{T} < 50 GeV/c");
	//		tex_hibin = new TLatex(0.21,0.62,Form("Cent. %d-%d%%",(int)ptBins_check[i]/2,(int)ptBins_check[i+1]/2));
			tex_y = new TLatex(0.21,0.69,"|y| < 2.4");

			//	    tex = new TLatex(0.21,0.72,Form("%.0f < p_{T} < %.0f GeV/c",_ptBins[i],_ptBins[i+1]));
		}
		if(varExp=="abs(By)"){
			tex_pt = new TLatex(0.21,0.75,"5 < p_{T} < 50 GeV/c");
	//		if(doubly==0)tex_hibin = new TLatex(0.21,0.62,"Cent. 0-90%");
	//		if(doubly==1)tex_hibin = new TLatex(0.21,0.62,"Cent. 0-30%");
	//		if(doubly==2)tex_hibin = new TLatex(0.21,0.62,"Cent. 30-90%");
			tex_y = new TLatex(0.21,0.69,Form("%.1f < |y| < %.1f",ptBins_check[i],ptBins_check[i+1]));
		}

		tex_pt->SetNDC();
		tex_pt->SetTextFont(42);
		tex_pt->SetTextSize(0.045);
		tex_pt->SetLineWidth(2);
		if(drawOpt ==1) tex_pt->Draw();
		tex_pt->Draw();

		/*
		tex_hibin->SetNDC();
		tex_hibin->SetTextFont(42);
		tex_hibin->SetTextSize(0.045);
		tex_hibin->SetLineWidth(2);
		if(drawOpt ==1) tex_hibin->Draw();
		tex_hibin->Draw();
		*/
		std::cout<<"CHEGUEI AQUI"<<std::endl;


		/*   	if(varExp=="abs(By)") tex = new TLatex(0.21,0.72,Form("%.1f < y < %.1f",_ptBins[i],_ptBins[i+1]));
				tex->SetNDC();
				tex->SetTextFont(42);
				tex->SetTextSize(0.045);
				tex->SetLineWidth(2);
				tex->Draw();*/

		//tex = new TLatex(0.75,0.80,"|y| < 2.4");
		//if(varExp=="abs(By)") tex = new TLatex(0.75,0.80,"15 < p_{T} < 50 GeV.c");
		tex_y->SetNDC();
		tex_y->SetTextFont(42);
		tex_y->SetTextSize(0.045);
		tex_y->SetLineWidth(2);
		if(drawOpt ==1)	tex_y->Draw();
		tex_y->Draw();

		std::cout<<"CHEGUEI AQUI"<<std::endl;
	
		CMS_lumi(c,19011,0);
		c->Update();
		std::cout<<"Canvas Updated"<<std::endl;
		c->SaveAs(Form("%s%s/%s_%s_%d%s_%s_%d%d_doubly%d_%.0f_%0.f",outplotf.Data(),_prefix.Data(),_isMC.Data(),_isPbPb.Data(),_count,_postfix.Data(),_varExp.Data(),(int)ptBins_check[i],(int)ptBins_check[i+1], doubly,centmin,centmax)+tree+".pdf");
		c->Update();
		c->SaveAs(Form("%s%s/%s_%s_%d%s_%s_%d%d_doubly%d_%.0f_%0.f",outplotf.Data(),_prefix.Data(),_isMC.Data(),_isPbPb.Data(),_count,_postfix.Data(),_varExp.Data(),(int)ptBins_check[i],(int)ptBins_check[i+1], doubly,centmin,centmax)+tree+".png");
		//c->SaveAs(Form("%s%s/%s_%s_%d%s_%s_%d%d_doubly%d_%.0f_%0.f",outplotf.Data(),_prefix.Data(),_isMC.Data(),_isPbPb.Data(),_count,_postfix.Data(),_varExp.Data(),(int)ptBins_check[i],(int)ptBins_check[i+1], doubly,centmin,centmax)+tree+".C");
		// c->SaveAs(Form("%s%s/%s_%s_%d%s_%d%d.png",outplotf.Data(),_prefix.Data(),_isMC.Data(),_isPbPb.Data(),_count,_postfix.Data(), (int)ptBins_check[i],(int)ptBins_check[i+1]));
		// c->SaveAs(Form("%s%s/%s_%s_%d%s_%d%d.C",outplotf.Data(),_prefix.Data(),_isMC.Data(),_isPbPb.Data(),_count,_postfix.Data(), (int)ptBins_check[i],(int)ptBins_check[i+1]));
		cMC->SaveAs(Form("%s%s/%s_%s_%d%s_%s_%d%d_doubly%d_%.0f_%0.f",outplotf.Data(),_prefix.Data(),"mc",_isPbPb.Data(),_count,_postfix.Data(),_varExp.Data(), (int)ptBins_check[i], (int)ptBins_check[i+1], doubly,centmin,centmax)+tree+".pdf");
		cMC->SaveAs(Form("%s%s/%s_%s_%d%s_%s_%d%d_doubly%d_%.0f_%0.f",outplotf.Data(),_prefix.Data(),"mc",_isPbPb.Data(),_count,_postfix.Data(),_varExp.Data(), (int)ptBins_check[i], (int)ptBins_check[i+1], doubly,centmin,centmax)+tree+".png");

		// cMC->SaveAs(Form("%s%s/%s_%s_%d%s_%d%d.png",outplotf.Data(),_prefix.Data(),"mc",_isPbPb.Data(),_count,_postfix.Data(), (int)ptBins_check[i], (int)ptBins_check[i+1]));
		//cMC->SaveAs(Form("%s%s/%s_%s_%d%s_%d%d.C",outplotf.Data(),_prefix.Data(),"mc",_isPbPb.Data(),_count,_postfix.Data(), (int)ptBins_check[i], (int)ptBins_check[i+1]));
		//        return;
		/*	    TH1* h = dh->createHistogram("BmassNew");
				h->GetEntries();
				h->Sumw2(kFALSE);
				h->SetBinErrorOption(TH1::kPoisson);roofitB.C
				TCanvas* cpull= new TCanvas(Form("cpull%d",_count),"",600,600);
				cpull->cd();
				TGraphAsymmErrors* pullgraph = new TGraphAsymmErrors();
				pullgraph->SetTitle("");
				pullgraph->SetMaximum(5);
				pullgraph->SetMinimum(-5);
				pullgraph->SetMarkerSize(1.55); pullgraph->SetMarkerStyle(20); pullgraph->SetLineColor(1); pullgraph->SetLineWidth(4);
				double x; double xfit; double y; double yfit;
				double binerr;
				for(int b = 0; b < h->GetNbinsX(); b++){
				modelcurve->GetPoint(modelcurve->findPoint(h->GetBinCenter(b+1)),xfit,yfit);
				binerr = h->GetBinContent(b+1) > yfit ? h->GetBinErrorLow(b+1) : h->GetBinErrorUp(b+1);
				pullgraph->SetPoint(b,h->GetBinCenter(b+1),(h->GetBinContent(b+1)-yfit)/binerr);
		//pullgraph->SetPointEYlow(b,h->GetBinErrorLow(b+1)/binerr);
		//pullgraph->SetPointEYhigh(b,h->GetBinErrorUp(b+1)/binerr);
		pullgraph->SetPointEYlow(b,1);
		pullgraph->SetPointEYhigh(b,1);
		}
		TLine* line = new TLine(5., 0., 6., 0.);
		line->SetLineStyle(9);
		line->SetLineWidth(4);
		line->SetLineColor(kGreen+1);
		line->Draw();
		pullgraph->Draw();
		cpull->SaveAs(Form("%s%s/%s_%s_%d%s_pull.pdf",outplotf.Data(),_prefix.Data(),_isMC.Data(),_isPbPb.Data(),_count,_postfix.Data()));*/

		RooCurve* modelcurve_back = new RooCurve();
		std::vector<double> back_variation; 
		std::vector<double> back_err;  
		RooCurve* modelcurve_signal = new RooCurve();
		std::vector<double> signal_variation; 
		std::vector<double> signal_err;
		std::vector<double> general_err;  
		double max_signal=0.;
		double max_back=0.;
		double full_err=0;



		

		if(syst==1){

			for(int j=0; j<background.size(); j++){
				RooFitResult* f_back = fit("background", background[j], tree, c, cMC, ds_cut, dsMC_cut, dh, dhMC, mass, frame, ptBins_check[i], ptBins_check[i+1], isMC, isPbPb, centmin, centmax, npfit);
				tex_pt->Draw();
			//	tex_hibin->Draw();
				tex_y->Draw();

				c->SaveAs(Form("%s%s/%s_%s_%d%s_%s_%d%d_%s_doubly_%d",outplotf.Data(),_prefix.Data(),_isMC.Data(),_isPbPb.Data(),_count,_postfix.Data(),_varExp.Data(),(int)ptBins_check[i],(int)ptBins_check[i+1],background[j].c_str(), doubly)+tree+".pdf");
				c->SaveAs(Form("%s%s/%s_%s_%d%s_%s_%d%d_%s_doubly_%d",outplotf.Data(),_prefix.Data(),_isMC.Data(),_isPbPb.Data(),_count,_postfix.Data(),_varExp.Data(),(int)ptBins_check[i],(int)ptBins_check[i+1],background[j].c_str(), doubly)+tree+".png");

				modelcurve_back = frame->getCurve(Form("model%d",_count));
				RooRealVar* fitYield_back = static_cast<RooRealVar*>(f_back->floatParsFinal().at(f_back->floatParsFinal().index(Form("nsig%d",_count))));
				back_variation.push_back(fitYield_back->getVal());
				back_err.push_back(abs(((yield-fitYield_back->getVal())/yield)*100));
				if(abs(((yield-fitYield_back->getVal())/yield)*100)>max_back) max_back=abs(((yield-fitYield_back->getVal())/yield)*100);

			}

			general_err.push_back(max_back);

			for(int j=0; j<signal.size(); j++){
				RooFitResult* f_signal = fit("signal", signal[j], tree, c, cMC, ds_cut, dsMC_cut, dh, dhMC, mass, frame, ptBins_check[i], ptBins_check[i+1], isMC, isPbPb, centmin, centmax, npfit);
				tex_pt->Draw();
			//	tex_hibin->Draw();
				tex_y->Draw();

				c->SaveAs(Form("%s%s/%s_%s_%d%s_%s_%d%d_%s_doubly_%d",outplotf.Data(),_prefix.Data(),_isMC.Data(),_isPbPb.Data(),_count,_postfix.Data(),_varExp.Data(),(int)ptBins_check[i],(int)ptBins_check[i+1],signal[j].c_str(), doubly)+tree+".pdf");
				cMC->SaveAs(Form("%s%s/%s_%s_%d%s_%s_%d%d_%s_doubly_%d",outplotf.Data(),_prefix.Data(),"mc",_isPbPb.Data(),_count,_postfix.Data(),_varExp.Data(), (int)ptBins_check[i], (int)ptBins_check[i+1],signal[j].c_str(), doubly)+tree+".pdf");
				c->SaveAs(Form("%s%s/%s_%s_%d%s_%s_%d%d_%s_doubly_%d",outplotf.Data(),_prefix.Data(),_isMC.Data(),_isPbPb.Data(),_count,_postfix.Data(),_varExp.Data(),(int)ptBins_check[i],(int)ptBins_check[i+1],signal[j].c_str(), doubly)+tree+".png");
				cMC->SaveAs(Form("%s%s/%s_%s_%d%s_%s_%d%d_%s_doubly_%d",outplotf.Data(),_prefix.Data(),"mc",_isPbPb.Data(),_count,_postfix.Data(),_varExp.Data(), (int)ptBins_check[i], (int)ptBins_check[i+1],signal[j].c_str(), doubly)+tree+".png");

				modelcurve_signal = frame->getCurve(Form("model%d",_count));
				RooRealVar* fitYield_signal = static_cast<RooRealVar*>(f_signal->floatParsFinal().at(f_signal->floatParsFinal().index(Form("nsig%d",_count))));
				signal_variation.push_back(fitYield_signal->getVal());
				signal_err.push_back(abs(((yield-fitYield_signal->getVal())/yield)*100));
				if(abs(((yield-fitYield_signal->getVal())/yield)*100)>max_signal) max_signal=abs(((yield-fitYield_signal->getVal())/yield)*100);


			}

			general_err.push_back(max_signal);

			full_err=sqrt(max_back*max_back+max_signal*max_signal);

			general_err.push_back(full_err);



		}


		if(syst==1){

			background_syst.push_back(back_variation);
			signal_syst.push_back(signal_variation);
			back_syst_rel_values.push_back(back_err);
			sig_syst_rel_values.push_back(signal_err);
			general_syst.push_back(general_err);

		}
	
		

	//validate_fit(outputw,1,Form("model%d",_count),Form("nsig%d",_count));
	}


	hMean->Write();
	hPt->Write();
	outf->Close();	

	std::ofstream myfile;
	myfile.open ("systematics_"+tree+".txt");


	if(syst==1){ 
		for(int i=0; i<nBins_check; i++){
			std::cout<<"pt bin = "<<ptBins_check[i]<<" "<<ptBins_check[i+1]<<std::endl;
			std::cout<<"nominal yield = "<<nominal_yields[i]<<std::endl;
			std::cout<<"count= "<<i<<std::endl;

			for(int j=0; j<background.size(); j++){
				std::cout<<" back sys in bin "<<i<<" ; with pdf "<< background[j] << " ="<<background_syst[i][j]<<std::endl;
				std::cout<<" back sys in bin "<<i<<" ; with pdf "<< background[j] << " ="<<back_syst_rel_values[i][j]<<" % "<<std::endl;
				myfile<<" back sys in bin "<<i<<" ; with pdf "<< background[j] << " ="<<back_syst_rel_values[i][j]<<" % "<<std::endl;

			}
			for(int j=0; j<signal.size(); j++){
				std::cout<<" signal sys in bin "<<i<<" ; with pdf "<< signal[j] << " ="<<signal_syst[i][j]<<std::endl;
				std::cout<<" signal sys in bin "<<i<<" ; with pdf "<< signal[j] << " ="<<sig_syst_rel_values[i][j]<<" % "<<std::endl;
				myfile<<" signal sys in bin "<<i<<" ; with pdf "<< signal[j] << " ="<<sig_syst_rel_values[i][j]<<" % "<<std::endl;

			}
		}

	}

	myfile.close();


	std::cout<<"CHEGUEI AQUI"<<std::endl;

	/*if(fitOnSaved == 1){
	  outf->Close();	
	  return;
	  }
	  ntMC->Project("hPtMC",varExp.Data(),TCut(weightmc)*(TCut(selmc.Data())&&"(BgenNew==23333)"));
	  divideBinWidth(hPtMC);
	//	ntGen->Project("hPtGen","Gpt",TCut(weightgen)*(TCut(selmcgen.Data())));
	//	divideBinWidth(hPtGen);*/

	TCanvas* cPt =  new TCanvas("cPt","",600,600);
	cPt->SetLogy();
	hPt->SetXTitle("B_{s} p_{T} (GeV/c)");
	hPt->SetYTitle("Uncorrected dN(B_{s})/dp_{T}");
	hPt->Draw();


	std::vector<std::string> labels_back = {"Linear", "2nd Poly", "3rd Poly"};
	std::vector<std::string> col_name_back;
	if(tree=="ntphi"&&_varExp=="BptNew") col_name_back= {"Background Model","5$<p_T<$10", "10$<p_T<$15", "15$<p_T<$20", "20$<p_T<$50"};
	if(_varExp=="HiBin") col_name_back= {"Background Model","0\\%$<Cent<$30\\%", "30\\%$<Cent<$90\\%"};
	if(tree=="ntKp" && _varExp=="BptNew") col_name_back= {"Background Model","5$<p_T<$7", "7$<p_T<$10", "10$<p_T<$15", "15$<p_T<$20", "20$<p_T<$30", "30$<p_T<$40", "40$<p_T<$50", "50$<p_T<$60" };

	std::vector<std::string> labels_signal = {"Triple Gaussian", "Fixed Mean", "Increased $\\sigma$ scaling factor", "Decreased $\\sigma$ scaling factor"};
	//		std::vector<std::string> labels_signal = {"Triple Gaussian", "Fixed Mean","Double Gaussian + $\\sigma$ scaling factor" , "Increased $\\sigma$ scaling factor", "Decreased $\\sigma$ scaling factor"};
	std::vector<std::string> col_name_signal;
	if(tree=="ntphi"&&_varExp=="BptNew") col_name_signal= {"Signal Model","5$<p_T<$10", "10$<p_T<$15", "15$<p_T<$20", "20$<p_T<$50"};
	if(_varExp=="HiBin") col_name_signal = {"Signal Model","0\\%$<Cent<$30\\%", "30\\%$<Cent<$90\\%"};
	if(tree=="ntKp" && _varExp=="BptNew") col_name_signal= {"Signal Model","5$<p_T<$7", "7$<p_T<$10", "10$<p_T<$15", "15$<p_T<$20", "20$<p_T<$30", "30$<p_T<$40", "40$<p_T<$50", "50$<p_T<$60" };



	std::vector<std::string> labels_general = {"Background", "Signal", "Total"};
	std::vector<std::string> col_name_general;
	if(tree=="ntphi"&&_varExp=="BptNew") col_name_general= {"Systematic Source","5$<p_T<$10", "10$<p_T<$15", "15$<p_T<$20", "20$<p_T<$50"};
	if(_varExp=="HiBin") col_name_general = {"Systematic Source","0\\%$<Cent<$30\\%", "30\\%$<Cent<$90\\%"};
	if(tree=="ntKp" && _varExp=="BptNew") col_name_general= {"Systematic Source","5$<p_T<$7", "7$<p_T<$10", "10$<p_T<$15", "15$<p_T<$20", "20$<p_T<$30", "30$<p_T<$40", "40$<p_T<$50", "50$<p_T<$60" };


	if(syst==1 && full==0){
		latex_table("background_systematics_table_"+std::string (_varExp.Data())+"_"+std::string (tree.Data()), nBins_check+1,  (int)(1+background.size()),  col_name_back,labels_back,back_syst_rel_values, "Background PDF Systematic Errors");
		latex_table("signal_systematics_table_"+std::string (_varExp.Data())+"_"+std::string (tree.Data()), nBins_check+1, (int)(1+signal.size()),    col_name_signal, labels_signal,sig_syst_rel_values, "Signal PDF Systematic Errors");
		latex_table("general_systematics_table_"+std::string (_varExp.Data())+"_"+std::string (tree.Data()),   nBins_check+1, 4 , col_name_general, labels_general, general_syst, "Overall PDF Variation Systematic Errors");	

	}



	cout << "Final Background = " << MyBackground << endl;

	cout << "Final Yield = " << yieldRec << endl;

}
