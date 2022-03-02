#include "TAxis.h"
//#include "uti.h"
#include "RooWorkspace.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooChi2Var.h"
#include "RooHist.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include <RooBifurGauss.h>
#include <fstream>
#include <string>
#include <iomanip>
#include "RooMCStudy.h"
#include <RooMinuit.h>


#include "RooMinimizer.h"
#include "RooProfileLL.h"
#include "RooMinuit.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooAbsReal.h"

using namespace RooFit;
using namespace std;

#define BSUBS_MASS 5.36682
#define BP_MASS 5.27915
#define JPSI_MASS 3.096916 


bool drawSup = 1;

void clean0 (TH1D* h);
int    drawOpt=0;
double _ErrCor=1;
double setparam0=100.;
double setparam2=0.02;
double setparam3=0.06;
double fixparam1=BSUBS_MASS;

double minhisto=2.3;
double maxhisto=3.8;
int nbinsmasshisto=50;
double binwidthmass=(maxhisto-minhisto)/nbinsmasshisto;

 double x_1 = 0.58;
 double x_2 = 0.95;
 double y_2 = 0.92;
 double y_1;
 double y_space = 0.04;
 double displace = 0.17;
TString seldata;
TString selmc;
TString selmcgen;
TString collisionsystem;
Float_t hiBinMin,hiBinMax,centMin,centMax;
Int_t _count=0;
RooWorkspace* inputw = new RooWorkspace();
RooWorkspace* outputw = new RooWorkspace("w");
RooWorkspace* w_val= new RooWorkspace("w_val");


Double_t yield;
Double_t yieldErr; 

Double_t FitMean;
Double_t FitMeanErr; 

Double_t FitWidth1; 
Double_t FitWidth1Err; 


Double_t FitWidth2; 
Double_t FitWidth2Err; 


RooAddPdf* model;
double SignalWidth = 0.20;


RooFitResult *fit(TString variation, TString pdf, TCanvas* c, RooDataSet* ds,  RooDataHist* dh, RooRealVar* mass, RooPlot* &outframe)
{

	cout << "Inside" << endl;
	cout<<"total data: "<<ds->numEntries()<<endl;
	TH1* h = dh->createHistogram("dimuon_mass");
	h->Sumw2(kFALSE);
	h->SetBinErrorOption(TH1::kPoisson);
	h->SetMarkerSize(1.55);
	h->SetMarkerStyle(20);
	h->SetLineColor(1);
	h->SetLineWidth(4);

	c->cd();
	RooPlot* frame = mass->frame();
	TPad *p1 = new TPad("p1","p1",0.,0.28,1.,0.99);
	// TPad *p1 = new TPad("p1","p1",0.05,0.05,0.99,0.99);
	p1->SetBorderMode(1); 
	p1->SetFrameBorderMode(0); 
	p1->SetBorderSize(2);
	p1->SetBottomMargin(0.10);
	p1->Draw(); 

	double shiftYpad = 0.08;

	TPad *p2 = new TPad("p2","p2",0.,0.10,1.,0.24+ shiftYpad);// 0.26 
	p2->SetTopMargin(0.);    
	p2->SetBorderMode(0);
	p2->SetBorderSize(2); 
	p2->SetFrameBorderMode(0); 
	p2->SetTicks(1,1); 
	//  p2->SetBottomMargin(0.2);
	p2->Draw();


	double mass_peak = 3.096916;

	double n_signal_initial = ds->sumEntries(TString::Format("abs(dimuon_mass-%g)<0.5",mass_peak)) - ds->sumEntries(TString::Format("abs(dimuon_mass-%g)<0.10&&abs(dimuon_mass-%g)>0.05",mass_peak,mass_peak));
	if(n_signal_initial<0)
		n_signal_initial=1;

	double n_combinatorial_initial = ds->sumEntries() - n_signal_initial;
	p1->cd();
	RooRealVar mean(Form("mean%d",_count),"",3.09,3.05,3.15) ;
	RooRealVar mean2(Form("mean2%d",_count),"",3.686,3.6,3.8) ;

	RooRealVar sigma1(Form("sigma1%d",_count),"",0.2,0.0,0.6) ;
	RooRealVar sigma2(Form("sigma2%d",_count),"",0.2,0.01,0.4) ;
	RooRealVar sigma3(Form("sigma3%d",_count),"",0.2,0.01,0.5) ;
	RooRealVar sigma4(Form("sigma4%d",_count),"",0.2,0.01,0.5) ;
	
	RooGaussian sig1(Form("sig1%d",_count),"",*mass,mean,sigma1);  
	RooGaussian sig2(Form("sig2%d",_count),"",*mass,mean,sigma2);  
	RooGaussian sig3(Form("sig3%d",_count),"",*mass,mean2,sigma3);  
	RooGaussian sig4(Form("sig4%d",_count),"",*mass,mean2,sigma4);  




	//Double Crystal Ball//


	RooRealVar cbmean(Form("cbmean%d",_count), Form("cbmean%d",_count) , 3.09, 2.9, 3.3) ;
	RooRealVar cbsigma1("cbsigma1", "cbsigma1" , 0.4, 0.0, 0.8);
	RooRealVar cbsigma2("cbsigma2", "cbsigma2" , 0.4, 0.0, 0.8) ;

	RooRealVar n1("n1","", 5.1,0,100);
	RooRealVar n2("n2","", 5.1,0,100);
	
	RooRealVar alpha1("alpha1","", 1.3,0,10);
	RooRealVar alpha2("alpha2","", 1.3,0,10);

	RooCBShape cball1("cball1", "crystal ball1", *mass, cbmean, cbsigma1, alpha1, n1);
	RooCBShape cball2("cball2", "crystal ball2", *mass, cbmean, cbsigma2, alpha2, n2);


	


	cout << "Pass Here Bro" << endl;


	RooRealVar c1(Form("c1%d",_count),"",1.,0.,5.);
	RooRealVar c2(Form("c2%d",_count),"",1.,0.,5.);


	//RooRealVar c2(Form("c2%d",_count),"",1.,0.,5.) ;
	RooGenericPdf sig1_gen(Form("sig1_gen%d",_count),"", Form("exp(-0.5*(((dimuon_mass-mean%d)*(dimuon_mass-mean%d))/((c1%d*sigma1%d)*(c1%d*sigma1%d))))",_count, _count, _count, _count, _count, _count), RooArgSet(*mass, mean, sigma1, c1));
	RooGenericPdf sig2_gen(Form("sig2_gen%d",_count),"", Form("exp(-0.5*(((dimuon_mass-mean%d)*(dimuon_mass-mean%d))/((c1%d*sigma2%d)*(c1%d*sigma2%d))))", _count, _count, _count, _count, _count, _count), RooArgSet(*mass, mean, sigma2, c1));

	RooGenericPdf sig3_gen(Form("sig3_gen%d",_count),"", Form("exp(-0.5*(((dimuon_mass-mean2%d)*(dimuon_mass-mean2%d))/((c2%d*sigma3%d)*(c2%d*sigma3%d))))",_count, _count, _count, _count, _count, _count), RooArgSet(*mass, mean2, sigma3, c2));
	RooGenericPdf sig4_gen(Form("sig4_gen%d",_count),"", Form("exp(-0.5*(((dimuon_mass-mean2%d)*(dimuon_mass-mean2%d))/((c2%d*sigma4%d)*(c2%d*sigma4%d))))", _count, _count, _count, _count, _count, _count), RooArgSet(*mass, mean2, sigma4, c2));

	//  RooGenericPdf sig1(Form("sig1%d",_count),"", Form("(1/(c1%d*sigma1%d*sqrt(2*pi)))*exp(-0.5*(((dimuon_mass-mean%d)*(dimuon_mass-mean%d))/((c1%d*sigma1%d)*(c1%d*sigma1%d))))",_count, _count, _count, _count, _count, _count, _count, _count), RooArgSet(*mass, mean, sigma1, c1));
	// RooGenericPdf sig2(Form("sig2%d",_count),"", Form("(1/(c1%d*sigma2%d*sqrt(2*pi)))*exp(-0.5*(((dimuon_mass-mean%d)*(dimuon_mass-mean%d))/((c1%d*sigma2%d)*(c1%d*sigma2%d))))", _count, _count, _count, _count, _count, _count, _count, _count), RooArgSet(*mass, mean, sigma2, c1));
	RooRealVar sig1frac(Form("sig1frac%d",_count),"",0.2,0.,1.);
	RooRealVar sig2frac(Form("sig2frac%d",_count),"",0.4,0.,1.);
	RooAddPdf* sig;
	RooAddPdf* sig_2;

	if(variation=="signal" && (pdf=="scal"|| pdf==""|| pdf=="scal-"|| pdf=="scal+")) sig = new RooAddPdf(Form("sig%d",_count),"",RooArgList(sig1_gen,sig2_gen),sig1frac);
	if((variation=="" && pdf=="") || (variation=="sigonly" && pdf=="")  || variation== "background"|| (variation=="signal" && (pdf=="1gauss"|| pdf=="fixed" || pdf=="merr" || pdf == "perr"))) sig = new RooAddPdf(Form("sig%d",_count),"",RooArgList(sig1,sig2),sig1frac);
	if(variation=="signal" && pdf=="3gauss") sig = new RooAddPdf(Form("sig%d",_count), "", RooArgList(sig1, sig2, sig3), RooArgList(sig1frac, sig2frac));
	//RooRealVar a0(Form("a0%d",_count),"a0",0.1,0.,1.);
	//RooRealVar a1(Form("a1%d",_count),"a1",0.1,0.,1.);
	//RooRealVar a2(Form("a2%d",_count),"a2",0.1,0.,1.);
	//RooChebychev bkg(Form("bkg%d",_count),"",*mass,RooArgSet(a0,a1,a2));
	RooRealVar a0(Form("a0%d",_count),"",0.1,-1e4,1.0);
	RooRealVar a1(Form("a1%d",_count),"",0.1,-1e4,1e4);
	RooRealVar a2(Form("a2%d",_count),"",1e0,-1e4,1e4);
	//	RooRealVar a3(Form("a3%d",_count),"",1e0,-1e4,1e4);
	RooPolynomial bkg_1st(Form("bkg%d",_count),"",*mass,RooArgSet(a0));//2nd orderpoly
	RooPolynomial bkg_2nd(Form("bkg%d",_count),"",*mass,RooArgSet(a0,a1));//2nd orderpoly
	RooPolynomial bkg_3rd(Form("bkg%d",_count),"",*mass,RooArgSet(a0,a1,a2));//2nd orderpoly
	RooRealVar lambda(Form("lambda%d", _count), "lambda",-2., -5., 5);
	RooExponential bkg(Form("bkg%d",_count),"",*mass,lambda);//2nd orderpoly
	// RooPolynomial bkg(Form("bkg%d",_count),"",*mass,a0);//linear
	//  RooPolynomial bkg(Form("bkg%d",_count),"",*mass,RooArgSet(a0,a1));//linear

	RooRealVar nsig(Form("nsig%d",_count),"",n_signal_initial,-1,ds->sumEntries()*3);
	RooRealVar nsig_2(Form("nsig_2%d",_count),"",400,-1,2000);

	RooRealVar nbkg(Form("nbkg%d",_count),"",n_combinatorial_initial,0.,ds->sumEntries()*10);
	//	RooRealVar nbkg(Form("nbkg%d",_count),"",n_combinatorial_initial,0.,200);
	std::cout<<"sumEntries: "<<ds->sumEntries()<<std::endl;
	std::cout<<"sumEntries: "<<ds->sumEntries()<<std::endl;
	std::cout<<"sumEntries: "<<ds->sumEntries()<<std::endl;
	std::cout<<"sumEntries: "<<ds->sumEntries()<<std::endl;
	std::cout<<"sumEntries: "<<ds->sumEntries()<<std::endl;
	std::cout<<"sumEntries: "<<ds->sumEntries()<<std::endl;
	std::cout<<"sumEntries: "<<ds->sumEntries()<<std::endl;
	std::cout<<"sumEntries: "<<ds->sumEntries()<<std::endl;
	std::cout<<"sumEntries: "<<ds->sumEntries()<<std::endl;
	std::cout<<"sumEntries: "<<ds->sumEntries()<<std::endl;
	std::cout<<"sumEntries: "<<ds->sumEntries()<<std::endl;
	std::cout<<"sumEntries: "<<ds->sumEntries()<<std::endl;
	std::cout<<"sumEntries: "<<ds->sumEntries()<<std::endl;

//	sig = new RooAddPdf(Form("sig%d",_count),"",RooArgList(sig1));


/*
	if(variation=="" && pdf=="") model = new RooAddPdf(Form("model%d",_count),"",RooArgList(*sig,bkg),RooArgList(nsig,nbkg));

	if(variation=="background" && pdf=="1st") model = new RooAddPdf(Form("model%d",_count),"",RooArgList(*sig,bkg_1st),RooArgList(nsig,nbkg));
	//	mean.setConstant();
	//sigma1.setConstant();
	if(pdf!="1gauss"){
//		sigma2.setConstant();
//		sig1frac.setConstant();
	}
	if(variation=="signal" && pdf=="3gauss"){
//		sigma3.setConstant();
//		sig2frac.setConstant();  
	}
	if(variation=="signal" && pdf=="fixed") mean.setConstant();

	if(variation=="signal" && pdf=="scal+") {
		c1.setConstant();
		c1.setVal(1.1);
		//c1.setVal(c1.getVal()+c1.getError());
		//fitResult = model->fitTo(*ds,Save(), Minos() , Extended(kTRUE));

	}
	if(variation=="signal" && pdf=="scal-") {
		c1.setConstant();
		c1.setVal(0.9);
		//c1.setVal(c1.getVal()- c1.getError());
		// fitResult = model->fitTo(*ds,Save(), Minos() , Extended(kTRUE));

	}
	*/
	//model = new RooAddPdf(Form("model%d",_count),"",RooArgList(*sig,bkg_3rd),RooArgList(nsig,nbkg));
//	//sig = new RooAddPdf(Form("sig%d",_count),"",RooArgList(sig1_gen));
//	sig1frac.setConstant();
//	sig1frac.setVal(1.0);


//	sig_2 = new RooAddPdf(Form("sig_2%d",_count),"",RooArgList(sig3,sig4),sig2frac);
//	sig = new RooAddPdf(Form("sig%d",_count), "", RooArgList(sig1, sig2, sig3,sig4), RooArgList(sig1frac, sig2frac));
//	model = new RooAddPdf(Form("model%d",_count),"",RooArgList(*sig,*sig_2,bkg),RooArgList(nsig,nsig_2,nbkg));
//	sig = new RooAddPdf(Form("sig%d",_count),"",RooArgList(sig1,sig2),sig1frac);
//	model = new RooAddPdf(Form("model%d",_count),"",RooArgList(*sig,bkg),RooArgList(nsig,nbkg));
//	sig1frac.setConstant();
//	sig1frac.setVal(0.0);


	sig = new RooAddPdf(Form("sig%d",_count),"",RooArgList(cball1,cball2),sig1frac);

	//	sig = new RooAddPdf(Form("sig%d",_count),"",RooArgList(cball1));
	//	sig = new RooAddPdf(Form("sig%d",_count),"",cball1);

	model = new RooAddPdf(Form("model%d",_count),"",RooArgList(*sig,bkg),RooArgList(nsig,nbkg));

	frame->SetMaximum(nsig.getVal()*0.9);

	/*
	// test removing bakcground
	//  RooRealVar* lambda =  (w->var(Form("lambda%d", _count)));
	//  RooRealVar* nbkg   =  (w->var(Form("nbkg%d", _count)));
	lambda.setConstant();
	nbkg.setVal(0.);
	nbkg.setConstant();


	RooFitResult* fitResult = model->fitTo(*ds,Save(), Minos() , Extended(kTRUE),Range(5.3,5.45));
	*/

	//	cout << "nbkg->getVal() = " << nbkg.getVal() << endl;

	//cout << "-----------------------------------------------------------------------" << endl;

	cout << "Begin FitTo AAAAAA" << endl;

	RooFitResult* fitResult = model->fitTo(*ds,Save(), Minos(),  RooFit::PrintLevel(0) , Extended(kTRUE));
	//RooFitResult* fitResult = model->fitTo(*ds,Save(), Minos() , Extended(kTRUE),Range(5.34,5.40));


	

	cout << "After FitTo BBBBB" << endl;


	cout << "sigma1.getVal() = " << sigma1.getVal() << endl; 
	

	//	fitResult->printArg();
	//	fitResult->printValue();
	//cout << "-----------------------------------------------------------------------" << endl;


	w_val->import(*model);
	w_val->import(nsig);
	w_val->import(cbmean);

	cout << "nbkg->getVal() = " << nbkg.getVal() << endl;
	


	RooAbsReal* nll = model->createNLL(*ds);
	double log_likelihood= nll->getVal();

	/*
	//Likehood Scan//

	TString rootoutputname = Form("signscan/root/sigscan_%.0f_%.0f_%.0f_%.0f.root",centmin,centmax,ptmin,ptmax);
	TFile* outf = new TFile(rootoutputname.Data(), "RECREATE");

	cout << "Pass Inside 2" << endl;


	TCanvas* cnll = new TCanvas("cnll", "", 600, 600);
	cnll->cd();
	cout << "Pass Inside 2.5" << endl;

	double MaxFrame;
	double MaxYield;

	if(ptmin == 5) MaxFrame = 50;
	if(ptmin == 15) MaxFrame = 200;
	if(ptmin == 20) MaxFrame = 200;
	if(ptmin == 10 && ptmax == 50 && centmin == 0 && centmax ==90) MaxFrame = 500;
	if(ptmin == 10 && ptmax == 50 && centmin == 0 && centmax ==30) MaxFrame = 300;
	if(ptmin == 10 && ptmax == 50 && centmin == 30 && centmax ==90) MaxFrame = 200;

	if(ptmin == 10 && ptmax == 15) MaxFrame = 200;


	if(ptmin == 5) MaxYield = 20;
	if(ptmin == 15) MaxYield = 70;
	if(ptmin == 20) MaxYield = 70;
	if(ptmin == 10 && ptmax == 50 && centmin == 0 && centmax ==90) MaxYield = 200;
	if(ptmin == 10 && ptmax == 50 && centmin == 0 && centmax ==30) MaxYield = 150;
	if(ptmin == 10 && ptmax == 50 && centmin == 30 && centmax ==90) MaxYield = 80;

	if(ptmin == 10 && ptmax == 15) MaxYield = 70;


	RooMinuit(*nll).migrad();
	RooPlot* frame1 = nsig.frame(RooFit::Bins(25), RooFit::Range(0., MaxYield));
	nll->plotOn(frame1, RooFit::ShiftToZero(), RooFit::LineColor(kBlack), RooFit::LineWidth(4), RooFit::Normalization(2., RooAbsReal::Raw));
	cout << "Pass Inside 3" << endl;




	cnll->cd();
	frame1->SetMinimum(0);
	frame1->SetMaximum(MaxFrame);
	frame1->GetXaxis()->SetTitle("N_{sig}^{(B^{0}_{s})}");
	frame1->GetYaxis()->SetTitle("-2log#lambda");
	frame1->GetXaxis()->SetTitleOffset(1.5);

	frame1->GetXaxis()->CenterTitle();
	frame1->GetYaxis()->CenterTitle();
	frame1->GetYaxis()->SetTitleOffset(1.6);
	frame1->SetTitle("");
	frame1->Draw();
	TString OutputName = Form("signscan/plots/sigscan_%.0f_%.0f_%.0f_%.0f.png",centmin,centmax,ptmin,ptmax);
	cout << "Pass Inside 4" << endl;

	cnll->SaveAs(OutputName.Data());



	std::cout << "DONE SCANNING" <<  std::endl;

	//DONE Scanning//
	*/

	//  a0.setVal(0.);
	// a0.setConstant();
	ds->plotOn(frame,Name("ds"),Binning(nbinsmasshisto),MarkerSize(1),MarkerStyle(20),MarkerColor(1),LineColor(1),LineWidth(4),LineColor(1));//draw an transparent hist
	std::cout<<"GETS hERE?"<<std::endl;

	model->plotOn(frame,Name(Form("bkg%d",_count)),Components(bkg),Normalization(1.0,RooAbsReal::RelativeExpected),Precision(1e-6),DrawOption("L"),LineStyle(7),LineColor(4),LineWidth(4));

	if(pdf!="1gauss"){
		model->plotOn(frame,Name(Form("sig%d",_count)),Components(*sig),Normalization(1.0,RooAbsReal::RelativeExpected),Precision(1e-6),DrawOption("L"),FillStyle(3002),FillColor(kOrange-3),LineStyle(7),LineColor(kOrange-3),LineWidth(4));
		std::cout<<"GETS hERE not 1gauss?"<<std::endl;
		model->plotOn(frame,Name(Form("sigF%d",_count)),Components(*sig),Normalization(1.0,RooAbsReal::RelativeExpected),Precision(1e-6),DrawOption("F"),FillStyle(3002),FillColor(kOrange-3),LineStyle(7),LineColor(kOrange-3),LineWidth(4));
		std::cout<<"GETS hERE not 1gauss?"<<std::endl;
	}
	else{
		model->plotOn(frame,Name(Form("sig%d",_count)),Components(sig1),Normalization(1.0,RooAbsReal::RelativeExpected),Precision(1e-6),DrawOption("L"),FillStyle(3002),FillColor(kOrange-3),LineStyle(7),LineColor(kOrange-3),LineWidth(4));
		std::cout<<"GETS hERE else?"<<std::endl;
		model->plotOn(frame,Name(Form("sigF%d",_count)),Components(sig1),Normalization(1.0,RooAbsReal::RelativeExpected),Precision(1e-6),DrawOption("F"),FillStyle(3002),FillColor(kOrange-3),LineStyle(7),LineColor(kOrange-3),LineWidth(4));
		std::cout<<"GETS hERE else?"<<std::endl;
	}
	//model->plotOn(frame,Name(Form("sig_2%d",_count)),Components(sig1),Normalization(1.0,RooAbsReal::RelativeExpected),Precision(1e-6),DrawOption("L"),FillStyle(3002),FillColor(kGreen-3),LineStyle(7),LineColor(kGreen-3),LineWidth(4));
	
	model->plotOn(frame,Name(Form("model%d",_count)),Normalization(1.0,RooAbsReal::RelativeExpected),Precision(1e-6),DrawOption("L"),LineColor(2),LineWidth(4));

	std::cout<<"GETS hERE?"<<std::endl;
	frame->SetMaximum(nsig.getVal()*0.9);
	//	frame->SetMaximum((h->GetBinContent(h->GetMaximumBin())+h->GetBinError(h->GetMaximumBin()))*1.8);
	model->paramOn(frame,Layout(x_2+0.5, x_2+0.5, y_1+0.16), Format("NEU",AutoPrecision(3)));
	 //model->paramOn(frame,Layout(0.65, x_2, y_1-0.06), Format("NEU",AutoPrecision(3)));
	
	
	frame->getAttText()->SetTextSize(0.00);
	/*  frame->getAttFill()->SetFillStyle(0);
		frame->getAttLine()->SetLineWidth(0);*/
	frame->SetTitle("");
	//frameMC->SetXTitle("m_{B} (GeV/c^{2})");
	//frame->SetXTitle("m_{J/#psi(#mu#mu)#it{#phi}(KK)} (GeV/c^{2})");
	//frame->SetYTitle("Events / (20 MeV/c^{2})");
	frame->SetYTitle("Events");
	
	frame->GetXaxis()->CenterTitle();
	frame->GetYaxis()->CenterTitle();
	frame->GetXaxis()->SetTitleOffset(1.1);
	frame->GetYaxis()->SetTitleOffset(1.4);
	frame->GetXaxis()->SetTitleSize(0.04);
	frame->GetXaxis()->SetTitleSize(0.0);
	frame->GetYaxis()->SetTitleSize(0.05);
	// frame->GetXaxis()->SetTitleFont(42);
	frame->GetXaxis()->SetLabelFont(42);
	frame->GetYaxis()->SetTitleFont(42);
	frame->GetYaxis()->SetLabelFont(42);
	frame->GetXaxis()->SetLabelSize(0.045);
	frame->GetXaxis()->SetLabelSize(0.0);
	frame->GetYaxis()->SetLabelSize(0.05);
	//	frame->SetXTitle("");
	//	frame->GetXaxis()->SetLabelSize(0.0);

	//  frameMC->GetYaxis()->SetMaxDigits(2);
	frame->SetStats(0);
	frame->SetMaximum(5000);
	frame->SetMaximum(nsig.getVal()*0.20);

	cout << "Pass Here Bro" << endl;

	(frame->GetXaxis())->SetRangeUser(minhisto,maxhisto);
	frame->GetXaxis()->SetNdivisions(-50205);	
	frame->Draw();
	cout << "Pass Here Before Pull" << endl;
	
	RooHist* pull_hist = frame->pullHist("ds",Form("model%d",_count));



	RooPlot* pull_plot = mass->frame();
	(pull_plot->GetXaxis())->SetRangeUser(minhisto,maxhisto);
	cout << "Here 2" << endl;

	RooRealVar x("x","",1e3,0,1e6);
	x.setVal(0.);
	cout << "Here 3" << endl;
	
	RooGenericPdf* line_ref = new RooGenericPdf("ref_0", "ref_0", "x", RooArgList(x));
	
	cout << "Here 4" << endl;
	line_ref->plotOn(pull_plot, LineStyle(7), LineColor(13), LineWidth(2));  

	cout << "Here bro" << endl;

	pull_plot->addPlotable(static_cast<RooPlotable*>(pull_hist),"XP");
	pull_plot->SetTitle("");

	pull_plot->SetXTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
	std::cout<<"pull start"<<std::endl;


	pull_plot->SetYTitle("Pull");
	pull_plot->GetYaxis()->SetTitleFont(42);  
	pull_plot->GetYaxis()->SetTitleFont(42);  
	pull_plot->GetYaxis()->SetTitleSize(0.17);
	pull_plot->GetXaxis()->SetTitleOffset(0.4);

	pull_plot->GetYaxis()->CenterTitle();
	pull_plot->GetXaxis()->SetTitleSize(0.20);
	pull_plot->GetXaxis()->SetTitleOffset(1.0);
	pull_plot->GetXaxis()->CenterTitle();
	//	pull_plot->GetXaxis()->SetTitleOffset(1.4);

	/*  pull_plot->GetYaxis()->SetLabelFont(42);
		pull_plot->GetYaxis()->SetLabelSize(0.13);
		pull_plot->GetYaxis()->SetLabelOffset(0.005);*/
	/*  pull_plot->SetTitleFont(42, "Y");
		pull_plot->SetTitleOffset(1.3, "Y");  
		pull_plot->SetTitleSize(0.5, "Y");*/
	pull_plot->GetXaxis()->SetLabelFont(42);
	pull_plot->GetXaxis()->SetLabelOffset(0.01);
	//pull_plot->GetXaxis()->SetTitleSize(0.2);
	//pull_plot->GetXaxis()->SetTitleOffset(1.);
	pull_plot->GetXaxis()->SetLabelFont(42);
	pull_plot->GetXaxis()->SetLabelSize(0.15);
	//pull_plot->GetXaxis()->SetTitleFont(42);
	pull_plot->GetXaxis()->SetTickLength(0.15);
	pull_plot->GetYaxis()->SetLabelFont(42);
	pull_plot->GetYaxis()->SetLabelOffset(0.01);
	pull_plot->GetYaxis()->SetLabelSize(0.15);
	//  pull_plot->GetYaxis()->SetTitleSize(0.1);
	pull_plot->GetYaxis()->SetNdivisions(305);
	pull_plot->GetXaxis()->SetNdivisions(-50205);
	//  c->SaveAs("dataset_afterfit_test.pdf");
	//	h->Draw("same e");
	// c->RedrawAxis();
	std::cout<<"pull done"<<std::endl;
	/*
	   mass->setRange("signal",5.2,5.6);
	   RooAbsReal* sigIntegral = sig->createIntegral(*mass,NormSet(*mass),Range("signal"));
	   double sigIntegralErr = sigIntegral->getPropagatedError(*fitResult);
	   RooAbsReal* bkgIntegral = bkg.createIntegral(*mass,NormSet(*mass),Range("signal"));
	   double bkgIntegralErr = bkgIntegral->getPropagatedError(*fitResult);
	   cout<<"nsig: "<<nsig.getVal()<<endl;
	//cout<<"nsig: "<<model->getParameters(*mass)->getRealValue("nsig")<<endl;
	cout<<"nsig error: "<<nsig.getError()<<endl;
	cout<<"sig integral: "<<sigIntegral->getVal()<<endl;
	cout<<"sig integral error: "<<sigIntegralErr<<endl;
	cout<<"nbkg: "<<nbkg.getVal()<<endl;
	//cout<<"nbkg: "<<model->getParameters(*mass)->getRealValue("nbkg")<<endl;
	cout<<"nbkg error: "<<nbkg.getError()<<endl;
	cout<<"bkg integral: "<<bkgIntegral->getVal()<<endl;
	cout<<"bkg integral error: "<<bkgIntegralErr<<endl;*/
	//cout<<"#expected events: "<<model->expectedEvents(RooArgSet(*mass))<<endl;
	//cout<<"model integral: "<<model->getNormIntegral(RooArgSet(*mass))->getVal()<<endl;
	//fitResult->Print("v");

	//RooHist* datahist = new RooHist();
	//datahist = frame->getHist("ds");
	//TGraphAsymmErrors* datagraph = static_cast<TGraphAsymmErrors*>(datahist);
	//RooCurve* modelcurve = new RooCurve();
	//modelcurve = frame->getCurve("model");

	yield = nsig.getVal();
	yieldErr = nsig.getError();
	
	FitMean = cbmean.getVal();
	FitMeanErr = cbmean.getError();
	 
	FitWidth1 = cbsigma1.getVal();
	FitWidth1Err = cbsigma1.getError();
	 
	FitWidth2 = cbsigma2.getVal();
	FitWidth2Err = cbsigma2.getError();


	

	TH1D* fh = (TH1D*)h->Clone("fh");
/*
	double nexpected = model->expectedEvents(RooArgSet(*mass));
	double dataArr[nbinsmasshisto]; double dataErrArr[nbinsmasshisto]; double fitArr[nbinsmasshisto];
	double val = 0;
	for(int i = 0; i < nbinsmasshisto; i++){
		dataArr[i] = h->GetBinContent(i+1);
		dataErrArr[i] = h->GetBinError(i+1);
		mass->setVal(h->GetBinCenter(i+1));
		val = model->getVal(RooArgSet(*mass))*1./nbinsmasshisto*nexpected;
		fitArr[i] = val;
		fh->SetBinContent(i+1, fitArr[i]);
		fh->SetBinError(i+1, sqrt(fitArr[i]));
	}
	cout<<"frame->chiSquare: "<<frame->chiSquare(Form("model%d",_count),Form("ds_cut%d",_count), fitResult->floatParsFinal().getSize())<<endl;
	cout<<"pars= "<<fitResult->floatParsFinal().getSize()<<endl;
	cout<<"chi2: "<<frame->chiSquare(Form("model%d",_count),Form("ds%d",_count))<<endl;
	//RooChi2Var chi2_lowstat("chi2_lowstat","chi2",*model,*dh);
	//cout<<chi2_lowstat.getVal()<<endl;


	double chiRoo = frame->chiSquare(Form("model%d",_count),Form("ds_cut%d",_count), fitResult->floatParsFinal().getSize());
	double chi2Std = 0;
	double chi2Neyman = 0;
	double chi2Peason = 0;
	double chi2BakerCousins = 0;
	chi2Cal(dataArr, dataErrArr, fitArr, nbinsmasshisto, chi2Std, chi2Neyman, chi2Peason, chi2BakerCousins);
	printf("chi2 Standard: %f\n",chi2Std);
	printf("chi2 Neyman: %f\n",chi2Neyman);
	printf("chi2 Peason: %f\n",chi2Peason);
	printf("chi2 Baker & Cousins: %f\n",chi2BakerCousins);
*/
	TLegend *leg = new TLegend(0.62,0.55,0.89,0.75,NULL,"brNDC"); 
	if(drawOpt == 1) {
		leg = new TLegend(0.60,0.7,0.89,0.89,NULL,"brNDC");
	}
	if(drawOpt == 0) {
		leg = new TLegend(0.50/1.20+0.10,0.42,0.89+0.10,0.86,NULL,"brNDC");
	}
		leg = new TLegend(0.60/1.20+0.10,0.42,0.89+0.10,0.86,NULL,"brNDC");


	leg->SetBorderSize(0);
	leg->SetTextSize(0.05);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	//leg->AddEntry(frame->findObject(Form("ds%d",_count),"Data","pl");
	leg->AddEntry(h,"Data","pl");
	leg->AddEntry(frame->findObject(Form("model%d",_count)),"Fit","l");
	leg->AddEntry(frame->findObject(Form("sig%d",_count)),"J/#psi Signal","f");
//	leg->AddEntry(frame->findObject(Form("sig_2%d",_count)),"#psi(2S) Signal","f");	
	//leg->AddEntry(frame->findObject(Form("bkg%d",_count)),"Combinatorial","l");
	leg->AddEntry(frame->findObject(Form("bkg%d",_count)),"Background","l");


	TLatex* texcms = new TLatex(0.21,0.88,"CMS");
	texcms->SetNDC();
	texcms->SetTextAlign(13);
	texcms->SetTextFont(62);
	texcms->SetTextSize(0.04);
	texcms->SetLineWidth(2);
	TLatex* texpre = new TLatex(0.30,0.88,"Preliminary");
	texpre->SetNDC();
	texpre->SetTextAlign(13);
	texpre->SetTextFont(52);
	texpre->SetTextSize(0.04);
	texpre->SetLineWidth(2);
	TLatex* texsup = new TLatex(0.21,0.815,"Supplementary");
	texsup->SetNDC();
	texsup->SetTextAlign(13);
	texsup->SetTextFont(52);
	texsup->SetTextSize(0.05);
	texsup->SetLineWidth(2);
	TLatex* texB;
	texB = new TLatex(0.21,0.7,"J/#psi");
	texB->SetNDC();
	texB->SetTextFont(62);
	texB->SetTextSize(0.075);
	texB->SetLineWidth(2);

	TLatex* texCol;
	if(collisionsystem=="pp"||collisionsystem=="PP"||collisionsystem=="ppInc") texCol= new TLatex(0.95,0.94, Form("28.0 pb^{-1} (%s 5.02 TeV)","pp"));
	else texCol= new TLatex(0.945,0.94, Form("%s 5.02 TeV (300 pb^{-1})","pp"));
	texCol->SetNDC();
	texCol->SetTextAlign(32);
	texCol->SetTextSize(0.055);
	texCol->SetTextFont(42);
/*
	float posit = 0.58;
	int nDOF = nbinsmasshisto-(fitResult->floatParsFinal().getSize());
	float nChi2 = chi2BakerCousins/(nbinsmasshisto-(fitResult->floatParsFinal().getSize()));
	int nDigit_chi2BakerCousins = 2;
	int nDigit_nChi2 = 2;
	chi2BakerCousins = roundToNdigit(chi2BakerCousins);
	nChi2 = roundToNdigit(nChi2);
	nDigit_chi2BakerCousins = sigDigitAfterDecimal(chi2BakerCousins);
	nDigit_nChi2 = sigDigitAfterDecimal(nChi2);
	TLatex* texChi = new TLatex(0.65,0.64, Form("#chi^{2}/nDOF = %2f", chiRoo));
	texChi->SetNDC();

	texChi->SetTextAlign(12);
	texChi->SetTextSize(0.03);
	texChi->SetTextFont(42);
*/

	nsig.setVal(0.);
	nsig.setConstant();
	RooFitResult* fitResult_nosig = model->fitTo(*ds,Save());
	RooAbsReal* nll_nosig = model->createNLL(*ds);
	double log_likelihood_nosig= nll_nosig->getVal();
	double real_significance = sqrt(2*(-log_likelihood+log_likelihood_nosig));
	std::cout<<"REAL SIGNIFICANCE= "<<real_significance<<std::endl;
	std::cout<<"***********************************************************************"<<std::endl;


	//RooFitResult* fitResult = model->fitTo(*ds,Save(),Minos());
	//ds->plotOn(frame,Name(Form("ds%d",_count)),Binning(nbinsmasshisto),MarkerSize(1.55),MarkerStyle(20),LineColor(1),LineWidth(4));
/*
	Double_t Significance =  real_significance;
	//	Double_t Significance =  yield/TMath::Sqrt(bkgd+yield);
	int nDigit_Significance = 3;
	Significance = roundToNdigit(Significance);
	nDigit_Significance = sigDigitAfterDecimal(Significance);
	TLatex* texSig = new TLatex(0.65,0.67,Form("Significance = %.*f", nDigit_Significance, Significance));
	cout<<"Significance = "<<Significance<<endl;
	texSig->SetNDC();
	texSig->SetTextFont(42);
	texSig->SetTextSize(0.03);
	texSig->SetLineWidth(2);

*/

	leg->Draw("same");
	//texcms->Draw();
	//texpre->Draw("same");
	texB->Draw();
	c->cd();
	p1->cd();
	leg->Draw("same");
	//	texcms->Draw("");
	//	texpre->Draw("same");
	//if(drawSup) texsup->Draw();
	texB->Draw();
	//texCol->Draw();
	//if(1) {

//	texSig->Draw("SAME");
		//texChi->Draw("same");
		//texYield->Draw("SAME");
	


	p2->cd();

	pull_plot->GetYaxis()->SetTitleOffset(0.4);
	pull_plot->SetYTitle("Pull");


	pull_plot->Draw();
	p1->cd();
	outframe = frame;
	outputw->import(*model);

//	float init_mean  = 3.096916;
//	SignalWidth = 0.2;

	mass->setRange("signal",JPSI_MASS-SignalWidth, JPSI_MASS+SignalWidth);
	RooAbsReal* RangeBakground = bkg.createIntegral(*mass,NormSet(*mass),Range("signal")); 
	RooAbsReal* RangeSig = sig->createIntegral(*mass,NormSet(*mass),Range("signal")); 



	//cout << "RangeBakground = " << RangeBakground << endl;
	yield = nsig.getVal() * RangeSig->getVal();
	yieldErr = nsig.getError() * RangeSig->getVal();

	double Calback = RangeBakground->getVal() * nbkg.getVal();
	double StatSig = yield/sqrt(yield + Calback);

	
	


	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.060);

	double JPsiMean = cbmean.getVal();
	double JPsiWidth = cbsigma1.getVal();

	lat->DrawLatex(0.23,0.47,Form("S = %.1f",yield));	
//	lat->DrawLatex(0.23,0.47,Form("S_{Err} = %.1f",yieldErr));		
	lat->DrawLatex(0.23,0.39,Form("B = %.1f",Calback));	
//	lat->DrawLatex(0.23,0.31,Form("Mean = %.2f GeV/c^{2}",JPsiMean));	
//	lat->DrawLatex(0.23,0.23,Form("Width = %.2f GeV/c^{2}f",JPsiWidth));	


	lat->Draw("SAME");

	return fitResult;

}

void clean0(TH1D* h)
{
	for (int i=1;i<=h->GetNbinsX();i++)
	{
		if(h->GetBinContent(i)==0) h->SetBinError(i,1);
	}
}

double ErrorOnSigma(double width, double errwidth, double smear, double errsmearing){
	double squarederroronsigma=(1+smear)*(1+smear)*errwidth*errwidth+width*width*errsmearing*errsmearing;
	double erroronsigma=TMath::Sqrt(squarederroronsigma);
	return erroronsigma;
}


void latex_table(std::string filename, int n_col, int n_lin, std::vector<std::string> col_name, std::vector<std::string> labels, 
		std::vector<std::vector<double> > numbers, std::string caption)
{
	std::ofstream file_check;
	std::ofstream file;

	//Begin Document                                                                                                                               

	file.open(filename + ".tex");
	file_check.open(filename + "_check.tex");

	file_check << "\\documentclass{article}" << std::endl;
	//file << "\\usepackage[utf8]{inputenc}" << std::endl;     
	file_check << "\\usepackage{rotating}" << std::endl;                                                                                   
	file_check << "\\usepackage{cancel}" << std::endl;
	file_check << "\\usepackage{geometry}" << std::endl;
	file_check << "\\usepackage{booktabs}" << std::endl;
	file_check << "\\geometry{a4paper, total={170mm,257mm}, left=20mm, top=20mm,}" << std::endl;

	file_check << "\\title{Bs/B+}" << std::endl;
	file_check << "\\author{Julia Silva}" << std::endl;
	file_check << "\\date{September 2019}" << std::endl;
	file_check << "\\begin{document}" << std::endl;
	file_check << "\\maketitle" << std::endl;

	// Create table                                                                                                                                
	//file_check << "\\begin{table}[!h]" << std::endl;
	//file_check << "\\centering" << std::endl;                                                                                                         
	//setup table size                                                                                                                             
	std::string col="c";
	//std::setprecision(3);

	for(int i=1; i<n_col; i++)
		col+="|c";

	file_check << "\\begin{sidewaystable}"<< std::endl;
	file_check << "\\begin{tabular}{"+col+"}" << std::endl;
	file_check << "\\toprule" << std::endl;
	file << "\\begin{tabular}{"+col+"}" << std::endl;
	file << "\\toprule" << std::endl;

	for(int c=0; c<n_col-1; c++){
		file << col_name[c] << " & ";
		file_check << col_name[c] << " & ";
	}

	file << col_name[n_col-1] << " \\\\ \\midrule" << std::endl;
	file_check << col_name[n_col-1] << " \\\\ \\midrule" << std::endl;


	for(int i=1; i<n_lin; i++)
	{
		file << labels[i-1] << " & ";
		file_check << labels[i-1] << " & ";

		for(int c=1; c<n_col-1; c++){
			file << std::setprecision(3)<< numbers[c-1][i-1]<< " \\% & ";
			file_check << std::setprecision(3) << numbers[c-1][i-1]<< " \\% & ";

		}

		file << std::setprecision(3)<< numbers[n_col-2][i-1]<< " \\% \\\\" << std::endl;
		file_check << std::setprecision(3)<< numbers[n_col-2][i-1]<< " \\% \\\\" << std::endl;
	}

	file << "\\bottomrule" << std::endl;
	file_check << "\\bottomrule" << std::endl;

	//End Table                                                                                                                                    
	file << "\\end{tabular}" << std::endl;
	file_check << "\\end{tabular}" << std::endl;
	file_check << "\\caption{"+caption+"}" << std::endl;
	file_check << "\\end{sidewaystable}"<< std::endl;

	//file_check << "\\end{table}" << std::endl;
	//End document                                                                                                                                 

	file_check << "\\end{document}" << std::endl;

	file.close();
	file_check.close();

	system(("pdflatex " + filename + "_check.tex").c_str());
	system(("open " + filename + "_check.pdf").c_str());
}


void validate_fit(RooWorkspace* w, int Opt)
{
	std::cout << "Now Perform Check on Fit" << std::endl;
	RooRealVar dimuon_mass = *(w->var("dimuon_mass"));
	RooAbsPdf* model  = w->pdf(Form("model%d",_count));
	//RooDataSet* data = (RooDataSet*) w->data("data");

	//model->fitTo(*data);

	std::vector<RooRealVar> params;
	params.push_back(*(w->var(Form("nsig%d",_count))));
	params.push_back(*(w->var(Form("cbmean%d",_count))));


	/*
	//Retrieving Parameters//

	RooRealVar* lambda =  (w->var(Form("lambda%d", _count)));
	RooRealVar* nbkg   =  (w->var(Form("nbkg%d", _count)));
	RooRealVar * mean  = (w->var(Form("mean%d", _count)));
	RooRealVar * sigma1  = (w->var(Form("sigma1%d", _count)));
	RooRealVar * sigma2  = (w->var(Form("sigma2%d", _count)));
	//	RooRealVar * c1  = (w->var(Form("c1%d", _count)));
	RooRealVar * sig1frac = (w->var(Form("sig1frac%d", _count)));




	lambda->setVal(lambda->getVal());
	nbkg->setVal(nbkg->getVal());
	mean->setVal(mean->getVal());
	sigma1->setVal(sigma1->getVal());
	sigma2->setVal(sigma2->getVal());
	//	c1->setVal(c1->getVal());
	sig1frac->setVal(sig1frac->getVal());
	*/




	/*
	// test removing bakcground
	RooRealVar* lambda =  (w->var(Form("lambda%d", _count)));
	RooRealVar* nbkg   =  (w->var(Form("nbkg%d", _count)));

	//	double step = 0.002;


	//	double nkgdValue = 51.15;
	//	double lambdaValue = -2.0236;


	cout << "BRO lambda = " << lambda->getVal() << endl;
	cout << "BRIS back = " << nbkg->getVal() << endl;

	double nkgdValue = nbkg->getVal();
	double lambdaValue = lambda->getVal();


	nbkg->setVal(50);
	lambda->setVal(lambdaValue);
	*/
	//	RooRealVar* nbkg   =  (w->var(Form("nbkg%d", _count)));
	//	nbkg->setVal(50);

	double n_signal_init = params[0].getVal();
	double n_signal_error_init = params[0].getError();

	cout << "N_signal initial value: " << n_signal_init << endl;
	cout << "N_signal initial Error value: " << n_signal_error_init << endl;



	int params_size = params.size();

	cout << "params_size " << params_size << endl;

	RooMCStudy* mcstudy = new RooMCStudy(*model, dimuon_mass,  Extended(), FitOptions(Save(kTRUE), PrintEvalErrors(0)));

	mcstudy->generateAndFit(500);

	cout << "DONE Generate and Fit " << endl;


//	TString XName[2] = {"nsig","mean"};
	TString XName[2] = {"J/#psi Signal Raw Yield","J/#psi Mean From Fit"};
	TString ParName[2] = {"nsig","mean"};
	
	vector<RooPlot*> framesPull, framesParam, framesError;

	for(int i = 0; i < params_size; ++i)
	{
		framesPull.push_back(mcstudy->plotPull(params.at(i),FrameBins(20),FrameRange(-3,3)));
		framesPull[i]->SetTitle("");
		framesParam.push_back(mcstudy->plotParam(params.at(i),FrameBins(20)));
		framesParam[i]->SetTitle("");
		framesError.push_back(mcstudy->plotError(params.at(i),FrameBins(20)));
		framesError[i]->SetTitle("");
	}

	vector<TGraph*> h1;
	vector<TGraph*> h2;
	vector<TGraph*> h3;

	for(int i = 0; i < params_size; ++i){
		h1.push_back(static_cast<TGraph*>(framesPull.at(i)->getObject(0)));
		h2.push_back(static_cast<TGraph*>(framesParam.at(i)->getObject(0)));
		h3.push_back(static_cast<TGraph*>(framesError.at(i)->getObject(0)));
	}

	gStyle->SetOptFit(0111);

	TCanvas* c_pull = new TCanvas("pulls", "pulls", 900, 800);

	
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.04);
	c_pull->SetBottomMargin(0.04);


	for(int i = 0; i < params_size; ++i){
		c_pull->cd();
		h1[i]->SetTitle("");
		h1[i]->Draw();
		c_pull->Update();
		h1[i]->Fit("gaus","","",-5,5);
		h1[i]->GetFunction("gaus")->SetLineColor(4);
		h1[i]->GetFunction("gaus")->SetLineWidth(5);
		h1[i]->GetXaxis()->SetTitle("Pull");
		h1[i]->GetYaxis()->SetTitle("");
		h1[i]->Draw("APsame");
	}


	TCanvas* c_params = new TCanvas("params", "params", 900, 800);
	for(int i = 0; i < params_size; ++i){
		c_params->cd();
		framesParam.at(i)->GetYaxis()->SetTitle("Toy MCs");
		framesParam.at(i)->GetYaxis()->SetTitleOffset(1.4);
		framesParam.at(i)->Draw();
	}


	/*TCanvas* c_params = new TCanvas("params", "params", 900, 800);

	  for(int i = 0; i < params_size; ++i){
	  c_params->cd();
	  h2[i]->SetTitle("");
	  h2[i]->Draw();
	  c_params->Update();
	  if(particle == 0){h2[i]->Fit("gaus","","",900,1200);}
	  else if(particle == 1){h2[i]->Fit("gaus","","",60, 120);}
	  h2[i]->GetFunction("gaus")->SetLineColor(4);
	  h2[i]->GetFunction("gaus")->SetLineWidth(5);
	//h2[i]->GetXaxis()->SetTitle();
	h2[i]->GetYaxis()->SetTitle("Toy MCs");
	h2[i]->Draw("same");
	}*/



	TCanvas* c_errors = new TCanvas("errors", "errors", 900, 800);
	gPad->SetLeftMargin(0.15);


	for(int i = 0; i < params_size; ++i){

		n_signal_init = params[i].getVal();
		n_signal_error_init = params[i].getError();

		c_pull->cd();
		h1[i]->SetTitle("");
		h1[i]->Draw();
		c_pull->Update();
		h1[i]->Fit("gaus","","",-5,5);
		h1[i]->GetFunction("gaus")->SetLineColor(4);
		h1[i]->GetFunction("gaus")->SetLineWidth(5);
		h1[i]->GetXaxis()->SetTitle(Form("%s Pull",XName[i].Data()));
		h1[i]->GetYaxis()->SetTitle("");
		h1[i]->Draw("APsame");
		c_errors->cd();

		h3[i]->SetTitle("");
		h3[i]->Draw();
		c_errors->Update();
		h3[i]->Fit("gaus","","",n_signal_error_init*0.7,n_signal_error_init*1.3);
		h3[i]->GetFunction("gaus")->SetLineColor(4);
		h3[i]->GetFunction("gaus")->SetLineWidth(5);
		h3[i]->GetXaxis()->SetTitle(Form("%s Error",XName[i].Data()));
		h3[i]->GetYaxis()->SetTitle("");
		h3[i]->Draw("APsame");



		c_params->cd();

		h2[i]->SetTitle("");
		h2[i]->Draw();
		c_params->Update();
		h2[i]->Fit("gaus","","",n_signal_init * 0.7,n_signal_init * 1.3);
		h2[i]->GetFunction("gaus")->SetLineColor(4);
		h2[i]->GetFunction("gaus")->SetLineWidth(5);
		h2[i]->GetXaxis()->SetTitle(Form("%s Mean",XName[i].Data()));
		h2[i]->GetYaxis()->SetTitle("");
		h2[i]->Draw("APsame");



		c_pull->SaveAs(Form("Plots/Pull/Pull_%s_%d_%d.png",ParName[i].Data(),Opt,_count));
		c_params->SaveAs(Form("Plots/Pull/Par_%s_%d_%d.png",ParName[i].Data(),Opt,_count));
		c_errors->SaveAs(Form("Plots/Pull/Error_%s_%d_%d.png",ParName[i].Data(),Opt,_count));

	}


}


