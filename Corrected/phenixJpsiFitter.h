#include <RooAbsReal.h>
#include <RooAddPdf.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooExponential.h>
#include <RooFitResult.h>
#include <RooGaussian.h>
#include <RooGenericPdf.h>
#include <RooHist.h>
#include <RooJohnson.h>
#include <RooMCStudy.h>
#include <RooMinuit.h>
#include <RooPlot.h>
#include <RooPolynomial.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>

#include <fstream>
#include <iomanip>
#include <string>

using namespace RooFit;
using namespace std;

float pdgMass = 3.096916;

void clean0(TH1F *h);
int drawOpt = 0;

float minhisto = 2.3;
float maxhisto = 3.8;
int nbinsmasshisto = 50;

float x_2 = 0.95;
float y_1 = 0;

TString collisionsystem = "";
Int_t _count = 0;

RooWorkspace *outputw = new RooWorkspace("w");
RooWorkspace *w_val = nullptr;

float yield = 0;
float yieldErr = 0;

float FitMean = 0;
float FitMeanErr = 0;

float FitWidth1 = 0;
float FitWidth1Err = 0;

RooAddPdf *model = nullptr;
float SignalWidth = 0.40;

TString DirName[2] = {"North", "South"};
TString MultName[3] = {"FVTXN", "FVTXS", "SVX"};

RooFitResult *fit(TString variation, TString pdf, TCanvas *c, RooDataSet *ds,
		RooDataHist *dh, RooRealVar *mass, RooPlot *&outframe, int Opt) {

	TH1 *h = dh->createHistogram("dimuon_mass");
	h->Sumw2(kFALSE);
	h->SetBinErrorOption(TH1::kPoisson);
	h->SetMarkerSize(1.55);
	h->SetMarkerStyle(20);
	h->SetLineColor(1);
	h->SetLineWidth(4);

	c->cd();
	RooPlot *frame = mass->frame();
	TPad *p1 = new TPad("p1", "p1", 0., 0.28, 1., 0.99);
	p1->SetBorderMode(1);
	p1->SetFrameBorderMode(0);
	p1->SetBorderSize(2);
	p1->SetBottomMargin(0.10);
	p1->Draw();

	float shiftYpad = 0.08;

	TPad *p2 = new TPad("p2", "p2", 0., 0.10, 1., 0.24 + shiftYpad); // 0.26
	p2->SetTopMargin(0.);
	p2->SetBorderMode(0);
	p2->SetBorderSize(2);
	p2->SetFrameBorderMode(0);
	p2->SetTicks(1, 1);
	p2->Draw();

	/*
	 * Signal Model
	 */
	p1->cd();

	RooRealVar meanMass("meanMass", "mean", pdgMass, pdgMass - 0.04,
			pdgMass + 0.2);

	RooRealVar sigma("sigma", "sigma", 0.1, 0.01, 0.2);
	RooRealVar a("a", "a", 2, 0, 12);
	RooRealVar n("n", "n", 1, 0, 200);
	RooCBShape signal(Form("signal%d", _count), "signal", *mass, meanMass, sigma,
			a, n);

	// RooGaussian signal("signal", "signal", mass, meanMass, sigma);

	RooRealVar mu = meanMass;
	RooRealVar lambda = sigma;
	RooRealVar gamma("gamma", "gamma", 1, 0, 0.1);
	RooRealVar delta("delta", "delta", 1, 0, 10);
	// RooJohnson signal("signal", "signal", mass, mu, lambda, gamma, delta);

	/*
	 * Background Model
	 */

	RooRealVar expConst("expConst", "expConst", -4, -5,
			-0.0); // creating the decay constant
	RooExponential background(
			Form("background%d", _count), "background", *mass,
			expConst); // creating the exponential background signal

	RooRealVar nSig("nSig", "nSig", 0.5 * ds->numEntries(),
			0.2 * ds->numEntries(), 2 * ds->numEntries());
	RooRealVar nBack("nBkg", "nBkg", 0.5 * ds->numEntries(),
			0.1 * ds->numEntries(), 2 * ds->numEntries());

	model =
		new RooAddPdf(Form("model%d", _count), "", RooArgList(signal, background),
				RooArgList(nSig, nBack));

	frame->SetMaximum(nSig.getVal() * 0.9);


	if(Opt == 0){
		a.setVal(2.63205);
		a.setConstant();
		n.setVal(132.520);
		n.setConstant();
	}
	if(Opt == 1){
		a.setVal(4.77014);
		a.setConstant();
		n.setVal(98.1205);
		n.setConstant();
	}

	// RooFitResult *fitResult = model->fitTo(*ds, Save(), Minos(),
	// RooFit::PrintLevel(0), Extended(kTRUE));
	RooFitResult *fitResult = model->fitTo(*ds);

	w_val->import(*model);
	w_val->import(nSig);
	w_val->import(meanMass);

	RooAbsReal *nll = model->createNLL(*ds);
	float log_likelihood = nll->getVal();

	ds->plotOn(frame, Name("ds"), Binning(nbinsmasshisto), MarkerSize(1),
			MarkerStyle(20), MarkerColor(1), LineColor(1), LineWidth(4),
			LineColor(1)); // draw an transparent hist

	model->plotOn(frame, Name(Form("bkg%d", _count)), Components(background),
			Normalization(1.0, RooAbsReal::RelativeExpected),
			Precision(1e-6), DrawOption("L"), LineStyle(7), LineColor(4),
			LineWidth(4));

	if (pdf != "1gauss") {
		model->plotOn(frame, Name(Form("sig%d", _count)), Components(signal),
				Normalization(1.0, RooAbsReal::RelativeExpected),
				Precision(1e-6), DrawOption("L"), FillStyle(3002),
				FillColor(kOrange - 3), LineStyle(7), LineColor(kOrange - 3),
				LineWidth(4));
		model->plotOn(frame, Name(Form("sigF%d", _count)), Components(signal),
				Normalization(1.0, RooAbsReal::RelativeExpected),
				Precision(1e-6), DrawOption("F"), FillStyle(3002),
				FillColor(kOrange - 3), LineStyle(7), LineColor(kOrange - 3),
				LineWidth(4));
	}

	model->plotOn(frame, Name(Form("model%d", _count)),
			Normalization(1.0, RooAbsReal::RelativeExpected),
			Precision(1e-6), DrawOption("L"), LineColor(2), LineWidth(4));


	ds->plotOn(frame, Name("ds"), Binning(nbinsmasshisto), MarkerSize(1),
			MarkerStyle(20), MarkerColor(1), LineColor(1), LineWidth(4),
			LineColor(1)); // draw an transparent hist

	frame->SetMaximum(nSig.getVal() * 0.9);
	model->paramOn(frame, Layout(x_2 + 0.5, x_2 + 0.5, y_1 + 0.16),
			Format("NEU", AutoPrecision(3)));

	frame->getAttText()->SetTextSize(0.00);
	frame->SetTitle("");
	frame->SetYTitle("Events");

	frame->GetXaxis()->CenterTitle();
	frame->GetYaxis()->CenterTitle();
	frame->GetXaxis()->SetTitleOffset(1.1);
	frame->GetYaxis()->SetTitleOffset(1.4);
	frame->GetXaxis()->SetTitleSize(0.04);
	frame->GetXaxis()->SetTitleSize(0.0);
	frame->GetYaxis()->SetTitleSize(0.05);
	frame->GetXaxis()->SetLabelFont(42);
	frame->GetYaxis()->SetTitleFont(42);
	frame->GetYaxis()->SetLabelFont(42);
	frame->GetXaxis()->SetLabelSize(0.045);
	frame->GetXaxis()->SetLabelSize(0.0);
	frame->GetYaxis()->SetLabelSize(0.05);
	frame->SetStats(0);
	frame->SetMaximum(5000);
	frame->SetMaximum(nSig.getVal() * 0.20);

	(frame->GetXaxis())->SetRangeUser(minhisto, maxhisto);
	frame->GetXaxis()->SetNdivisions(-50205);
	frame->Draw();

	RooHist *pull_hist = frame->pullHist("ds", Form("model%d", _count));

	RooPlot *pull_plot = mass->frame();
	(pull_plot->GetXaxis())->SetRangeUser(minhisto, maxhisto);

	RooRealVar x("x", "", 1e3, 0, 1e6);
	x.setVal(0.);

	RooGenericPdf *line_ref =
		new RooGenericPdf("ref_0", "ref_0", "x", RooArgList(x));

	line_ref->plotOn(pull_plot, LineStyle(7), LineColor(13), LineWidth(2));

	pull_plot->addPlotable(static_cast<RooPlotable *>(pull_hist), "XP");
	pull_plot->SetTitle("");

	pull_plot->SetXTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");

	pull_plot->SetYTitle("Pull");
	pull_plot->GetYaxis()->SetTitleFont(42);
	pull_plot->GetYaxis()->SetTitleFont(42);
	pull_plot->GetYaxis()->SetTitleSize(0.17);
	pull_plot->GetXaxis()->SetTitleOffset(0.4);

	pull_plot->GetYaxis()->CenterTitle();
	pull_plot->GetXaxis()->SetTitleSize(0.20);
	pull_plot->GetXaxis()->SetTitleOffset(1.0);
	pull_plot->GetXaxis()->CenterTitle();
	pull_plot->GetXaxis()->SetLabelFont(42);
	pull_plot->GetXaxis()->SetLabelOffset(0.01);
	pull_plot->GetXaxis()->SetLabelFont(42);
	pull_plot->GetXaxis()->SetLabelSize(0.15);
	pull_plot->GetXaxis()->SetTickLength(0.15);
	pull_plot->GetYaxis()->SetLabelFont(42);
	pull_plot->GetYaxis()->SetLabelOffset(0.01);
	pull_plot->GetYaxis()->SetLabelSize(0.15);
	pull_plot->GetYaxis()->SetNdivisions(305);
	pull_plot->GetXaxis()->SetNdivisions(-50205);

	yield = nSig.getVal();
	yieldErr = nSig.getError();

	FitMean = meanMass.getVal();
	FitMeanErr = meanMass.getError();

	FitWidth1 = sigma.getVal();
	FitWidth1Err = sigma.getError();

	TH1F *fh = (TH1F *)h->Clone("fh");

	TLegend *leg = new TLegend(0.62, 0.55, 0.89, 0.75, NULL, "brNDC");
	if (drawOpt == 1) {
		leg = new TLegend(0.60, 0.7, 0.89, 0.89, NULL, "brNDC");
	}
	if (drawOpt == 0) {
		leg =
			new TLegend(0.50 / 1.20 + 0.10, 0.42, 0.89 + 0.10, 0.86, NULL, "brNDC");
	}
	leg = new TLegend(0.60 / 1.20 + 0.10, 0.42, 0.89 + 0.10, 0.86, NULL, "brNDC");

	leg->SetBorderSize(0);
	leg->SetTextSize(0.05);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->AddEntry(h, "Data", "pl");
	leg->AddEntry(frame->findObject(Form("model%d", _count)), "Fit", "l");
	leg->AddEntry(frame->findObject(Form("sig%d", _count)), "J/#psi Signal", "f");
	leg->AddEntry(frame->findObject(Form("bkg%d", _count)), "Background", "l");

	TLatex *texcms = new TLatex(0.21, 0.88, "CMS");
	texcms->SetNDC();
	texcms->SetTextAlign(13);
	texcms->SetTextFont(62);
	texcms->SetTextSize(0.04);
	texcms->SetLineWidth(2);
	TLatex *texpre = new TLatex(0.30, 0.88, "Preliminary");
	texpre->SetNDC();
	texpre->SetTextAlign(13);
	texpre->SetTextFont(52);
	texpre->SetTextSize(0.04);
	texpre->SetLineWidth(2);
	TLatex *texsup = new TLatex(0.21, 0.815, "Supplementary");
	texsup->SetNDC();
	texsup->SetTextAlign(13);
	texsup->SetTextFont(52);
	texsup->SetTextSize(0.05);
	texsup->SetLineWidth(2);
	TLatex *texB;
	texB = new TLatex(0.21, 0.7, "J/#psi");
	texB->SetNDC();
	texB->SetTextFont(62);
	texB->SetTextSize(0.075);
	texB->SetLineWidth(2);

	TLatex *texCol;
	if (collisionsystem == "pp" || collisionsystem == "PP" ||
			collisionsystem == "ppInc")
		texCol = new TLatex(0.95, 0.94, Form("28.0 pb^{-1} (%s 5.02 TeV)", "pp"));
	else
		texCol = new TLatex(0.945, 0.94, Form("%s 5.02 TeV (300 pb^{-1})", "pp"));
	texCol->SetNDC();
	texCol->SetTextAlign(32);
	texCol->SetTextSize(0.055);
	texCol->SetTextFont(42);




	nSig.setVal(0.);
	nSig.setConstant();
	RooFitResult *fitResult_nosig = model->fitTo(*ds, Save());
	RooAbsReal *nll_nosig = model->createNLL(*ds);
	float log_likelihood_nosig = nll_nosig->getVal();
	float real_significance = sqrt(2 * (-log_likelihood + log_likelihood_nosig));
	std::cout << "REAL SIGNIFICANCE = " << real_significance << std::endl;

	leg->Draw("same");
	texB->Draw();
	c->cd();
	p1->cd();
	leg->Draw("same");
	texB->Draw();
	p2->cd();

	pull_plot->GetYaxis()->SetTitleOffset(0.4);
	pull_plot->SetYTitle("Pull");

	pull_plot->Draw();
	p1->cd();
	outframe = frame;
	outputw->import(*model);

	mass->setRange("signal", pdgMass - SignalWidth, pdgMass + SignalWidth);
	RooAbsReal *RangeBakground =
		background.createIntegral(*mass, NormSet(*mass), Range("signal"));
//	RooAbsReal *RangeSig =
//		signal.createIntegral(*mass, NormSet(*mass), Range("signal"));

//	yield = yield * RangeSig->getVal();
//	yieldErr = yieldErr * RangeSig->getVal();

	float Calback = RangeBakground->getVal() * nBack.getVal();
	float StatSig = yield / sqrt(yield + Calback);

	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.060);

	float JPsiMean = meanMass.getVal();
	float JPsiWidth = sigma.getVal();

	lat->DrawLatex(0.23, 0.47, Form("S = %.1f", yield));
	lat->DrawLatex(0.23, 0.39, Form("B = %.1f", Calback));

	lat->Draw("SAME");

	TLine *l1 = new TLine(pdgMass - SignalWidth, 0, pdgMass - SignalWidth,
			frame->GetMaximum());
	l1->SetLineStyle(2);
	l1->SetLineWidth(2);
	l1->SetLineColor(kGreen);

	TLine *l2 = new TLine(pdgMass + SignalWidth, 0, pdgMass + SignalWidth,
			frame->GetMaximum());
	l2->SetLineStyle(2);
	l2->SetLineWidth(2);
	l2->SetLineColor(kGreen);

	return fitResult;
}

void clean0(TH1F *h) {
	for (int i = 1; i <= h->GetNbinsX(); i++) {
		if (h->GetBinContent(i) == 0)
			h->SetBinError(i, 1);
	}
}

void validate_fit(RooWorkspace *w, int Arg1, int Arg2) {
	RooRealVar dimuon_mass = *(w->var("dimuon_mass"));
	RooAbsPdf *model = w->pdf(Form("model%d", _count));

	std::vector<RooRealVar> params;
	params.clear();
	params.push_back(*(w->var("nSig")));
	params.push_back(*(w->var("meanMass")));

	float n_signal_init = params[0].getVal();
	float n_signal_error_init = params[0].getError();

	int params_size = params.size();

	params_size = 2;

	RooMCStudy *mcstudy =
		new RooMCStudy(*model, dimuon_mass, Extended(),
				FitOptions(Save(kTRUE), PrintEvalErrors(0)));

	mcstudy->generateAndFit(1000);

	TString XName[2] = {"J/#psi Signal Raw Yield", "J/#psi Mean From Fit"};
	TString ParName[2] = {"nSig", "mean"};

	vector<RooPlot *> framesPull, framesParam, framesError;

	for (int i = 0; i < params_size; ++i) {
		framesPull.push_back(
				mcstudy->plotPull(params.at(i), FrameBins(20), FrameRange(-3, 3)));
		framesPull[i]->SetTitle("");
		framesParam.push_back(mcstudy->plotParam(params.at(i), FrameBins(20)));
		framesParam[i]->SetTitle("");
		framesError.push_back(mcstudy->plotError(params.at(i), FrameBins(20)));
		framesError[i]->SetTitle("");
	}

	vector<TGraph *> h1;
	vector<TGraph *> h2;
	vector<TGraph *> h3;

	for (int i = 0; i < params_size; ++i) {
		h1.push_back(static_cast<TGraph *>(framesPull.at(i)->getObject(0)));
		h2.push_back(static_cast<TGraph *>(framesParam.at(i)->getObject(0)));
		h3.push_back(static_cast<TGraph *>(framesError.at(i)->getObject(0)));
	}

	gStyle->SetOptFit(0111);

	TCanvas *c_pull = new TCanvas("pulls", "pulls", 900, 800);

	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.04);
	c_pull->SetBottomMargin(0.04);

	for (int i = 0; i < params_size; ++i) {
		c_pull->cd();
		h1[i]->SetTitle("");
		h1[i]->Draw();
		c_pull->Update();
		h1[i]->Fit("gaus", "", "", -5, 5);
		h1[i]->GetFunction("gaus")->SetLineColor(4);
		h1[i]->GetFunction("gaus")->SetLineWidth(5);
		h1[i]->GetXaxis()->SetTitle("Pull");
		h1[i]->GetYaxis()->SetTitle("");
		h1[i]->Draw("APsame");
	}

	TCanvas *c_params = new TCanvas("params", "params", 900, 800);
	for (int i = 0; i < params_size; ++i) {
		c_params->cd();
		framesParam.at(i)->GetYaxis()->SetTitle("Toy MCs");
		framesParam.at(i)->GetYaxis()->SetTitleOffset(1.4);
		framesParam.at(i)->Draw();
	}

	TCanvas *c_errors = new TCanvas("errors", "errors", 900, 800);
	gPad->SetLeftMargin(0.15);

	for (int i = 0; i < params_size; ++i) {

		n_signal_init = params[i].getVal();
		n_signal_error_init = params[i].getError();

		c_pull->cd();
		h1[i]->SetTitle("");
		h1[i]->Draw();
		c_pull->Update();
		h1[i]->Fit("gaus", "", "", -5, 5);
		h1[i]->GetFunction("gaus")->SetLineColor(4);
		h1[i]->GetFunction("gaus")->SetLineWidth(5);
		h1[i]->GetXaxis()->SetTitle(Form("%s Pull", XName[i].Data()));
		h1[i]->GetYaxis()->SetTitle("");
		h1[i]->Draw("APsame");
		c_errors->cd();

		c_params->cd();

		h2[i]->SetTitle("");
		h2[i]->Draw();
		c_params->Update();
		h2[i]->Fit("gaus", "", "", n_signal_init * 0.7, n_signal_init * 1.3);
		h2[i]->GetFunction("gaus")->SetLineColor(4);
		h2[i]->GetFunction("gaus")->SetLineWidth(5);
		h2[i]->GetXaxis()->SetTitle(Form("%s Mean", XName[i].Data()));
		h2[i]->GetYaxis()->SetTitle("");
		h2[i]->Draw("APsame");

		c_pull->SaveAs(Form("Plots/Pull/%s_%s/Pull_%s_%d.png", DirName[Arg1].Data(),
					MultName[Arg2].Data(), ParName[i].Data(), _count));
		c_params->SaveAs(Form("Plots/Pull/%s_%s/Par_%s_%d.png",
					DirName[Arg1].Data(), MultName[Arg2].Data(),
					ParName[i].Data(), _count));
		c_errors->SaveAs(Form("Plots/Pull/%s_%s/Error_%s_%d.png",
					DirName[Arg1].Data(), MultName[Arg2].Data(),
					ParName[i].Data(), _count));
	}
}
