#include <TCanvas.h>
#include <TChain.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH2Poly.h>
#include <TKey.h>
#include <TLine.h>
#include <TMath.h>
#include <TMatrix.h>
#include <TNtuple.h>
#include <TProfile.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVector3.h>
#include <TVectorD.h>
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "Math/DistFunc.h"
#include <RooRealSumPdf.h>
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooHistPdf.h"
#include <RooAddPdf.h>
#include <RooChebychev.h>
#include <RooExponential.h>
#include <RooProdPdf.h>
#include "RooDataHist.h"
#include "RooDecay.h"
#include "RooGaussModel.h"
#include "RooGamma.h"
void Fit_SPE_LG_doubleGamma(TString infile){
    // Input file
     TFile* fin = new TFile(infile,"READ");

    const int kMaxPulses = 1000;
	int number_of_samples1;
	float amplitude1[kMaxPulses];
	float charge[kMaxPulses];
	int npulses;
	bool fillt1 = true;
	float pl[kMaxPulses];
	float pr[kMaxPulses];
    TTree* tree = (TTree*)fin->Get("Events");
    tree->SetBranchAddress("nPulses", &npulses);
    tree->SetBranchAddress("fPulseCharge_pC",charge);
    tree->SetBranchAddress("bIsGood", &fillt1);
	tree->SetBranchAddress("fPulseRightEdge", pr);
		tree->SetBranchAddress("fPulseLeftEdge", pl);
    TH1D* hspe = new TH1D("hspe","hspe",200,0.0,30);
	hspe->Sumw2();
    // pC/V = 5
    for(int i=0; i<tree->GetEntries(); i++){
        tree->GetEntry(i);
        //if (height/charge<0.06)
	if(!fillt1)
		continue;
	for(int j=0;j<npulses;j++){
	    if(pl[j]<1200 && pl[j]>1180)
            	hspe->Fill(charge[j]);
	}
    }

	double qmin = 0.1;
	double qmax = 30;

    //setup fitting function
    RooRealVar q("q", "q", qmin, qmax);
	RooRealVar p0("p0", "p0", 0.4, 0.1, 0.9);
	RooRealVar k0("k0", "k0",18.0, 0.8, 25.0);
	RooRealVar q0("q0", "q0", 0.3, 0.01, 1.5);
	RooRealVar k1("k1", "k1", 2.0, 1.0, 300.0);
	RooRealVar q1("q1", "q1", 1.0, 0.6, 3);
	RooRealVar mu("mu", "mu", 0.0, -3.0, 3.0);
	mu.setConstant(true);
	RooGamma g0("g0", "g0", q, k0, q0, mu);
	RooGamma g1("g1", "g1", q, k1, q1, mu);
	RooRealVar ped_tau("ped_tau","ped_tau",0.1,0,2.2);
	RooRealVar prompt_mean("prompt_mean","prompt_mean",1.1, 0.1,1.6);
    RooRealVar prompt_sigma("prompt_sigma","prompt_sigma",0.1,0.0,0.6);
	RooRealVar prompt_frac("prompt_frac","prompt_frac",0.1,0.01,0.4);
	RooGaussModel gaussm("gaussm","gaussm",q,prompt_mean,prompt_sigma);
	RooDecay ped("ped","ped",q,ped_tau,gaussm,RooDecay::SingleSided);
	RooAddPdf* pdf = new  RooAddPdf("charge_pdf", "charge_pdf",
					RooArgList(ped,g0,g1), RooArgList(prompt_frac,p0));
	RooDataHist* data =
	  new RooDataHist("charge_data", "charge_data", q, hspe);
	RooAbsReal* nll = pdf->createNLL(*data);

  	//RooFitResult* fitresult = pdf->fitTo(*data,RooFit::Save());
	RooFitResult* fitresult = pdf->chi2FitTo(*data,RooFit::Save(),RooFit::SumW2Error(kTRUE));

    TCanvas* c1 = new TCanvas("c1","c1");
    RooPlot* frame = q.frame();
    frame->SetTitle("LG");
    data->plotOn(frame);
    pdf->plotOn(frame,RooFit::LineColor(kBlue));
    pdf->plotOn(frame,RooFit::Components("ped"),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
	pdf->plotOn(frame,RooFit::Components("g0"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen));
	pdf->plotOn(frame,RooFit::Components("g1"),RooFit::LineStyle(kDashed),RooFit::LineColor(kOrange));

    frame->Draw();


}
