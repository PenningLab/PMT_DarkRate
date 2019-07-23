#include <TMatrixDSymEigen.h>
//#include <libRATEvent.so>
#include <TAxis.h>
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

void Fit_SPE_LG_Gaussian(TString infile){
    // Read in the file	
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
    TH1D* hspe = new TH1D("hspe","hspe",250,0.0,25);
	hspe->Sumw2();
    // pC/V = 5
    for(int i=0; i<tree->GetEntries(); i++){
        tree->GetEntry(i);
        //if (height/charge<0.06)
	if(!fillt1)
		continue;
	for(int j=0;j<npulses;j++){
            	hspe->Fill(charge[j]);
	}
    }
    //double norm = hspe->Integral(1,hspe->GetNbinsX(),"width");
    //hspe->Scale(1 / norm , "width");

	double qmin = 0;
	double qmax = 25;

    RooRealVar q("q","q",0.0,25);
    //RooRealVar q1("q1","q1",0.5,4.0);
	RooRealVar mean("mean","mean",6.5,5,9);
	RooRealVar sigma("sigma","sigma",0.4,0.0,1.3);

    RooRealVar mean1("mean1","mean1",14.0,10,15);
	RooRealVar sigma1("sigma1","sigma1",0.5,0.0,4);

    RooRealVar mean2("mean2","mean2",20,14,25);
	RooRealVar sigma2("sigma2","sigma2",1.0,0.0,15.0);


	RooGaussian gau("gau","gau",q,mean,sigma);
    RooGaussian gau1("gau1","gau1",q,mean1,sigma1);
    RooGaussian gau2("gau2","gau2",q,mean2,sigma2);
	RooRealVar ped_tau("ped_tau","ped_tau",3,0.0,8);
	RooRealVar prompt_mean("prompt_mean","prompt_mean",1, 0.5,2.5);
    RooRealVar prompt_sigma("prompt_sigma","prompt_sigma",0.1,0.0,2.0);
	RooRealVar prompt_frac("prompt_frac","prompt_frac",0.05,0.0,0.1);
    RooRealVar gau_frac("gau_frac","gau_frac",0.2,0.0,0.6);
    RooRealVar gau1_frac("gau1_frac","gau1_frac",0.1,0.0,0.6);
    RooRealVar gau2_frac("gau2_frac","gau2_frac",0.05,0.0,0.3);
	RooGaussModel gaussm("gaussm","gaussm",q,prompt_mean,prompt_sigma);
	RooDecay ped("ped","ped",q,ped_tau,gaussm,RooDecay::SingleSided);
	RooAddPdf* pdf = new  RooAddPdf("charge_pdf", "charge_pdf",
					RooArgList(ped,gau,gau1,gau2), RooArgList(prompt_frac,gau_frac,gau1_frac,gau2_frac));
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
	pdf->plotOn(frame,RooFit::Components("gau"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen));
    pdf->plotOn(frame,RooFit::Components("gau1"),RooFit::LineStyle(kDashed),RooFit::LineColor(kMagenta));
    pdf->plotOn(frame,RooFit::Components("gau2"),RooFit::LineStyle(kDashed),RooFit::LineColor(kOrange));
	//pdf->plotOn(frame,RooFit::Components("g1"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGray));
	//pdf->plotOn(frame,ROOT.RooFit.Components("prompt_peak"),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kRed))
    //    time_pdf.plotOn(frame,ROOT.RooFit.Components("triplet_peak"),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kGreen))
	//pdf.plotOn(frame,ROOT.RooFit.Components("const_bg"),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kGray))
	//g0->plotOn(frame);
	//g1->plotOn(frame);
    frame->Draw();
    c1->cd(1)->SetLogy(1);
    c1->Update();


}
