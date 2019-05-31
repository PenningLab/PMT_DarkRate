void Fit_SPE_Gussian(){
    // Read in the file	
    TFile* fin = new TFile("/Users/wangbtc/Google Drive/Brandeis_code/PMT_SPE_long_run.root","READ");

    float charge=0;
    TTree* tree = fin->Get("pulse");
    tree->SetBranchAddress("pulseCharge",&charge);
    TH1D* hspe = new TH1D("hspe","hspe",100,0.0,10);
    // pC/V = 5
    for(int i=0; i<tree->GetEntries(); i++){
        tree->GetEntry(i);
        //if (height/charge<0.06)
            hspe->Fill(charge*5.0);
    }
    //double norm = hspe->Integral(1,hspe->GetNbinsX(),"width");
    //hspe->Scale(1 / norm , "width");

	double qmin = 0;
	double qmax = 50;

    RooRealVar q("q","q",0.0,10.0);
    //RooRealVar q1("q1","q1",0.5,4.0);
	RooRealVar mean("mean","mean",2.0,1.0,3.0);
	RooRealVar sigma("sigma","sigma",0.3,0.0,0.5);

    RooRealVar mean1("mean1","mean1",3.0,0.0,5.0);
	RooRealVar sigma1("sigma1","sigma1",0.5,0.0,0.1);

    RooRealVar mean2("mean2","mean2",6.0,3.0,8.0);
	RooRealVar sigma2("sigma2","sigma2",1.0,0.0,2.0);


	RooGaussian gau("gau","gau",q,mean,sigma);
    RooGaussian gau1("gau1","gau1",q,mean1,sigma1);
    RooGaussian gau2("gau2","gau2",q,mean2,sigma2);
	RooRealVar ped_tau("ped_tau","ped_tau",0.5,0.0,5.0);
	RooRealVar prompt_mean("prompt_mean","prompt_mean",0.0, 0.0,1.0);
    RooRealVar prompt_sigma("prompt_sigma","prompt_sigma",0.1,0.0,1.0);
	RooRealVar prompt_frac("prompt_frac","prompt_frac",0.5,0.1,1.0);
    RooRealVar gau_frac("gau_frac","gau_frac",0.4,0.0,1.0);
    RooRealVar gau1_frac("gau1_frac","gau1_frac",0.1,0.0,1.0);
    RooRealVar gau2_frac("gau2_frac","gau2_frac",0.05,0.0,1.0);
	RooGaussModel gaussm("gaussm","gaussm",q,prompt_mean,prompt_sigma);
	RooDecay ped("ped","ped",q,ped_tau,gaussm,RooDecay::SingleSided);
	RooAddPdf* pdf = new  RooAddPdf("charge_pdf", "charge_pdf",
					RooArgList(ped,gau,gau1,gau2), RooArgList(prompt_frac,gau_frac,gau1_frac,gau2_frac));
	RooDataHist* data =
	  new RooDataHist("charge_data", "charge_data", q, hspe);
	RooAbsReal* nll = pdf->createNLL(*data);

  	//RooFitResult* fitresult = pdf->fitTo(*data,RooFit::Save());
	RooFitResult* fitresult = pdf->fitTo(*data,RooFit::Save());

    TCanvas* c1 = new TCanvas("c1","c1");
    RooPlot* frame = q.frame();
    frame->SetTitle(hspe->GetName());
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
