void Fit_SPE_Brandeis_doubleGamma(){
    // Input file
    TFile* fin = new TFile("PMT_SPE_long_run.root","READ");

    float charge=0;
    TTree* tree = fin->Get("pulse");
    tree->SetBranchAddress("pulseCharge",&charge);
    TH1D* hspe = new TH1D("hspe","hspe",100,0.0,10);


    for(int i=0; i<tree->GetEntries(); i++){
        tree->GetEntry(i);
        hspe->Fill(charge*5.0);//convert to pC
    }

	double qmin = 0.1;
	double qmax = 10;

    //setup fitting function
    RooRealVar q("q", "q", qmin, qmax);
	RooRealVar p0("p0", "p0", 0.3, 0.05, 0.8);
	RooRealVar k0("k0", "k0",18.0, 1.0, 30.0);
	RooRealVar q0("q0", "q0", 0.3, 0.05, 2.0);
	RooRealVar k1("k1", "k1", 2.0, 1.0, 400.0);
	RooRealVar q1("q1", "q1", 1.0, 0.8, 2.5);
	RooRealVar mu("mu", "mu", 0.0, -1.0, 1.0);
	mu.setConstant(true);
	RooGamma g0("g0", "g0", q, k0, q0, mu);
	RooGamma g1("g1", "g1", q, k1, q1, mu);
	RooRealVar ped_tau("ped_tau","ped_tau",0.1,0,1);
	RooRealVar prompt_mean("prompt_mean","prompt_mean",0.0, 0.0,1.0);
    RooRealVar prompt_sigma("prompt_sigma","prompt_sigma",0.1,0.0,1.0);
	RooRealVar prompt_frac("prompt_frac","prompt_frac",0.3,0.2,0.8);
	RooGaussModel gaussm("gaussm","gaussm",q,prompt_mean,prompt_sigma);
	RooDecay ped("ped","ped",q,ped_tau,gaussm,RooDecay::SingleSided);
	RooAddPdf* pdf = new  RooAddPdf("charge_pdf", "charge_pdf",
					RooArgList(ped,g0,g1), RooArgList(prompt_frac,p0));
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
	pdf->plotOn(frame,RooFit::Components("g0"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen));
	pdf->plotOn(frame,RooFit::Components("g1"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGray));

    frame->Draw();


}
