#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>   // for check directory
#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include <time.h>
#include <numeric>
//#include <libRATEvent.so>
#include <TVector3.h>
#include <TString.h>
#include <TLine.h>
#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TH1F.h>
#include <TAxis.h>
#include <TKey.h>
#include <TF1.h>
#include <TChain.h>
#include <TProfile.h>
#include <TSystem.h>
#include <TH2F.h>
#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TVectorD.h>
#include <TMatrix.h>
#include <TH2Poly.h>
using namespace std;

void PlotSweep(TString filename,int sweep){
	TFile* infile = new TFile(filename.Data(),"READ");
	TTree* wforms = (TTree*)infile->Get("waveforms");
	int NEnts = wforms->GetEntries();

	if(sweep>=NEnts)
		cout<<"Sweep number too high. Only "<<NEnts<<" entries in file: "<<filename.Data()<<endl;
	else{
		int nos = infile->Get("Nsamples")->GetUniqueID();
		TGraph* t1 = new TGraph();
		vector<float> wformdata(nos);
		wforms->SetBranchAddress("pmt_waveforms",wformdata.data());
		wforms->GetEntry(sweep);
		for(int i=0;i<nos;i++){
			t1->SetPoint(i,i,wformdata[i]);
		}
		t1->SetTitle(";Sample (10ns);V");
		TCanvas * c1 = new TCanvas("sweep","sweep");
		t1->Draw("alp");
	}
}
