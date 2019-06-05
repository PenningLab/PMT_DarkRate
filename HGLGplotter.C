/////////////////////////////////////////////////////////////////////////////////////////
// To compile : g++ ana_combine_root.cc -o ana_combine_root -I${ROOTSYS}/include/ `root-config --cflags --libs` `python3.6m-config --cflags --ldflags`
// `python3.6m-config --cflags --ldflags` Requires an install of python 3.6, point to your location for
// python3.6m-config To execute (help infomation gives detail utility) : ./DDC10_data_readout -h
/* Revision log :
 *
    4/25/2018 RW : Code for reading scope's text file and do the pulse finding.
    7/23/2018 RW : Add option for user to identify each item easily.



*/
/////////////////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <dirent.h>
#include <fstream> // std::ifstream
#include <iostream>
#include <iostream> // std::cout
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <sys/stat.h> // for check directory
#include <sys/types.h>
#include <time.h>
#include <vector>

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

float shift = 0;

//#include <boost/date_time.hpp>
// using namespace std;

int HGLGplotter(TString fi1, TString fi2)
{
	std::cout << "\nCreating output file" << std::endl;
	//char out_path[320];
	//sprintf(out_path, "%s/%s", out_dir.c_str(), outfilename.c_str());

	TFile* fout = new TFile("HGLGSat.root", "RECREATE");

	//}
	//std::cout << "creating output ntuple" << std::endl;
	// variables
	const int kMaxPulses = 1000;
	int number_of_samples1;
	float amplitude1[kMaxPulses];
	float charge_v1[kMaxPulses];
	float amplitude_position1[kMaxPulses];
	float pl1[kMaxPulses];
	float pr1[kMaxPulses];

	float event_charge1;
	float event_charge_ten1;
	float event_baseline1;
	float event_rms1;
	float livetime1;
	float event_rate1;
	int npulses1 = 0;
	bool fillt1 = true;

	TChain* event1 = new TChain("Events");
	int checkker = event1->Add(fi1);
	if (checkker <= 0)
	{
		cout << " File is corrupted ! " << endl;
	}
	int nEvs1 = event1->GetEntries();
	std::cout << "\nCreating output file" << std::endl;
	event1->SetBranchAddress("nSamples", &number_of_samples1);
	event1->SetBranchAddress("fCharge_pC", &event_charge1);
	event1->SetBranchAddress("fChargePrompt_pC", &event_charge_ten1);
	event1->SetBranchAddress("fBaseline_V", &event_baseline1);
	event1->SetBranchAddress("bIsGood", &fillt1);

	event1->SetBranchAddress("nPulses", &npulses1);
	event1->SetBranchAddress("fPulseHeight_V", amplitude1);
	event1->SetBranchAddress("fPulseRightEdge", pr1);
	event1->SetBranchAddress("fPulseLeftEdge", pl1);
	event1->SetBranchAddress("fPulseCharge_pC", charge_v1);

	int number_of_samples2;
	float amplitude2[kMaxPulses];
	float charge_v2[kMaxPulses];
	float amplitude_position2[kMaxPulses];
	float pl2[kMaxPulses];
	float pr2[kMaxPulses];

	float event_charge2;
	float event_charge_ten2;
	float event_baseline2;
	float event_rms2;
	float livetime2;
	float event_rate2;
	int npulses2 = 0;
	bool fillt2 = true;

	TChain* event2 = new TChain("Events");
	checkker = event2->Add(fi2);
	if (checkker <= 0)
	{
		cout << " File is corrupted ! " << endl;
	}
	int nEvs2 = event2->GetEntries();
	std::cout << "\nCreating output file" << std::endl;
	event2->SetBranchAddress("nSamples", &number_of_samples2);
	event2->SetBranchAddress("fCharge_pC", &event_charge2);
	event2->SetBranchAddress("fChargePrompt_pC", &event_charge_ten2);
	event2->SetBranchAddress("fBaseline_V", &event_baseline2);
	event2->SetBranchAddress("bIsGood", &fillt2);

	event2->SetBranchAddress("nPulses", &npulses2);
	event2->SetBranchAddress("fPulseHeight_V", amplitude2);
	event2->SetBranchAddress("fPulseRightEdge", pr2);
	event2->SetBranchAddress("fPulseLeftEdge", pl2);
	event2->SetBranchAddress("fPulseCharge_pC", charge_v2);
	//    short int sweep=0;
	// int run=0;
	if (nEvs1 != nEvs2)
	{
		std::cout << "Files do not have same number of waveforms. exiting" << std::endl;
		return -1;
	}

	TH2F* h_hgvlg = new TH2F("h_hgvlg", "HG vs LG;HG Charge (pC);LG Charge (pC)", 1000, 0, 6500, 500, 0, 5000);
	TH2F* h_hg_qvh = new TH2F("h_hg_qvh","HG Charge vs Height;Charge (pC);Height (V)",1000,0,6500,1000,0,2.1);
	TH2F* h_lg_qvh = new TH2F("h_lg_qvh","LG Charge vs Height;Charge (pC);Height (V)",500,0,5000,1000,0,2.1);
	TH2F* h_lg_hoq = new TH2F("h_lg_hoq","LG Charge/Height",500,0,5000,500,0,2500);
	TH2F* h_hg_hoq = new TH2F("h_hg_hoq","HG Charge/Height",500,0,6500,500,0,3500);
	cout << "start looping" << endl;

	double numavp = 0;
	
	// int nument = tree->GetEntries();
	for (int j = 0; j < nEvs1; j++)
	{
		event1->GetEntry(j);
		event2->GetEntry(j);
		// loop throught pulses
		if (j % 100 == 0)
		{
			cout << " This is sweep : " << j << endl;
		}
		if (!fillt1 || !fillt2)
			continue;
		//std::vector<std::pair<int, float>> ch1;
		//std::vector<std::pair<int, float>> ch2;
		float ch1=0,ch2=0;
		int nch1=0,nch2=0;
		for (int k = 0; k < npulses1; k++)
		{
			int initialsam = pl1[k];
			int finalsam = pr1[k];

			if ((initialsam > 1220 && initialsam < 1180) || charge_v1[k]<6.3)
				continue;
			//std::pair<int, float> dum = std::make_pair(k, );
			ch1+=charge_v1[k];
			nch1++;
			h_hg_qvh->Fill(charge_v1[k],amplitude1[k]);
			h_hg_hoq->Fill(charge_v1[k],charge_v1[k]/amplitude1[k]);
			// Fill plot here
		}
		for (int k = 0; k < npulses2; k++)
		{
			int initialsam = pl2[k];
			int finalsam = pr2[k];

			if ((initialsam > 1220 && initialsam < 1180) || charge_v2[k]<0.63)
				continue;
			//std::pair<int, float> dum = std::make_pair(k, );
			ch2+=charge_v1[k];
			nch2++;
			h_lg_qvh->Fill(charge_v2[k],amplitude2[k]);
			h_lg_hoq->Fill(charge_v2[k],charge_v2[k]/amplitude2[k]);
			// Fill plot here
		}
		h_hgvlg->Fill(ch1/(double)nch1,ch2/(double)nch2);
	}
	fout->cd();
	h_hgvlg->Write();
	h_lg_qvh->Write();
	h_hg_qvh->Write();
	h_hg_hoq->Write();
	h_lg_hoq->Write();
	fout->Write();
	fout->Close();

	return 0;
}
