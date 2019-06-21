/////////////////////////////////////////////////////////////////////////////////////////
// To compile : g++ ana_combine_root.cc -o ana_combine_root -I${ROOTSYS}/include/ `root-config --cflags --libs`
// To execute (help infomation gives detail utility) : ./DDC10_data_readout -h
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
using namespace std;
static void show_usage(string name)
{
	cout << " Usage : ./ana_combine_root [-co] file1 " << name << " Options:\n"
	     << " -o : Name of output file.\n"
	     << " -i : Name of input file.\n"
	     << " -n : Number of root files\n"
	     << " -wd : Working directory\n"
	     << " -od : Output directory\n"
	     << " -t : Use temperature measurements\n"
	     << " -s : Number of sweeps\n"
	     << " -debug : Get in the debugging mode.\n"
	     << " -h or --help : Show the usage\n"
	     << " Enjoy ! -Ryan Wang" << endl;
}
int main(int argc, char* argv[])
{
	std::string filename;
	std::string outfilename;
	std::string working_dir;
	std::string out_dir = "";
	int number_of_files = 0;
	bool event_tree_enable = false;
	bool use_temp = false;
	bool trig_pmt = false;
	bool use_frac = false;
	bool use_trigger = false;
	bool calcrate = false;
	bool stored_livetime = false;
	double sample_time = 10e-9;
	int frac_time = 0;
	int frac_start = 0;
	float charge_threshold = -1;
	double num_sweeps = -1;
	std::string pydir;
	int initial_run = 0;
	if (argc < 2)
	{
		show_usage(argv[0]);
		return 1;
	}
	cout << "loading in arguments" << endl;
	for (int i = 1; i < argc; ++i)
	{
		string arg = argv[i];
		if ((arg == "-h") || (arg == "--help"))
		{
			show_usage(argv[0]);
			return 0;
		}
		else if (arg == "-wd")
		{
			working_dir = argv[i + 1];
		}
		else if (arg == "-od")
		{
			out_dir = argv[i + 1];
		}
		else if (arg == "-n")
		{
			number_of_files = atoi(argv[i + 1]);
		}
		else if (arg == "-i")
		{
			filename = argv[i + 1];
		}
		else if (arg == "-o")
		{
			outfilename = argv[i + 1];
		}
		else if (arg == "-init")
		{
			initial_run = atoi(argv[i + 1]);
		}
		else if (arg == "-e")
		{
			event_tree_enable = true;
		}
		else if (arg == "-pmt")
		{
			trig_pmt = true;
		}
		else if (arg == "-ltime")
		{
			stored_livetime = true;
		}
		else if (arg == "-frac")
		{
			use_frac = true;
			frac_time = atof(argv[i + 1]);
			frac_start = atof(argv[i + 2]);
		}
	}
	if (out_dir == "")
		out_dir = working_dir;
	std::ifstream temp_file;
	if (use_temp)
	{
		char in_temp[320];
		sprintf(in_temp, "%s/%s", working_dir.c_str(), "temptimes.txt");
		temp_file.open(in_temp, std::ifstream::in);
	}

	std::cout << "\nCreating output file" << std::endl;
	char out_path[320];
	sprintf(out_path, "%s/%s", out_dir.c_str(), outfilename.c_str());

	TFile* fout = new TFile(out_path, "RECREATE");

	//}
	std::cout << "creating output ntuple" << std::endl;
	// variables
	const int kMaxPulses = 1000;
	int number_of_samples;
	float amplitude[kMaxPulses];
	float charge_v[kMaxPulses];
	float amplitude_position[kMaxPulses];
	float pl[kMaxPulses];
	float pr[kMaxPulses];
	float biggeststep[kMaxPulses];
	// reconstruction params
	float pulse_length99[kMaxPulses];
	float pulse_length95[kMaxPulses];
	float pulse_length90[kMaxPulses];
	float pulse_length80[kMaxPulses];
	float pulse_length75[kMaxPulses];
	float pulse_length50[kMaxPulses];
	float pulse_length25[kMaxPulses];
	float pulse_length5[kMaxPulses];
	float pulse_length1[kMaxPulses];
	float pulse_length05[kMaxPulses];
	float CalibratedTime[kMaxPulses];
	// float found_pulses[kMaxPulses][kMaxPulseSamples];
	// int found_pulses_nsamples[kMaxPulses];
	// float windowratio;
	// float pulsebaseline_rms;

	float triggerHeight;
	float triggerStart;
	float triggerStartSam;
	float triggerRising1;
	float triggerRising5;
	float triggerPosition;
	float triggerWidth;

	float event_charge;
	float event_charge_ten;
	float event_baseline;
	float event_rms;
	double livetime;
	float event_rate;
	int npulses = 0;
	bool isgood;

	TTree* event = new TTree("event", "event");
	std::cout << "\nCreating output file" << std::endl;
	event->Branch("nSamples", &number_of_samples, "number_of_samples/I");
	if (stored_livetime)
	{
		event->Branch("dLiveTime_s", &livetime, "livetime/D");
		event->Branch("dEventRate_Hz", &event_rate, "event_rate/D");
	}
	event->Branch("fCharge_pC", &event_charge, "event_charge/F");
	event->Branch("fChargePrompt_pC", &event_charge_ten, "event_charge_ten/F");
	event->Branch("fBaseline_V", &event_baseline, "event_baseline/F");
	event->Branch("bIsGood", &isgood, "isgood/O");
	event->Branch("nPulses", &npulses, "npulses/I");
	event->Branch("fPulseHeight_V", amplitude, "amplitude[npulses]/F");
	event->Branch("fPulseRightEdge", pr, "pr[npulses]/F");
	event->Branch("fPulseLeftEdge", pl, "pl[npulses]/F");
	event->Branch("fPulseCharge_pC", charge_v, "charge_v[npulses]/F");
	event->Branch("fPulsePeakTime", amplitude_position, "amplitude_position[npulses]/F");
	event->Branch("fCalibratedTime", CalibratedTime, "CalibratedTime[npulses]/F");
	event->Branch("fBigStep", biggeststep, "biggeststep[npulses]/F");
	event->Branch("fPulseLength05", pulse_length05, "pulse_length05[npulses]/F");
	event->Branch("fPulseLength1", pulse_length1, "pulse_length1[npulses]/F");
	event->Branch("fPulseLength5", pulse_length5, "pulse_length5[npulses]/F");
	event->Branch("fPulseLength25", pulse_length25, "pulse_length25[npulses]/F");
	event->Branch("fPulseLength50", pulse_length50, "pulse_length50[npulses]/F");
	event->Branch("fPulseLength75", pulse_length75, "pulse_length75[npulses]/F");
	event->Branch("fPulseLength90", pulse_length90, "pulse_length90[npulses]/F");
	event->Branch("fPulseLength95", pulse_length95, "pulse_length95[npulses]/F");
	event->Branch("fPulseLength99", pulse_length99, "pulse_length99[npulses]/F");

	if (use_trigger)
	{
		event->Branch("fTriggerTime", &triggerPosition, "triggerPosition/F");
		event->Branch("fTriggerHeight_V", &triggerHeight, "triggerHeight/F");
		event->Branch("fTriggerWidth", &triggerWidth, "triggerWidth/F");
	}

	//    short int sweep=0;
	// int run=0;
	vector<double> dark_count;
	vector<double> dark_count_error;
	TH1F* h_avgphd = new TH1F("h_avgphd", "Average Pulse;Sample (10ns);ADC counts", 8192, 0, 8192);
	h_avgphd->Sumw2();

	TH1D* h_sum = new TH1D("ADC_sum_waveform", ("#font[132]{WFD SumWaveForm}"), 10000, 0, 10000);
	h_sum->SetXTitle("#font[132]{Sample (10ns)}");
	h_sum->GetXaxis()->SetLabelFont(132);
	h_sum->GetYaxis()->SetLabelFont(132);

	cout << "start looping" << endl;
	for (int i = initial_run; i < number_of_files; i++)
	{
		char root_file_name[320];
		sprintf(root_file_name, "%s/%u_%s", working_dir.c_str(), i, filename.c_str());
		cout << "Reading in file " << i << endl;
		TTree* tree;
		TFile* fin = new TFile(root_file_name, "READ");
		if (fin == NULL || fin->IsZombie())
		{
			cout << " File is corrupted ! " << endl;
			dark_count.push_back(-1);
			dark_count_error.push_back(-1);
			continue;
		}
		else
		{
			tree = (TTree*)fin->Get("event");
		}

		// double nos = fin->Get("Nsamples")->GetUniqueID();
		TH1F* dark_hit = (TH1F*)fin->Get("dark_rate");
		dark_count.push_back(dark_hit->GetMean());
		dark_count_error.push_back(dark_hit->GetMeanError());
		float waveforms[8192];
		cout << " Loading tree branches" << endl;
		tree->SetBranchAddress("nSamples", &number_of_samples);
		tree->SetBranchAddress("raw_waveforms", waveforms);
		if (stored_livetime)
			tree->SetBranchAddress("dLiveTime_s", &livetime);
		tree->SetBranchAddress("fCharge_pC", &event_charge);
		tree->SetBranchAddress("fChargePrompt_pC", &event_charge_ten);
		tree->SetBranchAddress("fBaseline_V", &event_baseline);
		tree->SetBranchAddress("bIsGood", &isgood);
		tree->SetBranchAddress("nPulses", &npulses);
		tree->SetBranchAddress("fPulseHeight_V", amplitude);
		tree->SetBranchAddress("fPulseRightEdge", pr);
		tree->SetBranchAddress("fPulseLeftEdge", pl);
		tree->SetBranchAddress("fPulseCharge_pC", charge_v);
		tree->SetBranchAddress("fPulsePeakTime", amplitude_position);
		tree->SetBranchAddress("fCalibratedTime", CalibratedTime);
		tree->SetBranchAddress("fBigStep", biggeststep);
		tree->SetBranchAddress("fPulseLength5", pulse_length5);
		tree->SetBranchAddress("fPulseLength25", pulse_length25);
		tree->SetBranchAddress("fPulseLength50", pulse_length50);
		tree->SetBranchAddress("fPulseLength75", pulse_length75);
		tree->SetBranchAddress("fPulseLength90", pulse_length90);
		tree->SetBranchAddress("fPulseLength95", pulse_length95);
		tree->SetBranchAddress("fPulseLength99", pulse_length99);

		if (use_trigger)
		{
			tree->SetBranchAddress("fTriggerTime", &triggerPosition);
			tree->SetBranchAddress("fTriggerHeight_V", &triggerHeight);
			tree->SetBranchAddress("fTriggerWidth", &triggerWidth);
		}
		// if livetime recorded and internal trigger get rate using that
		// if external trigger and external trigger time provided use nsamples to calculate rate
		// tree->SetBranchAddress("triggerpulseHeight",&triggerpulseHeight);
		// tree->SetBranchAddress("triggerpulseWidth",&triggerpulseWidth);
		// tree->SetBranchAddress("triggerpulsePeakTime",&triggerpulsePeakTime);
		cout << " processing root file No. " << i << endl;

		int nument = tree->GetEntries();
		for (int j = 0; j < nument; j++)
		{
			tree->GetEntry(j);
			// loop throught pulses
			if (stored_livetime)
				event_rate = (double)npulses / livetime;
			event->Fill();
			if (!isgood)
				continue;

			for (int sam = 0; sam < number_of_samples; sam++)
			{
				h_sum->Fill(sam, waveforms[sam]);
				// waveforms[sam] = raw_waveform[sam] - thisbase;
			}
			for (int k = 0; k < npulses; k++)
			{
				int initialsam = pl[k];
				int finalsam = pr[k];
				for (int ns = fmax(initialsam - 3, 0); ns < fmin(finalsam + 3, number_of_samples); ns++)
				{
					h_avgphd->Fill(ns, waveforms[ns]);
				}
			}
		}
		std::cout << " Finished processing file No. " << i << std::endl;
		fin->Close();
	} // main for loop
	fout->cd();

	TGraphErrors* dark_plot = new TGraphErrors();
	int backcount = 0;
	for (int h = 0; h < dark_count.size(); h++)
	{
		double temp_dark_rate = dark_count[h] * 1e8;
		double temp_dark_rate_error = dark_count_error[h] * 1e8;
		if (temp_dark_rate < 0)
		{
			backcount++;
			continue;
		}
		dark_plot->SetPoint(h - backcount, h, temp_dark_rate);
		dark_plot->SetPointError(h - backcount, 0, temp_dark_rate_error);
	}

	dark_plot->SetName("dark_plot");
	dark_plot->SetTitle(";Run (100s);Dark Rate (Hz)");
	TCanvas* cdark = new TCanvas("cdark", "cdark");
	dark_plot->SetMarkerStyle(24);
	dark_plot->SetMarkerColor(2);
	dark_plot->Draw("AP");

	h_avgphd->Write();
	dark_plot->Write();
	cdark->Write();
	event->Write();
	fout->Write();
	fout->Close();
	return 0;
}
