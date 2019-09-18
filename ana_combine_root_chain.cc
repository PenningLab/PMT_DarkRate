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
	float chargecut = 0;
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
		else if (arg == "-ccut")
		{
			chargecut = atof(argv[i + 1]);
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
	float triggerRising05;
	float triggerPosition;
	float triggerWidth;

	float event_charge;
	float event_windowcharge;
	float event_charge_ten;
	float event_baseline;
	float event_rms;
	float livetime;
	float event_rate;
	int npulses = 0;
	bool fillt = true;

	TTree* event = new TTree("Events", "Events");
	std::cout << "\nCreating output file" << std::endl;
	event->Branch("nSamples", &number_of_samples, "number_of_samples/I");
	if (stored_livetime)
	{
		event->Branch("dLiveTime_s", &livetime, "livetime/D");
		event->Branch("dEventRate_Hz", &event_rate, "event_rate/D");
	}
	event->Branch("fCharge_pC", &event_charge, "event_charge/F");
	event->Branch("fWindowCharge_pC", &event_windowcharge, "event_windowCharge/F");
	event->Branch("fChargePrompt_pC", &event_charge_ten, "event_charge_ten/F");
	event->Branch("fBaseline_V", &event_baseline, "event_baseline/F");
	event->Branch("bIsGood", &fillt, "fillt/O");
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

	TTree* filtertree = new TTree("filter", "filterTree");
	float filterResponse, charge;
	int n;
	filtertree->Branch("ADCCns", &filterResponse, "filterResponse/F");
	filtertree->Branch("n", &n, "n/I");
	filtertree->Branch("Charge", &charge, "charge/F");

	if (use_trigger)
	{
		// event->Branch("nTriggers", &nTrigs, "nTriggers/I");
		event->Branch("fTriggerStart", &triggerStart, "triggerStart/F");
		event->Branch("fTriggerStartSam", &triggerStartSam, "triggerStartSam/F");
		event->Branch("fTriggerRising05", &triggerRising05, "triggerRising05/F");
		event->Branch("fTriggerRising1", &triggerRising1, "triggerRising1/F");
		event->Branch("fTriggerRising5", &triggerRising5, "triggerRising5/F");
		event->Branch("fTriggerTime", &triggerPosition, "triggerPosition/F");
		event->Branch("fTriggerHeight_V", &triggerHeight, "triggerHeight/F");
		event->Branch("fTriggerWidth", &triggerWidth, "triggerWidth/F");
	}

	//    short int sweep=0;
	// int run=0;
	vector<double> dark_count;
	vector<double> dark_count_error;
	TH1F* h_avgphd = new TH1F("h_avgphd", "Average Pulse;Sample (10ns);ADC counts", 500, -250, 250);
	h_avgphd->Sumw2();
	TH2F* h_responsepop = new TH2F("h_responsepop",
	    "Filter peak response density (SPE);n [samples];Peak Filter Response [ADCC*samples];Fraction of Pulses", 20, 0, 20, 2000, 0, 2000);
	TH2F* h_responsepopno = new TH2F("h_responsepopno",
	    "Filter peak response density (Noise);n [samples];Peak Filter Response [ADCC*samples];Fraction of Pulses", 20, 0, 20, 2000, 0, 2000);
	TH2F* h_responseEff
	    = new TH2F("h_responseEff", "Filter peak response efficiency;n [samples];Threshold [ADCC*samples];Efficiency", 20, 0, 20, 2000, 0, 2000);
	TH2F* h_responseEffno
	    = new TH2F("h_responseEffno", "Filter peak response efficiency;n [samples];Threshold [ADCC*samples];Efficiency", 20, 0, 20, 2000, 0, 2000);
	TH1D* h_sum = new TH1D("ADC_sum_waveform", ("#font[132]{WFD SumWaveForm}"), 10000, 0, 10000);
	h_sum->SetXTitle("#font[132]{Sample (10ns)}");
	h_sum->GetXaxis()->SetLabelFont(132);
	h_sum->GetYaxis()->SetLabelFont(132);

	cout << "start looping" << endl;
	// for (int i = initial_run; i < number_of_files; i++)
	//{
	char root_file_name[320];
	sprintf(root_file_name, "%s/%s", working_dir.c_str(), filename.c_str());
	cout << "Reading in files " << root_file_name << endl;
	TChain* tree = new TChain("event");
	int checkker = tree->Add(root_file_name);
	int numspes = 0, popscale = 0, popscaleno = 0;
	// TFile* fin = new TFile(root_file_name, "READ");
	if (checkker <= 0 || tree->GetEntries() <= 0)
	{
		cout << " File is corrupted ! " << endl;
	}
	else
	{
		// double nos = fin->Get("Nsamples")->GetUniqueID();
		// TH1F* dark_hit = (TH1F*)fin->Get("dark_hits");
		// dark_count.push_back(dark_hit->GetMean());
		// dark_count_error.push_back(dark_hit->GetMeanError());
		float waveforms[8192];
		cout << " Loading tree branches" << endl;
		tree->SetBranchAddress("nSamples", &number_of_samples);
		tree->SetBranchAddress("raw_waveforms", waveforms);
		if (stored_livetime)
			tree->SetBranchAddress("dLiveTime_s", &livetime);
		tree->SetBranchAddress("fCharge_pC", &event_charge);
		tree->SetBranchAddress("fChargePrompt_pC", &event_charge_ten);
		tree->SetBranchAddress("fBaseline_V", &event_baseline);
		tree->SetBranchAddress("bIsGood", &fillt);
		tree->SetBranchAddress("fWindowCharge_pC", &event_windowcharge);

		tree->SetBranchAddress("nPulses", &npulses);
		tree->SetBranchAddress("fPulseHeight_V", amplitude);
		tree->SetBranchAddress("fPulseRightEdge", pr);
		tree->SetBranchAddress("fPulseLeftEdge", pl);
		tree->SetBranchAddress("fPulseCharge_pC", charge_v);
		tree->SetBranchAddress("fPulsePeakTime", amplitude_position);
		tree->SetBranchAddress("fCalibratedTime", CalibratedTime);
		tree->SetBranchAddress("fBigStep", biggeststep);
		tree->SetBranchAddress("fPulseLength05", pulse_length05);
		tree->SetBranchAddress("fPulseLength1", pulse_length1);
		tree->SetBranchAddress("fPulseLength5", pulse_length5);
		tree->SetBranchAddress("fPulseLength25", pulse_length25);
		tree->SetBranchAddress("fPulseLength50", pulse_length50);
		tree->SetBranchAddress("fPulseLength75", pulse_length75);
		tree->SetBranchAddress("fPulseLength90", pulse_length90);
		tree->SetBranchAddress("fPulseLength95", pulse_length95);
		tree->SetBranchAddress("fPulseLength99", pulse_length99);

		if (use_trigger)
		{
			// event->Branch("nTriggers", &nTrigs, "nTriggers/I");
			tree->SetBranchAddress("fTriggerStart", &triggerStart);
			tree->SetBranchAddress("fTriggerStartSam", &triggerStartSam);
			tree->SetBranchAddress("fTriggerRising05", &triggerRising05);
			tree->SetBranchAddress("fTriggerRising1", &triggerRising1);
			tree->SetBranchAddress("fTriggerRising5", &triggerRising5);
			tree->SetBranchAddress("fTriggerTime", &triggerPosition);
			tree->SetBranchAddress("fTriggerHeight_V", &triggerHeight);
			tree->SetBranchAddress("fTriggerWidth", &triggerWidth);
		}
		// if livetime recorded and internal trigger get rate using that
		// if external trigger and external trigger time provided use nsamples to calculate rate
		// tree->SetBranchAddress("triggerpulseHeight",&triggerpulseHeight);
		// tree->SetBranchAddress("triggerpulseWidth",&triggerpulseWidth);
		// tree->SetBranchAddress("triggerpulsePeakTime",&triggerpulsePeakTime);
		// cout << " processing root file No. " << i << endl;
		double numavp = 0;
		// flaot plotscale = 0;
		int nument = tree->GetEntries();
		for (int j = 0; j < nument; j++)
		{
			tree->GetEntry(j);
			// loop throught pulses
			// fillt = true;
			if (j % 100 == 0)
			{
				cout << " This is sweep : " << j << ", filter = " << filterResponse << endl;
			}
			if (stored_livetime)
				event_rate = (double)npulses / livetime;

			for (int sam = 0; sam < number_of_samples; sam++)
			{
				h_sum->Fill(sam, waveforms[sam]);
				// waveforms[sam] = raw_waveform[sam] - thisbase;
			}
			int passedpulses = 0;
			for (int k = 0; k < npulses; k++)
			{
				int initialsam = (use_trigger ? (pl[k] + triggerStartSam) : pl[k]);
				int finalsam = (use_trigger ? (pr[k] + triggerStartSam) : pr[k]);
				if (fillt)
				{
					// fillt = false;
					break;
				}
				if (amplitude_position[k] > 1240 && amplitude_position[k] < 1150)
					continue;
				for (int ns = fmax(initialsam - 20, 0); ns < fmin(finalsam + 20, number_of_samples); ns++)
				{
					h_avgphd->Fill(ns - amplitude_position[k], waveforms[ns]);
				}
				if (k > passedpulses)
				{
					amplitude[passedpulses] = amplitude[k];
					charge_v[passedpulses] = charge_v[k];
					amplitude_position[passedpulses] = amplitude_position[k];
					pl[passedpulses] = pl[k];
					pr[passedpulses] = pr[k];
					biggeststep[passedpulses] = biggeststep[k];
					// reconstruction params
					pulse_length99[passedpulses] = pulse_length99[k];
					pulse_length95[passedpulses] = pulse_length95[k];
					pulse_length90[passedpulses] = pulse_length90[k];
					// pulse_length80[passedpulses] = pulse_length80[k];
					pulse_length75[passedpulses] = pulse_length75[k];
					pulse_length50[passedpulses] = pulse_length50[k];
					pulse_length25[passedpulses] = pulse_length25[k];
					pulse_length5[passedpulses] = pulse_length5[k];
					pulse_length1[passedpulses] = pulse_length1[k];
					pulse_length05[passedpulses] = pulse_length05[k];
					CalibratedTime[passedpulses] = CalibratedTime[k];
				}
				numavp++;
				passedpulses++;

				charge = charge_v[k];
				if (charge > chargecut)
					popscale += 1;
				if (charge <= chargecut)
					popscaleno += 1;
				for (n = 1; n < 17; ++n)
				{
					float tempres;
					filterResponse = 0;
					// float descale = 1 / (8192);
					for (int ns = fmax(initialsam, 0); ns < fmin(finalsam, number_of_samples); ns++)
					{
						tempres = 0;
						unsigned int endPoint = ns;
						unsigned int startPoint = ((int)ns - n + 1 > 0 ? ns - n + 1 : 0);
						for (unsigned int i = startPoint; i <= endPoint; ++i)
							tempres -= ((double)waveforms[i] * 8192) * 0.5;

						endPoint = (startPoint <= 1 ? 0 : startPoint - 1);
						startPoint = ((int)endPoint - n + 1 > 0 ? endPoint - n + 1 : 0);
						for (unsigned int i = startPoint; i <= endPoint; ++i)
							tempres += ((double)waveforms[i] * 8192) * 1;

						endPoint = (startPoint <= 1 ? 0 : startPoint - 1);
						startPoint = ((int)endPoint - n + 1 > 0 ? endPoint - n + 1 : 0);
						for (unsigned int i = startPoint; i <= endPoint; ++i)
							tempres -= ((double)waveforms[i] * 8192) * 0.5;
						if (tempres > filterResponse)
							filterResponse = tempres;
					}
					filtertree->Fill();
					if (charge > chargecut)
					{
						h_responsepop->Fill(n, filterResponse);
					}

					if (charge <= chargecut)
					{
						h_responsepopno->Fill(n, filterResponse);
					}
				}
				numspes++;
			}
			if (fillt)
				npulses = passedpulses;
			event->Fill();
		}
		// fin->Close();
		//} // main for loop
		fout->cd();
		h_avgphd->Scale(1.0 / numavp);
		h_avgphd->Write();
		event->Write();
		filtertree->Write();
		h_responsepop->Scale(1.0 / (double)popscale);
		h_responsepopno->Scale(1.0 / (double)popscaleno);
		for (int id = 0; id < h_responsepop->GetNbinsX(); id++)
		{
			float totEff = 0;
			float totEffNo = 0;
			int maxth = h_responsepop->GetNbinsY();
			for (int ith = maxth; ith > 0; ith--)
			{
				totEff += h_responsepop->GetBinContent(id, ith);
				h_responseEff->SetBinContent(id, ith, totEff);
				totEffNo += h_responsepopno->GetBinContent(id, ith);
				h_responseEffno->SetBinContent(id, ith, totEffNo);
			}
		}
		std::cout << numspes << " SPEs in tree" << std::endl;
		h_responsepop->Write();
		h_responseEff->Write();
		h_responsepopno->Write();
		h_responseEffno->Write();
	}

	fout->Write();
	fout->Close();
	return 0;
}
