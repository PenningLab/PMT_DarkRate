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
#include <TH3F.h>
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
double SimpsIntegral(const float samples[], double baseline, int start, int end)
{
	int len;
	double qsum = 0.0;
	if ((end - start) % 2 == 0)
	{
		/* If there are an even number of samples, then there are an odd
		number of intervals; but Simpson's rule works only on an even
		number of intervals. Therefore we use Simpson's method on the
		all but the final sample, and integrate the last interval
		using the trapezoidal rule */
		len = end - start - 1;
		qsum += (samples[end - 1] + samples[end - 2] - 2 * baseline) / 2.0;
	}
	else
		len = end - start;

	double qsimps;
	qsimps = samples[start] - baseline;
	for (int i = start; i < start + len; i += 2)
		qsimps += (samples[i] - baseline) * 4;
	for (int i = start + 1; i < len + start - 1; i += 2)
		qsimps += (samples[i] - baseline) * 2;
	qsimps += samples[start + len - 1] - baseline;
	qsimps /= 3.0;

	qsum += qsimps;
	return qsum;
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
	float basecut = 1;
	bool isS2 = false;
	bool Qdep = false;
	bool isSPE = false;
	double resistance = 0.005;
	double phd = 6.2;
	bool watch = false;
	int n;
	float thresh;
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
		else if (arg == "-bcut")
		{
			basecut = atof(argv[i + 1]);
		}
		else if (arg == "-spe")
		{
			isSPE = true;
		}
		else if (arg == "-s2")
		{
			isS2 = true;
		}
		else if (arg == "-qdep")
		{
			Qdep = true;
			n = atoi(argv[i + 1]);
			thresh = atof(argv[i + 2]);
		}
		else if (arg == "-phd")
		{
			phd = atof(argv[i + 1]);
		}
		else if (arg == "-wa")
		{
			watch = true;
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
	float pulseFWHM[kMaxPulses];
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
	float event_basecharge;
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
	event->Branch("fBaselinerms_V", &event_rms, "event_rms/F");
	event->Branch("bIsGood", &fillt, "fillt/O");
	event->Branch("fBaseQ_rms", &event_basecharge, "event_basecharge/F");
	event->Branch("nPulses", &npulses, "npulses/I");
	event->Branch("fPulseHeight_V", amplitude, "amplitude[npulses]/F");
	event->Branch("fPulseRightEdge", pr, "pr[npulses]/F");
	event->Branch("fPulseLeftEdge", pl, "pl[npulses]/F");
	event->Branch("fPulseCharge_pC", charge_v, "charge_v[npulses]/F");
	event->Branch("fPulsePeakTime", amplitude_position, "amplitude_position[npulses]/F");
	event->Branch("fCalibratedTime", CalibratedTime, "CalibratedTime[npulses]/F");
	event->Branch("fFWHM", pulseFWHM, "pulseFWHM[npulses]/F");
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
	float filterResponse, charge, pfwhm;
	filtertree->Branch("ADCCns", &filterResponse, "filterResponse/F");
	filtertree->Branch("n", &n, "n/I");
	filtertree->Branch("Charge", &charge, "charge/F");
	// filtertree->Branch("FWHM", &pfwhm, "pfwhm/F");

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
	int nmin = (isS2 ? 10 : 1);
	int nmax = (isS2 ? 100 : 17);
	int nbins = nmax - nmin + 2;
	int maxth = 1000;
	int minth = 0;
	cout << "defining hists" << endl;
	TH1F* h_avgphd = new TH1F("h_avgphd", "Average Pulse;Sample (10ns);ADC counts", 500, -250, 250);
	h_avgphd->Sumw2();
	cout << "defining pop hists" << endl;

	TH2F* h_responsepop;
	TH2F* h_responsepopno;
	TH2F* h_responseEff;
	TH2F* h_responseEffno;
	TH2F* h_responseEffComb;
	if (isSPE)
	{
		h_responsepop
		    = new TH2F("h_responsepop", "Filter peak response density (SPE);n [samples];Peak Filter Response [ADCC*samples];Fraction of Pulses",
		        nbins, nmin - 1, nmax + 1, maxth - minth, minth, maxth);
		cout << "boo" << endl;
		h_responsepopno
		    = new TH2F("h_responsepopno", "Filter peak response density (Noise);n [samples];Peak Filter Response [ADCC*samples];Fraction of Pulses",
		        nbins, nmin - 1, nmax + 1, maxth - minth, minth, maxth);
		cout << "defining eff hists" << endl;
		h_responseEff = new TH2F("h_responseEff", "Filter peak response efficiency;n [samples];Threshold [ADCC*samples];Efficiency", nbins, nmin - 1,
		    nmax + 1, maxth - minth, minth, maxth);
		h_responseEffno = new TH2F("h_responseEffno", "Filter peak response efficiency;n [samples];Threshold [ADCC*samples];Efficiency", nbins,
		    nmin - 1, nmax + 1, maxth - minth, minth, maxth);
		h_responseEffComb = new TH2F("h_responseEffComb", "Filter Combined efficiency;n [samples];Threshold [ADCC*samples];Efficiency", nbins,
		    nmin - 1, nmax + 1, maxth - minth, minth, maxth);
	}
	TH3F* h_responsepopQ;
	TH1F* h_QFrac;
	TH1F* h_responseEffQ;
	if (Qdep)
	{
		/*h_responsepopQ = new TH3F("h_responsepopQ", "Filter efficiency;n [samples];Threshold [ADCC*samples];Charge [phd];Efficiency", nbins, nmin -
		   1, nmax + 1, 600, 0, 300, maxth - minth, minth, maxth);*/
		h_QFrac = new TH1F("h_Qfrac", "Filter efficiency;Charge [phd];Efficiency", 3000, 0, 300);
		h_responseEffQ = new TH1F("h_responseEffQ", "Filter efficiency;Charge [phd];Efficiency", 3000, 0, 300);
	}

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
		tree->SetBranchAddress("fBaselinerms_V", &event_rms);
		tree->SetBranchAddress("bIsGood", &fillt);
		tree->SetBranchAddress("fBaseQ_rms", &event_basecharge);
		tree->SetBranchAddress("fWindowCharge_pC", &event_windowcharge);

		tree->SetBranchAddress("nPulses", &npulses);
		tree->SetBranchAddress("fPulseHeight_V", amplitude);
		tree->SetBranchAddress("fPulseRightEdge", pr);
		tree->SetBranchAddress("fPulseLeftEdge", pl);
		tree->SetBranchAddress("fPulseCharge_pC", charge_v);
		tree->SetBranchAddress("fPulsePeakTime", amplitude_position);
		tree->SetBranchAddress("fCalibratedTime", CalibratedTime);
		tree->SetBranchAddress("fFWHM", pulseFWHM);
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
		float A = 1, B = 0.5;
		if (isS2)
		{
			A = 0.5;
			B = 1;
		}
		double numavp = 0;
		int numQ = 0;
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
			if (!fillt && (event_basecharge) > 0.1 * phd)
			{
				// fillt = false;
				continue;
			}
			for (int k = 0; k < npulses; k++)
			{
				int initialsam = (use_trigger ? (pl[k] + triggerStartSam) : pl[k]);
				int finalsam = (use_trigger ? (pr[k] + triggerStartSam) : pr[k]);
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
				// pfwhm = pulseFWHM[k];

				numspes++;
				if (Qdep && charge_v[k] > 0.25 * phd)
				{
					charge = charge_v[k] / phd;
					// n = ;
					float tempres;
					filterResponse = 0;
					int m = n;
					if (isS2)
						m = 4 * n;
					// float descale = 1 / (8192);
					float lastval = 0;
					int ns = amplitude_position[k] + (n + m + n) / 2;
					tempres = 0;
					unsigned int endPoint = ns;
					unsigned int startPoint = ((int)ns - n + 1 > 0 ? ns - n + 1 : 0);
					for (unsigned int i = startPoint; i <= endPoint; ++i)
						filterResponse -= ((double)waveforms[i] * 8192) * B;

					endPoint = (startPoint <= 1 ? 0 : startPoint - 1);
					startPoint = ((int)endPoint - m + 1 > 0 ? endPoint - m + 1 : 0);
					for (unsigned int i = startPoint; i <= endPoint; ++i)
						filterResponse += ((double)waveforms[i] * 8192) * A;

					endPoint = (startPoint <= 1 ? 0 : startPoint - 1);
					startPoint = ((int)endPoint - n + 1 > 0 ? endPoint - n + 1 : 0);
					for (unsigned int i = startPoint; i <= endPoint; ++i)
						filterResponse -= ((double)waveforms[i] * 8192) * B;

					// filtertree->Fill();
					h_QFrac->Fill(charge_v[k] / phd);
					if (filterResponse > thresh)
						h_responseEffQ->Fill(charge_v[k] / phd);
					numQ++;
				}
			}
			charge = 1000 * SimpsIntegral(waveforms, 0, 1180, 1250) / 5.0;
			event_windowcharge = charge;

			int inits = 1180;
			int fins = 1300;

			if (watch && (charge / phd) > chargecut)
			{

				// for (n = nmin; n <= nmax; ++n)
				//{
				n = 2;
				char plotname[30];
				sprintf(plotname, "n = %d", n);
				TCanvas* c1 = new TCanvas(plotname);
				TGraph* t11 = new TGraph();
				TGraph* t22 = new TGraph();
				float tempres;

				int m = n;
				if (isS2)
					m = 4 * n;
				// float descale = 1 / (8192);
				float lastval = 0;
				// int ns = winpeak + (n + m + n) / 2;
				// tempres = 0;
				for (int ns = inits; ns <= fins; ++ns)
				{
					filterResponse = 0;
					unsigned int endPoint = ns;
					unsigned int startPoint = ((int)ns - n + 1 > 0 ? ns - n + 1 : 0);
					for (unsigned int i = startPoint; i <= endPoint; ++i)
						filterResponse -= ((double)waveforms[i] * 8192) * B;

					endPoint = (startPoint <= 1 ? 0 : startPoint - 1);
					startPoint = ((int)endPoint - m + 1 > 0 ? endPoint - m + 1 : 0);
					for (unsigned int i = startPoint; i <= endPoint; ++i)
						filterResponse += ((double)waveforms[i] * 8192) * A;

					endPoint = (startPoint <= 1 ? 0 : startPoint - 1);
					startPoint = ((int)endPoint - n + 1 > 0 ? endPoint - n + 1 : 0);
					for (unsigned int i = startPoint; i <= endPoint; ++i)
						filterResponse -= ((double)waveforms[i] * 8192) * B;

					t22->SetPoint(ns - inits, ns, filterResponse);
				}
				for (int ns = inits - 2 * n - m; ns <= fins; ++ns)
				{
					t11->SetPoint(ns - inits, ns, waveforms[ns] * 8192);
				}
				t11->SetMarkerStyle(2);
				t22->SetMarkerColor(2);
				t22->SetMarkerStyle(3);
				t22->SetMarkerSize(3);

				c1->cd(0);

				t11->Draw("alp");

				t22->Draw("p");

				c1->Draw();
				c1->Modified();
				c1->Update();
				c1->Print("stuf.png", "png");
				std::cout << "new plot" << std::endl;
				std::getchar();
				//}
			}
			if (isSPE)
			{
				int winpeak = 0;
				float winmax = 0;
				for (int wp = 1180; wp < 1250; wp++)
				{
					if (waveforms[wp] > waveforms[winpeak])
						winpeak = wp;
				}
				for (n = nmin; n <= nmax; ++n)
				{
					float tempres;
					filterResponse = 0;
					int m = n;
					if (isS2)
						m = 4 * n;
					// float descale = 1 / (8192);
					float lastval = 0;
					int ns = winpeak + (n + m + n) / 2;
					// tempres = 0;
					unsigned int endPoint = ns;
					unsigned int startPoint = ((int)ns - n + 1 > 0 ? ns - n + 1 : 0);
					for (unsigned int i = startPoint; i <= endPoint; ++i)
						filterResponse -= ((double)waveforms[i] * 8192) * B;

					endPoint = (startPoint <= 1 ? 0 : startPoint - 1);
					startPoint = ((int)endPoint - m + 1 > 0 ? endPoint - m + 1 : 0);
					for (unsigned int i = startPoint; i <= endPoint; ++i)
						filterResponse += ((double)waveforms[i] * 8192) * A;

					endPoint = (startPoint <= 1 ? 0 : startPoint - 1);
					startPoint = ((int)endPoint - n + 1 > 0 ? endPoint - n + 1 : 0);
					for (unsigned int i = startPoint; i <= endPoint; ++i)
						filterResponse -= ((double)waveforms[i] * 8192) * B;

					filtertree->Fill();
					if (charge > 0.25 * phd)
					{
						h_responsepop->Fill(n, filterResponse);
					}
					else
					{
						h_responsepopno->Fill(n, filterResponse);
					}
				}
				if (charge > 0.25 * phd)
					popscale += 1;
				else
					popscaleno += 1;
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
		if (isSPE)
		{
			for (int id = 0; id < h_responsepop->GetNbinsX(); id++)
			{
				float totEff = 0;
				float totEffNo = 0;

				int maxth = h_responsepop->GetNbinsY() + 1;
				for (int ith = maxth; ith > 0; ith--)
				{
					totEff += h_responsepop->GetBinContent(id, ith) / (double)popscale;
					totEffNo += h_responsepopno->GetBinContent(id, ith) / (double)popscaleno;
					h_responseEff->SetBinContent(id, ith, totEff);
					h_responseEffno->SetBinContent(id, ith, 1 - totEffNo);
					h_responseEffComb->SetBinContent(id, ith, sqrt((1.0 - totEffNo) * totEff));
				}
			}
			h_responsepop->Scale(1.0 / (double)popscale);
			h_responsepopno->Scale(1.0 / (double)popscaleno);

			h_responsepop->Write();
			h_responseEff->Write();
			h_responsepopno->Write();
			h_responseEffno->Write();
			h_responseEffComb->Write();

			for (int i = nmin; i <= nmax; ++i)
			{
				int lobebin = h_responseEff->GetXaxis()->FindBin(i);
				char canamae[80];
				sprintf(canamae, "(n = %d)", i);
				TCanvas* c1 = new TCanvas(canamae, canamae);
				TGraph* g1 = new TGraph();
				TGraph* g2 = new TGraph();
				for (int j = 1; j <= maxth - minth; j++)
				{
					g1->SetPoint(j - 1, h_responseEff->GetYaxis()->GetBinCenter(j), h_responseEff->GetBinContent(lobebin, j));
					g2->SetPoint(j - 1, h_responseEffno->GetYaxis()->GetBinCenter(j), h_responseEffno->GetBinContent(lobebin, j));
				}
				c1->cd();
				g1->SetTitle("Efficiency");
				g2->SetTitle("Noise Rejection Efficiency");
				g1->SetMarkerColor(2);
				g1->SetLineColor(2);
				g2->SetMarkerColor(4);
				g2->SetLineColor(4);
				g1->GetYaxis()->SetRangeUser(0.996, 1);
				g2->GetYaxis()->SetRangeUser(0.996, 1);
				g1->Draw("alp");
				g2->Draw("lp");
				fout->cd();
				c1->Write();
			}
		}
		if (Qdep)
		{
			/*
			for (int id = 0; id < h_responsepopQ->GetNbinsX(); id++)
			{

			    for (int iq = 0; iq < h_responsepopQ->GetNbinsY(); iq++)
			    {
			        float totEffQ = 0;
			        int maxth = h_responsepopQ->GetNbinsZ() + 1;
			        for (int ith = maxth; ith > 0; ith--)
			        {
			            totEffQ += h_responsepopQ->GetBinContent(id, iq, ith) / (double)numQ;
			            h_responseEffQ->SetBinContent(id, iq, ith);
			        }
			    }
			}
			h_responsepopQ->Scale(1.0 / (double)numQ);
			h_responseEffQ->Write();*/
			h_responseEffQ->Divide(h_QFrac);
			h_responseEffQ->Write();
		}
		std::cout << popscale << " SPEs in tree " << numQ << std::endl;
	}

	fout->Write();
	fout->Close();
	return 0;
}
