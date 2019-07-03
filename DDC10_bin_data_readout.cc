/////////////////////////////////////////////////////////////////////////////////////////
// To compile : g++ -o DDC10_bin_data_readout DDC10_bin_data_readout.cc -I${ROOTSYS}/include/ `root-config --cflags --libs`
// To execute (help infomation gives detail utility) : ./DDC10_data_readout -h
/* Revision log :
 *
        4/25/2018 RW : Code for reading scope's text file and do the pulse
 finding. 7/23/2018 RW : Add option for user to identify each item easily.



*/
/////////////////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>
#include <dirent.h>
#include <fstream> // std::ifstream
#include <iostream> // std::cout
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <sys/stat.h> // for check directory
#include <sys/types.h>
#include <time.h>
#include <vector>

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
#include <TObject.h>
#include <TProfile.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVector3.h>
#include <TVectorD.h>
using namespace std;
ifstream fin;
bool signal_start = false;
bool debug_mode = false;

bool lim_sams = false;

bool smoothing = false;
bool rolling = false;
bool triangle = false;
bool boxsmoothing = false;
bool use_basefile = false;

bool isgood = true;

// define pulse finding parameters
double pulseThresh = 5.0;
double trigPulseThresh = 3.5;
double windowSize = 3.0;
double edgeThresh = 3.0;
double lookforward = 3.0;
int current_sweep = 0;
double number_of_peaks = 0.0;
int adc_per_Volt = 8192;
double resistance = 0.005;
int baseline_samples_set = 160;
int promptwindow = 300;
// smoothing parameters
int MovingWindowSize;
int iteration = 0.0;
int pth = 3.5;
int windowstart = 18;
int windowfin = 30;

// input parameters
int number_of_samples;
int Nchannels;
int Nevts;
short int* buff;

// internal run tracking
vector<float> raw_waveform;

vector<float> startv;
vector<float> endv;
vector<float> start_pointv;
vector<float> end_pointv;
vector<float> baseline;
vector<float> baselinev;
vector<float> trigbaselinev;
vector<float> pulse_left_edge;
vector<float> pulse_right_edge;

vector<float> smoothingv;
vector<float> smoothedBaseLine;

// output parameters
const int kMaxPulses = 1000;

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
float CalibratedTime[kMaxPulses];
// float windowratio;
// float pulsebaseline_rms;

float triggerHeight = 0;
float triggerPosition;
float triggerWidth;

float event_charge;
float event_charge_ten;
float event_baseline;
float event_rms;
float event_windowCharge;
int npulses = 0;
vector<double> vlivetime;

// Simpson Integral
double SimpsIntegral(const vector<float>& samples, double baseline, int start, int end)
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

// Pulse Finding
void extract_event(vector<float>& v, double b, double rms, int nos, int trigger = 0, bool trig = false)
{

	double pThresh = (trig ? trigPulseThresh : pulseThresh) * rms * windowSize;
	double eThresh = edgeThresh * rms * windowSize;
	isgood=true;
	event_windowCharge =0;
	if (trig)
		triggerHeight = 0;
	double temp_charge = 0;
	double temp_ten_charge = 0;
	double pulse_height_thresh = pth * rms;
	// cout<<" vector size is : "<<v.size()<<endl;
	// getchar();
	// Let's looking for the Pulses
	for (int i = 0; i < nos - windowSize; i++)
	{
		// std::cout << "Sample " << i << std::endl;
		double integral = SimpsIntegral(v, b, i, i + windowSize);
		int left = 0;
		int right = 0;
		int temp_peak = 0;
		double temp_startv = 0;
		double temp_endv = 0;

		double temp_bigstep = 0;

		if ((v[i] - b) > pulse_height_thresh)
		{
			if (debug_mode)
			{
				cout << " This is sample : " << i << "; integral value is : " << integral << "; pThresh is : " << pThresh << endl;
			}
			left = i;
			integral = 1.0e9;
			while ((integral) > eThresh && left > windowSize)
			{
				left--;
				integral = SimpsIntegral(v, b, left, left + windowSize);
				temp_startv = v[left];
				if (debug_mode)
				{

					cout << "Left is  : " << left << " integral is : " << integral << "; eThresh is : " << eThresh << endl;
					getchar();
				}
			}

			if (debug_mode)
				cout << "Starting right edge search" << endl;
			integral = 1.0e9;
			right = i + windowSize;
			double thischarge = 0;
			bool end = false;
			while (!end)
			{
				while ((integral) > eThresh && right < nos - 1)
				{
					right++;
					integral = SimpsIntegral(v, b, right - windowSize, right);
					temp_endv = v[right];
					if (debug_mode)
					{

						cout << "right window is  : " << right << "; integral is : " << integral << "; eThresh is : " << eThresh << endl;
						getchar();
					}
				}

				end = true;
				int r = right;
				while (r < fmin(fmin((int)v.size() - 1, nos - 1), right + lookforward))
				{
					r++;
					integral = SimpsIntegral(v, b, r - windowSize, r);
					if ((integral) > pThresh)
					{
						right = r;
						end = false;
						if (debug_mode)
						{

							cout << "right lookforward is  : " << right << "; integral is : " << integral << "; eThresh is : " << pThresh << endl;
							getchar();
						}
						break;
					}
				}
			}
			double max = -1.0e9;
			double totalq = SimpsIntegral(v, b, left, right);
			double tempq5 = 0, tempq25 = 0, tempq50 = 0, tempq75 = 0, tempq90 = 0, tempq95 = 0, tempq99 = 0;
			for (int j = left; j < right; j++)
			{
				double s = v[j] - b;
				if ((s) > max)
				{
					max = s;
					temp_peak = j;
					if (j > 0 && (v[j] - v[j - 1]) > temp_bigstep)
						temp_bigstep = v[j] - v[j - 1];
				}
				double ratio = SimpsIntegral(v, b, left, j) / totalq;
				if (ratio <= 0.05 && tempq5 < (j - left))
					tempq5 = j - left;
				if (ratio <= 0.25 && tempq25 < (j - left))
					tempq25 = j - left;
				if (ratio <= 0.5 && tempq50 < (j - left))
					tempq50 = j - left;
				if (ratio <= 0.75 && tempq75 < (j - left))
					tempq75 = j - left;
				if (ratio <= 0.9 && tempq90 < (j - left))
					tempq90 = j - left;
				if (ratio <= 0.95 && tempq95 < (j - left))
					tempq95 = j - left;
				if (ratio <= 0.99 && tempq99 < (j - left))
					tempq99 = j - left;
			}
			if (right > nos)
				continue;

			/*// noise veto
			float width = (right - left);
			float nwidth = 2 * width;
			int nright = right + nwidth;
			int nleft = left - nwidth;
			if (nright > nos)
			{
			    nright = nos;
			    nleft -= (nleft > nwidth ? nwidth : 0);
			}
			if (nleft < 0)
			{
			    nleft = 0;
			    nright += (nright < (nos - nwidth) ? nwidth : nos);
			}
			double dratio
			    = SimpsIntegral(v, b, nleft, nright) / (((nwidth + 1) / width) * SimpsIntegral(v, b, left, right));
			*/
			// cout<<" Peak is : "<<temp_peak<<" max is : "<<max<<endl;

			// if (temp_peak>0 &&temp_peak<8000){
			// if (thischarge<1.0)
			// if (trig)
			//	break;

			//}
			// cout<<" This is sample : "<<i<<" Charge integral is :
			// "<<SimpsIntegral(v,b,left,right)<<endl;

			if (SimpsIntegral(v, b, left, right) <= eThresh)
			{
				i = right - 1;
				continue;
			}
			i = right;

			// if (max < pulse_height_thresh)
			//	continue;

			startv.push_back(temp_startv);
			endv.push_back(temp_endv);
			if (signal_start)
			{
				npulses++;
				amplitude[npulses - 1] = max;
				amplitude_position[npulses - 1] = temp_peak;
				pl[npulses - 1] = left;
				pr[npulses - 1] = right;
				charge_v[npulses - 1] = SimpsIntegral(v, b, left, right) / resistance;

				CalibratedTime[npulses - 1] = temp_peak - trigger;
				// windowratio[npulses - 1] = dratio;
				// pulsebaseline_rms[npulses - 1] = rms;
				// event_n[npulses-1]=current_sweep;
				biggeststep[npulses - 1] = temp_bigstep;
				pulse_length5[npulses - 1] = tempq5;
				pulse_length25[npulses - 1] = tempq25;
				pulse_length50[npulses - 1] = tempq50;
				pulse_length75[npulses - 1] = tempq75;
				// pulse_length80[npulses-1]=tempq80;
				pulse_length90[npulses - 1] = tempq90;
				pulse_length95[npulses - 1] = tempq95;
				pulse_length99[npulses - 1] = tempq99;
				if (left < baseline_samples_set)
					isgood = false;
				temp_charge += charge_v[npulses - 1];
				if (i < promptwindow)
					temp_ten_charge += charge_v[npulses - 1];
				pulse_left_edge.push_back(left);
				pulse_right_edge.push_back(right);
			}
			else if (triggerHeight < max)
			{
				triggerHeight = max;
				triggerPosition = left;
				triggerWidth = right - left;
				pulse_left_edge.push_back(left);
				pulse_right_edge.push_back(right);
			}

		} // if statement
	}
	// getchar();
	if (!trig)
	{
		event_charge_ten = temp_ten_charge;
		event_charge = temp_charge;
		event_windowCharge = SimpsIntegral(v, b, trigger + windowstart, trigger + windowfin) / resistance;
		//if(isgood)
		//	std::cout<<event_windowCharge<<std::endl;
	}
}
// Find the baseline
double baseline_rms(vector<float>& v, vector<float>& sample, double* irms)
{
	double rms = 0;
	double temp_base = 0;
	// double baseline_samples = accumulate(v.begin(),v.end(),0);
	double baseline_samples = 0;
	for (int k = 0; k < v.size(); k++)
	{
		baseline_samples += v[k];
	}
	// if (signal_start&&debug_mode)
	//    cout<<" baseline_samples is  : "<<baseline_samples<<endl;
	baseline_samples /= v.size();
	for (int i = 0; i < v.size(); i++)
	{
		//    if (signal_start&&debug_mode)
		//        cout<<" baseline sample is  : "<<v[i]<<endl;
		rms += pow(v[i] - baseline_samples, 2);
	}
	rms = sqrt(rms / v.size());
	if (signal_start && debug_mode)
	{
		cout << " rms is : " << rms << " ;baseline is : " << baseline_samples << endl;
		getchar();
	}

	//}
	irms[0] = rms;
	return baseline_samples;
}

// Trigger analysis
double Trigger_info(vector<float> waveform)
{
	double rms_trigger = 0;
	double base_trigger = baseline_rms(trigbaselinev, waveform,
	    &rms_trigger); // Calculate baseline and rms then pass to pulse finder
	extract_event(waveform, base_trigger, rms_trigger, number_of_samples, 0, true);
	double time;
	baselinev.clear();
	if (pulse_left_edge.size() == 0)
	{
		time = 0;
		cout << " Can not find the trigger pulse for this event ! " << endl;
	}
	else
	{
		time = (pulse_left_edge.back());
	}
	pulse_left_edge.clear();
	pulse_right_edge.clear();
	trigbaselinev.clear();
	// waveform.clear();
	return time;
}

void getwaveform(vector<float>& v, int channel, int numread, float mult = 1, bool trig = false)
{
	int starti = ((current_sweep - numread) * Nchannels + channel) * (4 + 2 + number_of_samples) + 4;
	/*if(current_sweep%100000==0){
	        cout<<"starting at "<<starti<<" in buffer"<<endl;
	        cout<<"Evt "<<(current_sweep-numread)<<" of this buffer"<<endl;
	}*/
	double datum;
	for (int i = 0; i < number_of_samples; i++)
	{
		datum = (double)buff[i + starti] * mult / (double)adc_per_Volt;
		v.push_back(datum);
		if (i < baseline_samples_set)
		{
			if (!trig)
				baselinev.push_back(datum);
			else
				trigbaselinev.push_back(datum);
		}
	}
}
int calcnumchannels(int mask)
{
	int numchans = 0;
	while (mask > 0)
	{
		numchans += mask % 2;
		mask /= 2;
	}
	return numchans;
}
bool readlivetime(char datafilename[])
{
	char linea[250];
	cout << "Attempting to read log file" << endl;
	TString filenameform = datafilename;
	filenameform.ReplaceAll(".bin", ".log");
	ifstream loginfile;
	loginfile.open(filenameform.Data(), ios::in);
	if (!loginfile.is_open())
		return false;
	else
	{
		// first five dummy lines
		loginfile.getline(linea, 250);
		loginfile.getline(linea, 250);
		loginfile.getline(linea, 250);
		loginfile.getline(linea, 250);
		loginfile.getline(linea, 250);
	}

	vlivetime.resize(Nevts);
	int evtcounter = 0;
	while (loginfile.getline(linea, 250) && evtcounter < Nevts)
	{
		evtcounter++;
		std::stringstream blank(linea);
		int iev = -1, dum1 = -1;
		long int liveentry = -1;
		char c = ',';
		// bool isfail = false;
		// cout << linea << endl;
		std::string item;
		if (std::getline(blank, item, c))
			iev = std::atoi(item.data());
		if (std::getline(blank, item, c))
			dum1 = std::atoi(item.data());
		if (std::getline(blank, item, c))
			liveentry = std::atol(item.data());

		vlivetime[Nevts - evtcounter] = (double)liveentry / 6e8;
		// cout << iev << ", live entry is: " << liveentry << ", livetime is: " << vlivetime[Nevts - evtcounter] << endl;
		if (liveentry == -1 || dum1 == -1 || iev == -1)
		{
			std::cout << "Reading logfile failed skipping" << std::endl;
			vlivetime.clear();
			return false;
		}
		// vlivetime[Nevts - evtcounter] = (double)liveentry / 6e8;
		// cout << iev << ", live entry is: " << liveentry << ", livetime is: " << vlivetime[Nevts - evtcounter] << endl;
	}
	loginfile.close();
	return true;
}

static void show_usage(string name)
{
	cout << " Usage : ./DDC10_data_readout [-co] file1 " << name << " Options:\n"
	     << " -o : Name of output file.\n"
	     << " -i : Name of input file.\n"
	     << " -wd : Working directory\n"
	     << " -wform : Waveform channel\n"
	     << " -t : trigger channel\n"
	     << " -invert: invert waveform\n"
	     << " -trigger : invert trigger pulse \n"
	     << " -pt : pulse threshold\n"
	     << " -tri : Traiangle smoothing enabled\n"
	     << " -debug : Get in the debugging mode.\n"
	     << " -sit : number of smoothing interation\n"
	     << " -h or --help : Show the usage\n"
	     << " Enjoy ! -Ryan Wang" << endl;
}
int main(int argc, char* argv[])
{
	string filename;
	string outfilename;
	string working_dir;
	string baseline_file;

	bool use_trigger = false;
	bool trigger_inversion = false;
	bool invert_waveform = false;
	bool readlogs = false;

	int num_sams = 0;
	int trig_channel;
	int wform_channel = 1;

	ifstream logfile;

	if (argc < 4)
	{
		show_usage(argv[0]);
		return 1;
	}
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
		else if (arg == "-i")
		{
			filename = argv[i + 1];
		}
		else if (arg == "-o")
		{
			outfilename = argv[i + 1];
		}
		else if (arg == "-t")
		{
			trig_channel = atoi(argv[i + 1]);
			use_trigger = true;
		}
		else if (arg == "-wform")
		{
			wform_channel = atoi(argv[i + 1]);
		}
		else if (arg == "-bs")
		{
			baseline_samples_set = atoi(argv[i + 1]);
		}
		else if (arg == "-bf")
		{
			baseline_file = argv[i + 1];
			use_basefile = true;
		}
		else if (arg == "-pt")
		{
			pulseThresh = atof(argv[i + 1]);
		}
		else if (arg == "-win")
		{
			windowSize = atof(argv[i + 1]);
		}
		else if (arg == "-trigger")
		{
			trigger_inversion = true;
		}
		else if (arg == "-invert")
		{
			invert_waveform = true;
		}
		else if (arg == "-sams")
		{
			lim_sams = true;
			num_sams = atof(argv[i + 1]);
		}
		else if (arg == "-debug")
		{
			debug_mode = true;
		}
		else if (arg == "-readlogs")
		{
			readlogs = true;
			// pydir = argv[i + 1];
		}
		else if (arg == "-spe")
		{
			windowstart = atoi(argv[i + 1]);
			windowfin = atoi(argv[i + 2]);
		}
	}

	double fixedbase;
	double fixedrms;
	if (use_basefile)
	{
		TFile* basefile = new TFile(baseline_file.c_str(), "READ");
		if (basefile == NULL)
		{
			cout << "Couldn't open baseline file | Using standard methods instead" << endl;
			use_basefile = false;
		}
		else
		{
			TVectorD* fbase = (TVectorD*)basefile->Get("baseline");
			TVectorD* frms = (TVectorD*)basefile->Get("rms");
			fixedbase = fbase[0][0];
			fixedrms = frms[0][0];
			delete fbase;
			delete frms;
		}
		basefile->Close();
	}

	short int datum;
	int dummy;
	int mask, size;
	char open_filename[200]; // full input filename
	sprintf(open_filename, "%s/%s", working_dir.c_str(), filename.c_str());
	cout << " We are opening " << open_filename << endl;
	fin.open(open_filename, ios::binary | ios::in | ios::ate);

	if (fin.is_open())
	{
		// memblock.resize(size);
		size = fin.tellg();
		fin.seekg(0, ios::beg);
		// cout<<"size: "<<size<<endl;
		// numevts = new char [5];
		fin.read((char*)&Nevts, sizeof(Nevts));
		cout << Nevts << " Events" << endl;

		fin.read((char*)&number_of_samples, sizeof(number_of_samples));
		cout << number_of_samples << " Samples" << endl;

		fin.read((char*)&mask, sizeof(mask));
		cout << "Mask: " << mask << endl;
		Nchannels = calcnumchannels(mask);
		cout << Nchannels << " channels" << endl;

		fin.read((char*)&dummy, sizeof(dummy));
	}
	else
	{
		cout << "Failed to open file" << endl;
		return -1;
	}
	int predsize = Nevts * Nchannels * (2 * 4 + 2 * number_of_samples + 4) + 4 * 4;
	if (size < predsize)
	{
		cout << "Warning::Size predicted from header is greater than actual size" << endl;
		return -1;
	}

	if (wform_channel >= Nchannels || wform_channel < 0 || ((trig_channel >= Nchannels || trig_channel < 0) && use_trigger))
	{
		cout << "Channel numbers given do not match file header, double check your "
		        "inputs."
		     << endl;
		return -1;
	}
	int evtsize = Nchannels * (2 * 4 + 2 * number_of_samples + 4);
	int buffsize;
	if ((Nevts * evtsize) > 335544320)
		buffsize = 335544320 / evtsize;
	else
		buffsize = Nevts;
	cout << "Using buffer of " << buffsize << " events" << endl;

	if (readlogs)
		if (!readlivetime(open_filename))
		{
			readlogs = false;
		}

	// Plots for debugging pulse finding algorithm
	TGraph* t11;
	TGraph* t22;
	TGraph* t33;
	TGraph* t44;
	TGraph* t55;

	// char linea[200];// temp char for line in the file
	// TString* buff=new TString();
	char out_filename[200]; // full output filename
	int iana = 0; // TGraph counter
	int pcount = 0; // Pulse counter
	double temp_sum = 0;
	double rms_value;

	sprintf(out_filename, "%s/%s", working_dir.c_str(), outfilename.c_str());
	cout << " Out put filename is : " << out_filename << endl;
	TFile* fout = new TFile(out_filename, "RECREATE");

	TH1D* h_sum = new TH1D(("ADC_sum_waveform" + filename).c_str(), ("#font[132]{WFD " + filename + " SumWaveForm}").c_str(), 10000, 0, 10000);
	h_sum->SetXTitle("#font[132]{Sample (2ns)}");
	h_sum->GetXaxis()->SetLabelFont(132);
	h_sum->GetYaxis()->SetLabelFont(132);
	// Tetsing the dark hit counter
	TH1F* dark_hits = new TH1F("dark_rate", "dark_rate", 50000, 0, 50000);
	// dark_hits->SetBit(TH1::kCanRebin);

	// Create Tree to store properties of pulses found by pulse finder
	float waveforms[8192];
	float trigger_t;
	double livetime;
	TTree* event = new TTree("event", "Event tree");
	// event->Branch("npulses",&npulses,"npulses/I");
	event->Branch("nSamples", &number_of_samples, "number_of_samples/I");
	event->Branch("raw_waveforms", waveforms, "waveforms[number_of_samples]/F");

	if (readlogs)
		event->Branch("dLiveTime_s", &livetime, "livetime/D");

	event->Branch("fCharge_pC", &event_charge, "event_charge/F");
	event->Branch("fChargePrompt_pC", &event_charge_ten, "event_charge_ten/F");
	event->Branch("fBaseline_V", &event_baseline, "event_baseline/F");
	event->Branch("fBaselinerms_V", &event_rms, "event_rms/F");
	event->Branch("bIsGood", &isgood, "isgood/O");
	event->Branch("nPulses", &npulses, "npulses/I");
	event->Branch("fPulseHeight_V", amplitude, "amplitude[npulses]/F");
	event->Branch("fPulseRightEdge", pr, "pr[npulses]/F");
	event->Branch("fPulseLeftEdge", pl, "pl[npulses]/F");
	event->Branch("fPulseCharge_pC", charge_v, "charge_v[npulses]/F");
	event->Branch("fPulsePeakTime", amplitude_position, "amplitude_position[npulses]/F");
	event->Branch("fCalibratedTime", CalibratedTime, "CalibratedTime[npulses]/F");
	// event->Branch("windowratio", &windowRatio);
	event->Branch("fBigStep", biggeststep, "biggeststep[npulses]/F");
	event->Branch("fPulseLength5", pulse_length5, "pulse_length5[npulses]/F");
	event->Branch("fPulseLength25", pulse_length25, "pulse_length25[npulses]/F");
	event->Branch("fPulseLength50", pulse_length50, "pulse_length50[npulses]/F");
	event->Branch("fPulseLength75", pulse_length75, "pulse_length75[npulses]/F");
	// event->Branch("fPulseLength80", pulse_length80,"pulse_length[npulses]/F");
	event->Branch("fPulseLength90", pulse_length90, "pulse_length90[npulses]/F");
	event->Branch("fPulseLength95", pulse_length95, "pulse_length95[npulses]/F");
	event->Branch("fPulseLength99", pulse_length99, "pulse_length99[npulses]/F");

	if (use_trigger)
	{
		event->Branch("fTriggerTime", &trigger_t, "trigger_t/F");
		event->Branch("fTriggerHeight_V", &triggerHeight, "triggerHeight/F");
		event->Branch("fTriggerWidth", &triggerWidth, "triggerWidth/F");
		event->Branch("fWindowCharge_pC", &event_windowCharge, "event_windowCharge/F");
	}
	// Store the waveform plot for debugging
	TCanvas* waveplot;
	vector<float> baseline_sweep;
	vector<float> trigwaveform;

	int skip = 0;
	int readin = 0;
	int lastadd = 0;
	for (int sweep = 0; sweep < Nevts; sweep++)
	{
		while ((readin) < (sweep + 1))
		{
			int arrsize = buffsize * evtsize / 2;
			if ((Nevts - readin) < buffsize)
			{
				arrsize = (Nevts - readin) * evtsize / 2;
			}
			buff = new short int[arrsize];
			fin.read((char*)&buff[0], arrsize * sizeof(buff[0]));
			lastadd = readin;
			readin += arrsize * 2 / evtsize;

			cout << "Read in " << readin << " events so far" << endl;
		}

		current_sweep = sweep;

		if (use_trigger)
		{
			signal_start = false;
			getwaveform(trigwaveform, trig_channel, lastadd, (trigger_inversion ? -1.0 : 1.0), true);
			trigger_t = Trigger_info(trigwaveform);
			trigwaveform.clear();
		}

		if (sweep % 1000 == 0)
		{
			cout << " This is sweep : " << sweep << endl;
			cout << "Trigger time: " << trigger_t << endl;
		}

		signal_start = true;
		getwaveform(raw_waveform, wform_channel, lastadd, (invert_waveform ? -1.0 : 1.0));

		std::copy(raw_waveform.begin(), raw_waveform.end(), waveforms);
		// wforms_tree->Fill();
		npulses = 0.0;
		double thisbase;
		rms_value = (use_basefile ? fixedrms : 0);

		thisbase = (use_basefile ? fixedbase : baseline_rms(baselinev, raw_waveform, &rms_value));
		extract_event(raw_waveform, thisbase, rms_value, (lim_sams ? num_sams : number_of_samples), (use_trigger ? trigger_t : 0));

		if (debug_mode)
		{
			cout << " basline is  : " << thisbase << " rms is : " << rms_value << endl;
			getchar();
		}
		if (readlogs)
		{
			livetime = vlivetime[Nevts - sweep - 1];
			vlivetime.pop_back();
		}
		event_baseline = thisbase;
		event_rms = rms_value;
		event->Fill();
		// event_time.push_back(number_of_samples);
		baseline_sweep.push_back(thisbase); // save baseline for checking baseline shifting
		double drate = npulses;

		dark_hits->Fill(drate);
		// npeaks.push_back(npulses);
		baselinev.clear();
		// fill canvases
		for (int sam = 0; sam < number_of_samples; sam++)
		{
			h_sum->Fill(sam, raw_waveform[sam] - thisbase);
		}
		// fill canvases
		if (sweep < 100)
		{
			t11 = new TGraph();
			t22 = new TGraph();
			t33 = new TGraph();
			double rmax = 0, rmin = 1e16;
			for (int j = 0; j < (int)pulse_left_edge.size(); j++)
			{
				t22->SetPoint(j, pulse_left_edge[j], startv[j]);
				t33->SetPoint(j, pulse_right_edge[j], endv[j]);
			}
			for (int sam = 0; sam < number_of_samples; sam++)
			{
				t11->SetPoint(sam, sam, raw_waveform[sam]);
				if (std::fabs(raw_waveform[sam]) > rmax)
					rmax = std::fabs(raw_waveform[sam]);
				if (raw_waveform[sam] < rmin)
					rmin = raw_waveform[sam];
			}
			char plotname[30];
			sprintf(plotname, "waveform%d", sweep);
			waveplot = new TCanvas(plotname);
			TLine* line = new TLine(0, thisbase, number_of_samples, thisbase);
			TLine* line2 = new TLine(0, -rms_value + thisbase, number_of_samples, -rms_value + thisbase);
			TLine* line3 = new TLine(0, rms_value + thisbase, number_of_samples, rms_value + thisbase);

			line->SetLineColor(3);
			line->SetLineStyle(3);
			line->SetLineWidth(3);
			line2->SetLineColor(5);
			line2->SetLineStyle(3);
			line2->SetLineWidth(3);
			line3->SetLineColor(5);
			line3->SetLineStyle(3);
			line3->SetLineWidth(3);
			t22->SetMarkerColor(2);
			t22->SetMarkerStyle(3);
			t22->SetMarkerSize(3);
			t33->SetMarkerColor(4);
			t33->SetMarkerStyle(3);
			t33->SetMarkerSize(3);
			t11->Draw("alp");
			t22->Draw("p");
			t33->Draw("p");
			line->Draw("");
			line2->Draw("");
			line3->Draw("");
			if (use_trigger)
			{
				TLine* line4 = new TLine(trigger_t, rmax, trigger_t, rmin * 1.1);
				line4->SetLineColor(46);
				line4->SetLineStyle(3);
				line4->SetLineWidth(3);
				line4->Draw("");
			}
			waveplot->Write();
			cout << " plotting the waveform, this is sweep : " << sweep << endl;
			delete line;
			delete line2;
			delete line3;
			waveplot->Close();
		}
		pulse_left_edge.clear();
		pulse_right_edge.clear();
		startv.clear();
		endv.clear();
		raw_waveform.clear();
	}

	TGraph* baseline_plot = new TGraph();
	for (int i = 0; i < baseline_sweep.size(); i++)
	{
		baseline_plot->SetPoint(i, i, baseline_sweep[i]);
	}

	// Baseline plot
	TCanvas* bplot = new TCanvas("bplot", "bplot");
	baseline_plot->SetMarkerStyle(22);
	baseline_plot->Draw("AP");

	cout << " Total sweeps is : " << Nevts << endl;
	h_sum->Scale(1.0 / (double)Nevts);
	event->Write();
	bplot->Write();
	h_sum->Write();
	fout->Write();
	fout->Close();

	return 0;
}
