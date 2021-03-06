/////////////////////////////////////////////////////////////////////////////////////////
/*
To compile :
g++ DDC10_pkl_data_readout.cc -o DDC10_pkl_data_readout -I${ROOTSYS}/include/ `root-config --cflags --libs` `python3-config --cflags --ldflags --libs`
-fPIC
 */
// To execute (help infomation gives detail utility) : ./DDC10_data_readout -h
/* Revision log :
 *
        4/25/2018 RW : Code for reading scope's text file and do the pulse
 finding. 7/23/2018 RW : Add option for user to identify each item easily.



*/
/////////////////////////////////////////////////////////////////////////////////////////
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <cstdio>
#include <cstdlib>
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
//#include <libRATEvent.so>
#include <Python.h>
#include <numpy/arrayobject.h>
#include <numpy/ndarrayobject.h>

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
bool goodevent = false;
bool lim_sams = false;

bool smoothing = false;
bool rolling = false;
bool triangle = false;
bool boxsmoothing = false;
bool use_basefile = false;

// define pulse finding parameters
double pulseThresh = 8.0;
double tPulseThresh = 10.0;
double windowSize = 3.0;
double edgeThresh = 3.0;
double lookforward = 10.0;
int current_sweep = 0;
double number_of_peaks = 0.0;
int adc_per_Volt = 8192;
const double ns = 1e-9;
double start_time = 0;
double timescale = 10e-9 / ns;
double resistance = 50;
int baseline_samples_set = 160;
int MovingWindowSize;
int iteration = 0.0;
int pth = 5.0;
int promptwindow = 300;
// input parameters
int number_of_samples;
int Nchannels;
int Nevts;
short int* buff;

// internal run tracking
vector<float> raw_waveform;

vector<float> trigwaveform;

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
const int kMaxPulses = 800;

float amplitude[kMaxPulses];
float charge_v[kMaxPulses];
float amplitude_position[kMaxPulses];
float pl[kMaxPulses];
float pr[kMaxPulses];
float biggeststep[kMaxPulses];
// reconstruction params
float pulse_peak20[kMaxPulses];
float pulse_length99[kMaxPulses];
float pulse_length95[kMaxPulses];
float pulse_length90[kMaxPulses];
float pulse_length80[kMaxPulses];
float pulse_length75[kMaxPulses];
float pulse_length50[kMaxPulses];
float pulse_length25[kMaxPulses];
// float pulse_length10[kMaxPulses];
float pulse_length5[kMaxPulses];
float pulse_length1[kMaxPulses];
float pulse_length05[kMaxPulses];
float CalibratedTime[kMaxPulses];
// float windowratio;
// float pulsebaseline_rms;

float triggerHeight;
float triggerStart;
float triggerStartSam;
float triggerRising05;
float triggerRising1;
float triggerRising5;
float triggerPosition;
float triggerWidth;
int nTrigs = 0;

float event_charge;
float event_charge_ten;
float event_baseline;
float event_rms;
int npulses = 0;
// vector<double> vlivetime;

TH1D* h_sum;
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
// Pulse Finding
void extract_event(vector<float>& v, double b, double rms, int nos, int trigger = 0, bool trig = false)
{

	double pThresh = (trig ? tPulseThresh : pulseThresh) * rms * windowSize;
	double eThresh = edgeThresh * rms * windowSize;

	double temp_charge = 0;
	double temp_ten_charge = 0;
	double pulse_height_thresh = pth * rms;
	// cout<<" vector size is : "<<v.size()<<endl;
	// getchar();
	// Let's looking for the Pulses
	npulses = 0;
	goodevent = true;
	for (int i = 0; i < v.size() - windowSize; i++)
	{
		// std::cout << "Sample " << i << std::endl;
		double integral = SimpsIntegral(v, b, i, i + windowSize);
		int left = 0;
		int right = 0;
		int temp_peak = 0;
		double temp_startv = 0;
		double temp_endv = 0;

		double temp_bigstep = 0;

		if ((integral) > pThresh)
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
			double tempq05 = 0, tempq1 = 0, tempq5 = 0, tempq25 = 0, tempq50 = 0, tempq75 = 0, tempq90 = 0, tempq95 = 0, tempq99 = 0;
			for (int j = left; j < right; j++)
			{
				double s = v[j] - b;
				if (s > max)
				{
					max = s;
					temp_peak = j;
					if (j > 0 && (v[j] - v[j - 1]) > temp_bigstep)
						temp_bigstep = v[j] - v[j - 1];
				}
				// if((right-j)<(temp_peak-left))
				double ratio = SimpsIntegral(v, b, left, j) / totalq;
				if (ratio <= 0.005 && tempq05 < (j - left))
					tempq05 = j - left;
				if (ratio <= 0.01 && tempq1 < (j - left))
					tempq1 = j - left;
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

			if ((SimpsIntegral(v, b, left, right)) <= eThresh)
			{
				i = right - 1;
				continue;
			}
			i = right;

			if ((max) < pulse_height_thresh)
				continue;

			npulses++;
			if (signal_start)
			{
				startv.push_back(temp_startv);
				endv.push_back(temp_endv);

				amplitude[npulses - 1] = max;
				amplitude_position[npulses - 1] = temp_peak * timescale;
				pl[npulses - 1] = (left - triggerStartSam) * timescale;
				pr[npulses - 1] = (right - triggerStartSam) * timescale;
				charge_v[npulses - 1] = 1e3 * timescale * SimpsIntegral(v, b, left, right) / resistance;

				CalibratedTime[npulses - 1] = timescale * (temp_peak - triggerStartSam);
				// windowratio[npulses-1] =dratio;
				// pulsebaseline_rms[npulses-1] =rms;
				// event_n[npulses-1]=current_sweep;
				biggeststep[npulses - 1] = temp_bigstep;
				pulse_length05[npulses - 1] = tempq05 * timescale;
				pulse_length1[npulses - 1] = tempq1 * timescale;
				pulse_length5[npulses - 1] = tempq5 * timescale;
				pulse_length25[npulses - 1] = tempq25 * timescale;
				pulse_length50[npulses - 1] = tempq50 * timescale;
				pulse_length75[npulses - 1] = tempq75 * timescale;
				// pulse_length80[npulses-1]=tempq80;
				pulse_length90[npulses - 1] = tempq90 * timescale;
				pulse_length95[npulses - 1] = tempq95 * timescale;
				pulse_length99[npulses - 1] = tempq99 * timescale;

				temp_charge += charge_v[npulses - 1];
				if (i - trigger < promptwindow)
					temp_ten_charge += charge_v[npulses - 1];
				if (left < baseline_samples_set)
					goodevent = false;
				double temp_20start = left, temp_20fin = right;
				for (int j = temp_peak; j > left; j--)
				{
					double s = v[j] - b;
					if (s >= 0.2 * max)
					{
						temp_20start = j;
					}
					else
						break;
				}
				for (int j = temp_peak; j < right; j++)
				{
					double s = v[j] - b;
					if (s >= 0.2 * max)
					{
						temp_20fin = j;
					}
					else
						break;
				}
				pulse_peak20[npulses - 1] = (temp_20fin - temp_20start) * timescale;
			}
			else
			{
				triggerHeight = max;
				triggerStart = left * timescale;
				triggerStartSam = left;
				triggerRising05 = tempq05 * timescale;
				triggerRising1 = tempq1 * timescale;
				triggerRising5 = tempq5 * timescale;
				triggerPosition = (temp_peak - left) * timescale;
				triggerWidth = timescale * (right - left);
			}
			if (!signal_start && npulses > 1)
				continue;
			pulse_left_edge.push_back(left);
			pulse_right_edge.push_back(right);

		} // if statement
	}
	// getchar();
	if (!trig)
	{
		event_charge_ten = temp_ten_charge;
		event_charge = temp_charge;
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
		time = (pulse_left_edge[0]);
	}
	nTrigs = npulses;
	pulse_left_edge.clear();
	pulse_right_edge.clear();
	trigbaselinev.clear();
	// waveform.clear();
	return time;
}

void getwaveform(vector<float>& v, int evt, PyArrayObject* arr, float mult = 1, bool trig = false)
{
	/*if(current_sweep%100000==0){
	        cout<<"starting at "<<starti<<" in buffer"<<endl;
	        cout<<"Evt "<<(current_sweep-numread)<<" of this buffer"<<endl;
	}*/
	double datum;
	for (int i = 0; i < number_of_samples; i++)
	{
		// cout << "reading sample " << i << endl;
		// std::cout << "hi" << std::endl;
		datum = *((npy_double*)(PyArray_GETPTR2(arr, evt, i)));
		// std::cout << "bi" << std::endl;
		datum *= mult;
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
	std::string pydir;

	bool use_trigger = true;
	bool trigger_inversion = true;
	bool invert_waveform = true;

	int num_sams = 0;
	int trig_channel;
	int wform_channel = 1;

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
		else if (arg == "-pydir")
		{
			pydir = argv[i + 1];
			// pydir.append("/readlogs.py");
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

	char open_filename[200]; // full input filename
	sprintf(open_filename, "%s/%s", working_dir.c_str(), filename.c_str());
	cout << " We are opening " << open_filename << endl;

	std::cout << "Initializing python" << std::endl;
	Py_Initialize();
	import_array();
	std::cout << "Updating syspath with pydir: " << pydir << std::endl;
	PyObject* sysPath = PySys_GetObject("path");
	PyList_Append(sysPath, PyUnicode_FromFormat("%s", pydir.c_str()));

	std::cout << "Loading module" << std::endl;
	// Load module
	PyObject* pName = PyUnicode_FromString("ReadNpyfiles");
	PyObject* pModule = PyImport_Import(pName);
	PyArrayObject* pData = NULL;
	PyArrayObject* pXData = NULL;
	PyArrayObject* pTrig = NULL;
	Py_DECREF(pName);
	Py_DECREF(sysPath);
	int nd1 = 0, nd2 = 0;
	if (pModule != NULL)
	{
		std::cout << "Py Module Found" << std::endl;

		// Get function from module
		PyObject* pFunc = PyObject_GetAttrString(pModule, "loadpklfile");
		if (pFunc && PyCallable_Check(pFunc))
		{
			PyObject* pargs = PyTuple_New(2);
			PyObject* pval = NULL;
			pval = PyUnicode_FromFormat("%s", open_filename);
			PyTuple_SetItem(pargs, 0, pval);
			pval = PyLong_FromLong(wform_channel);
			PyTuple_SetItem(pargs, 1, pval);
			if (pargs != NULL)
				pval = PyObject_CallObject(pFunc, pargs);
			else
				pval = 0;
			if (pval != NULL && PyTuple_Check(pval))
			{
				// PyObject *pTemp = PyTuple_GetItem(pval, 0);
				PyObject* p1 = 0;
				p1 = PyTuple_GetItem(pval, 0);
				// pXData = reinterpret_cast<PyArrayObject*>(p1);
				pXData = reinterpret_cast<PyArrayObject*>(PyArray_FromObject(p1, NPY_FLOAT64, 1, 1));
				p1 = PyTuple_GetItem(pval, 1);
				pTrig = reinterpret_cast<PyArrayObject*>(PyArray_FromObject(p1, NPY_FLOAT64, 1, 2));
				p1 = PyTuple_GetItem(pval, 2);
				pData = reinterpret_cast<PyArrayObject*>(PyArray_FromObject(p1, NPY_FLOAT64, 1, 2));
				// pData2 = reinterpret_cast<PyArrayObject*>(PyTuple_GetItem(pval, 3));
				if (pTrig != NULL && pXData != NULL && pData != NULL)
				{
					npy_intp* pshape = PyArray_SHAPE(pXData);
					int indx = pshape[0] - 1;
					timescale = *((double*)PyArray_GETPTR1(pXData, indx));
					start_time = *((double*)PyArray_GETPTR1(pXData, 0));
					// std::cout << "Start Time is " << start_time / ns << " ns" << std::endl;
					// std::cout << "End Time at index " << indx << " is " << timescale / ns << " ns" << std::endl;
					timescale -= start_time;
					timescale *= 1.0 / (double)pshape[0];
					timescale *= 1.0 / ns;
					start_time *= 1.0 / ns;
					// std::cout << "Start Time is " << start_time << " ns" << std::endl;
					std::cout << "Time Scale is " << timescale << " ns per sample" << std::endl;
					std::cout << "X Array has " << PyArray_NDIM(pXData) << " dims and " << pshape[0] << " Elements" << std::endl;
					int ndim;

					pshape = PyArray_SHAPE(pTrig);
					ndim = PyArray_NDIM(pTrig);
					std::cout << "Trig Array has " << ndim << " dims and " << pshape[0];
					for (int idims = 1; idims < ndim; idims++)
						std::cout << " x " << pshape[idims];
					std::cout << " Elements" << std::endl;
					number_of_samples = pshape[1];
					Nevts = pshape[0];
					pshape = PyArray_SHAPE(pData);
					ndim = PyArray_NDIM(pData);

					std::cout << "Data Array has " << ndim << " dims and " << pshape[0];
					for (int idims = 1; idims < ndim; idims++)
						std::cout << " x " << pshape[idims];
					std::cout << " Elements" << std::endl;

					Py_INCREF(pData);
					Py_INCREF(pTrig);
					// Py_XDECREF(pXData);
				}
				// cout << "Array is int? " << PyArray_ISINTEGER(pData) << endl;
				// Py_XDECREF(pTemp);
			}
			else
			{
				PyErr_Print();
				std::cout << "Failed to get result" << std::endl;
			}

			Py_XDECREF(pval);
			Py_DECREF(pargs);
		}
		else
		{
			if (PyErr_Occurred())
				PyErr_Print();
			std::cout << "Couldn't find loadnpyfile " << std::endl;
		}
		Py_XDECREF(pFunc);
		Py_XDECREF(pModule);
	}
	else
	{
		PyErr_Print();
		std::cout << "Failed to load module" << std::endl;
	}
	// cout << "Typechecking pData" << endl;
	if (pData == NULL || pTrig == NULL)
	{
		std::cout << "Failed to read data file. Exiting" << std::endl;
		Py_FinalizeEx();
		return -1;
	}
	//	Py_XDECREF(pXData);

	// Plots for debugging pulse finding algorithm
	TGraph* t11;
	TGraph* t22;
	TGraph* t33;
	TGraph* t44;
	TGraph* t55;

	// char linea[200];// temp char for line in the file
	// TString* buff=new TString();
	char out_filename[200]; // full output filename
	double rms_value;

	sprintf(out_filename, "%s/%s", working_dir.c_str(), outfilename.c_str());
	cout << " Out put filename is : " << out_filename << endl;
	TFile* fout = new TFile(out_filename, "RECREATE");

	h_sum = new TH1D(("ADC_sum_waveform" + filename).c_str(), ("#font[132]{WFD " + filename + " SumWaveForm}").c_str(), 10000, 0, 10000);
	h_sum->SetXTitle("#font[132]{Sample (2ns)}");
	h_sum->GetXaxis()->SetLabelFont(132);
	h_sum->GetYaxis()->SetLabelFont(132);
	// Tetsing the dark hit counter
	TH1F* dark_hits = new TH1F("dark_rate", "dark_rate", 25000, 0, 50000);
	// dark_hits->SetBit(TH1::kCanRebin);

	// Create Tree to store properties of pulses found by pulse finder
	float waveforms[10001];
	float trigger_t;
	double livetime;
	TTree* event = new TTree("event", "Event tree");
	// event->Branch("npulses",&npulses,"npulses/I");
	event->Branch("nSamples", &number_of_samples, "number_of_samples/I");
	event->Branch("raw_waveforms", waveforms, "waveforms[number_of_samples]/F");

	event->Branch("bIsGood", &goodevent, "goodevent/O");
	event->Branch("fCharge_pC", &event_charge, "event_charge/F");
	event->Branch("fChargePrompt_pC", &event_charge_ten, "event_charge_ten/F");
	event->Branch("fBaseline_V", &event_baseline, "event_baseline/F");

	event->Branch("nPulses", &npulses, "npulses/I");
	event->Branch("fPulseHeight_V", amplitude, "amplitude[npulses]/F");
	event->Branch("fPulseRightEdge", pr, "pr[npulses]/F");
	event->Branch("fPulseLeftEdge", pl, "pl[npulses]/F");
	event->Branch("fPulseCharge_pC", charge_v, "charge_v[npulses]/F");
	event->Branch("fPulsePeakTime", amplitude_position, "amplitude_position[npulses]/F");
	event->Branch("fCalibratedTime", CalibratedTime, "CalibratedTime[npulses]/F");
	event->Branch("fPulsePeak20", pulse_peak20, "pulse_peak20[npulses]/F");
	// event->Branch("windowratio", &windowRatio);
	event->Branch("fBigStep", biggeststep, "biggeststep[npulses]/F");
	event->Branch("fPulseLength05", pulse_length05, "pulse_length05[npulses]/F");
	event->Branch("fPulseLength1", pulse_length1, "pulse_length1[npulses]/F");
	event->Branch("fPulseLength5", pulse_length5, "pulse_length5[npulses]/F");
	// event->Branch("fPulseLength10", pulse_length10, "pulse_length10[npulses]/F");
	event->Branch("fPulseLength25", pulse_length25, "pulse_length25[npulses]/F");
	event->Branch("fPulseLength50", pulse_length50, "pulse_length50[npulses]/F");
	event->Branch("fPulseLength75", pulse_length75, "pulse_length75[npulses]/F");
	// event->Branch("fPulseLength80", pulse_length80,"pulse_length[npulses]/F");
	event->Branch("fPulseLength90", pulse_length90, "pulse_length90[npulses]/F");
	event->Branch("fPulseLength95", pulse_length95, "pulse_length95[npulses]/F");
	event->Branch("fPulseLength99", pulse_length99, "pulse_length99[npulses]/F");

	if (use_trigger)
	{
		event->Branch("nTriggers", &nTrigs, "nTriggers/I");
		event->Branch("fTriggerStart", &triggerStart, "triggerStart/F");
		event->Branch("fTriggerStartSam", &triggerStartSam, "triggerStartSam/F");
		event->Branch("fTriggerRising05", &triggerRising05, "triggerRising05/F");
		event->Branch("fTriggerRising1", &triggerRising1, "triggerRising1/F");
		event->Branch("fTriggerRising5", &triggerRising5, "triggerRising5/F");
		event->Branch("fTriggerTime", &triggerPosition, "triggerPosition/F");
		event->Branch("fTriggerHeight_V", &triggerHeight, "triggerHeight/F");
		event->Branch("fTriggerWidth", &triggerWidth, "triggerWidth/F");
	}
	// Store the waveform plot for debugging
	TCanvas* waveplot;
	vector<float> baseline_sweep;

	for (int sweep = 0; sweep < Nevts; sweep++)
	{
		/*
		if ((readin) < (sweep + 1)) {
		  int arrsize = buffsize * evtsize / 2;
		  if ((Nevts - readin) < buffsize) {
		    arrsize = (Nevts - readin) * evtsize / 2;
		  }
		  buff = new short int[arrsize];
		  fin.read((char *)&buff[0], arrsize * sizeof(buff[0]));
		  lastadd = readin;
		  readin += arrsize * 2 / evtsize;

		  cout << "Read in " << readin << " events so far" << endl;
		}
		*/
		if (sweep % 100 == 0)
		{
			cout << " This is sweep : " << sweep << endl;
		}
		current_sweep = sweep;

		if (use_trigger)
		{
			signal_start = false;
			getwaveform(trigwaveform, sweep, pTrig, (trigger_inversion ? -1.0 : 1.0), true);
			trigger_t = Trigger_info(trigwaveform);
			trigwaveform.clear();
		}

		signal_start = true;

		// cout << "reading Waveform" << endl;
		getwaveform(raw_waveform, sweep, pData, (invert_waveform ? -1.0 : 1.0));
		// cout << "Waveform read" << endl;

		std::copy(raw_waveform.begin(), raw_waveform.end(), waveforms);
		// wforms_tree->Fill();
		// nPulses = 0.0;
		double thisbase;
		rms_value = (use_basefile ? fixedrms : 0);

		thisbase = (use_basefile ? fixedbase : baseline_rms(baselinev, raw_waveform, &rms_value));
		extract_event(raw_waveform, thisbase, rms_value, number_of_samples, (use_trigger ? trigger_t : 0));

		if (debug_mode)
		{
			cout << " basline is  : " << thisbase << " rms is : " << rms_value << endl;
			getchar();
		}

		event_baseline = thisbase;
		event_rms = rms_value;
		event->Fill();
		// event_time.push_back(number_of_samples);
		baseline_sweep.push_back(thisbase); // save baseline for checking baseline shifting
		double drate = npulses;
		drate *= 1.0 / (double)(timescale * ns * number_of_samples);
		dark_hits->Fill(drate);
		baselinev.clear();

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
				t22->SetPoint(j, pulse_left_edge[j] * timescale, startv[j]);
				t33->SetPoint(j, pulse_right_edge[j] * timescale, endv[j]);
			}
			for (int sam = 0; sam < number_of_samples; sam++)
			{
				t11->SetPoint(sam, sam * timescale, raw_waveform[sam]);
				if (std::fabs(raw_waveform[sam]) > rmax)
					rmax = std::fabs(raw_waveform[sam]);
				if (raw_waveform[sam] < rmin)
					rmin = raw_waveform[sam];
			}
			char plotname[30];
			sprintf(plotname, "waveform%d", sweep);
			waveplot = new TCanvas(plotname);
			TLine* line = new TLine(0, thisbase, number_of_samples * timescale, thisbase);
			TLine* line2 = new TLine(0, -rms_value + thisbase, number_of_samples * timescale, -rms_value + thisbase);
			TLine* line3 = new TLine(0, rms_value + thisbase, number_of_samples * timescale, rms_value + thisbase);

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
				TLine* line4 = new TLine(triggerStartSam * timescale, rmax, triggerStartSam * timescale, rmin * 1.1);
				line4->SetLineColor(46);
				line4->SetLineStyle(3);
				line4->SetLineWidth(3);
				line4->Draw("");
			}
			waveplot->Write();
			cout << " plotting the waveform, this is sweep : " << sweep << endl;
		}
		pulse_left_edge.clear();
		pulse_right_edge.clear();
		startv.clear();
		endv.clear();
		raw_waveform.clear();
		// skip += ;
	}
	// PyArray_XDECREF(pData);
	// PyArray_XDECREF(pXData);
	// PyArray_XDECREF(pTrig);
	// Py_XDECREF(pData);
	// Py_XDECREF(pXData);
	// Py_XDECREF(pTrig);
	Py_FinalizeEx();

	// cout << " after tree fill ! " << endl;
	TGraph* baseline_plot = new TGraph();
	for (int i = 0; i < (int)baseline_sweep.size(); i++)
	{
		baseline_plot->SetPoint(i, i, baseline_sweep[i]);
	}
	// Baseline plot
	TCanvas* bplot = new TCanvas("bplot", "bplot");
	baseline_plot->SetMarkerStyle(22);
	baseline_plot->Draw("AP");

	h_sum->Scale(1.0 / (double)Nevts);
	// TObject* nosinfo = new TObject();
	// nosinfo->SetUniqueID(number_of_samples);

	cout << " Total sweeps is : " << Nevts << endl;

	event->Write();
	bplot->Write();
	h_sum->Write();
	fout->Write();
	fout->Close();

	return 0;
}
