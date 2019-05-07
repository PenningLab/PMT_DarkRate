/////////////////////////////////////////////////////////////////////////////////////////
// To compile : g++ -I${ROOTSYS}/include/ `root-config --cflags --libs` -o
// DDC10_bin_data_readout DDC10_bin_data_readout.cc To execute (help infomation
// gives detail utility) : ./DDC10_data_readout -h
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
#include "numpy/arrayobject.h"
#include <Python.h>
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

bool lim_sams = false;

bool smoothing = false;
bool rolling = false;
bool triangle = false;
bool boxsmoothing = false;
bool use_basefile = false;

// define pulse finding parameters
double pulseThresh = 8.0;
double tpulseThresh = 8.0;
double windowSize = 3.0;
double edgeThresh = 3.0;
double lookforward = 3.0;
int current_sweep = 0;
double number_of_peaks = 0.0;
int adc_per_Volt = 8192;
double resistance = 0.005;
int baseline_samples_set = 160;
int MovingWindowSize;
int iteration = 0.0;
int pth = 5.0;

// vector<int> run_info;
vector<float> raw_waveform;
vector<float> amplitude;
vector<float> charge_v;
vector<float> startv;
vector<float> endv;
vector<float> start_pointv;
vector<float> end_pointv;
vector<float> baseline;
vector<float> baselinev;
vector<float> trigbaselinev;
vector<float> pulse_left_edge;
vector<float> pulse_right_edge;
vector<float> amplitude_position;
vector<float> pl;
vector<float> pr;
vector<float> npeaks;

vector<float> biggeststep;
vector<float> smoothingv;
// lzap like params
vector<float> pulse_length99;
vector<float> pulse_length95;
vector<float> pulse_length90;
vector<float> pulse_length80;
vector<float> pulse_length75;
vector<float> pulse_length50;
vector<float> pulse_length25;
vector<float> pulse_length5;

// vector<float> trigger_time;
vector<float> CalibratedTime;
vector<float> windowratio;
vector<float> pulsebaseline_rms;
vector<short int> event_n;

vector<float> smoothedBaseLine;

vector<float> triggerHeight;
vector<float> triggerPosition;
vector<float> triggerWidth;

vector<float> event_charge;
vector<float> event_charge_ten;
vector<float> event_baseline;
vector<float> event_rms;
vector<float> event_time;

int number_of_samples;
int Nchannels;
int Nevts;
short int *buff;

// Simpson Integral
double SimpsIntegral(const vector<float> &samples, double baseline, int start,
                     int end) {
  int len;
  double qsum = 0.0;
  if ((end - start) % 2 == 0) {
    /* If there are an even number of samples, then there are an odd
    number of intervals; but Simpson's rule works only on an even
    number of intervals. Therefore we use Simpson's method on the
    all but the final sample, and integrate the last interval
    using the trapezoidal rule */
    len = end - start - 1;
    qsum += (samples[end - 1] + samples[end - 2] - 2 * baseline) / 2.0;
  } else
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
void extract_event(vector<float> &v, double b, double rms, int nos, int trigger,
                   bool trig = false) {

  double pThresh = (trig ? tpulseThresh : pulseThresh) * rms * windowSize;
  double eThresh = edgeThresh * rms * windowSize;

  double temp_charge = 0;
  double temp_ten_charge = 0;
  double pulse_height_thresh = pth * rms;
  // cout<<" vector size is : "<<v.size()<<endl;
  // getchar();
  // Let's looking for the Pulses
  for (int i = 0; i < (int)v.size(); i++) {
    double integral = SimpsIntegral(v, b, i, i + windowSize);
    int left = 0;
    int right = 0;
    int temp_peak = 0;
    double temp_startv = 0;
    double temp_endv = 0;

    double temp_bigstep = 0;

    if (integral > pThresh) {
      if (debug_mode) {
        cout << " This is sample : " << i << " integral value is : " << integral
             << " pThresh is : " << pThresh << endl;
      }
      left = i;
      integral = 1.0e9;
      while (integral > eThresh && left > windowSize) {
        left--;
        integral = SimpsIntegral(v, b, left, left + windowSize);
        temp_startv = v[left];
        if (debug_mode) {

          cout << "Left is  : " << left << " integral is : " << integral
               << " eThresh is : " << eThresh << endl;
          getchar();
        }
      }

      integral = 1.0e9;
      right = i + windowSize;
      double thischarge = 0;
      bool end = false;
      while (!end) {
        while (integral > eThresh && right < number_of_samples - 1) {
          right++;
          integral = SimpsIntegral(v, b, right - windowSize, right);
          temp_endv = v[right];
        }
        end = true;
        int r = right;
        while (r < fmin((int)v.size() - 1, right + lookforward)) {
          r++;
          integral = SimpsIntegral(v, b, r - windowSize, r);
          if (integral > pThresh) {
            right = r;
            end = false;
            break;
          }
        }
      }

      double max = -1.0e9;
      double totalq = SimpsIntegral(v, b, left, right);
      double tempq5 = 0, tempq25 = 0, tempq50 = 0, tempq75 = 0, tempq80 = 0,
             tempq90 = 0, tempq95 = 0, tempq99 = 0;
      for (int j = left; j < right; j++) {
        double s = v[j] - b;
        if (s > max) {
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
        if (ratio <= 0.8 && tempq80 < (j - left))
          tempq80 = j - left;
        if (ratio <= 0.9 && tempq90 < (j - left))
          tempq90 = j - left;
        if (ratio <= 0.95 && tempq95 < (j - left))
          tempq95 = j - left;
        if (ratio <= 0.99 && tempq99 < (j - left))
          tempq99 = j - left;
      }
      if (right > nos)
        continue;

      // noise veto
      float width = (right - left);
      float nwidth = 2 * width;
      int nright = right + nwidth;
      int nleft = left - nwidth;
      if (nright > nos) {
        nright = nos;
        nleft -= (nleft > nwidth ? nwidth : 0);
      }
      if (nleft < 0) {
        nleft = 0;
        nright += (nright < (nos - nwidth) ? nwidth : nos);
      }
      double dratio =
          SimpsIntegral(v, b, nleft, nright) /
          (((nwidth + 1) / width) * SimpsIntegral(v, b, left, right));

      // cout<<" Peak is : "<<temp_peak<<" max is : "<<max<<endl;

      // if (temp_peak>0 &&temp_peak<8000){
      // if (thischarge<1.0)
      if (trig)
        break;

      //}
      // cout<<" This is sample : "<<i<<" Charge integral is :
      // "<<SimpsIntegral(v,b,left,right)<<endl;

      if (SimpsIntegral(v, b, left, right) <= eThresh) {
        i = right - 1;
        continue;
      }
      i = right;

      if (max < pulse_height_thresh)
        continue;

      number_of_peaks++;
      if (signal_start) {
        amplitude.push_back(max);
        amplitude_position.push_back(temp_peak);
        pl.push_back(left);
        pr.push_back(right);
        charge_v.push_back(SimpsIntegral(v, b, left, right) / resistance);
        startv.push_back(temp_startv);
        endv.push_back(temp_endv);
        CalibratedTime.push_back(temp_peak - trigger);
        windowratio.push_back(dratio);
        pulsebaseline_rms.push_back(rms);
        event_n.push_back(current_sweep);
        biggeststep.push_back(temp_bigstep);
        pulse_length5.push_back(tempq5);
        pulse_length25.push_back(tempq25);
        pulse_length50.push_back(tempq50);
        pulse_length75.push_back(tempq75);
        pulse_length80.push_back(tempq80);
        pulse_length90.push_back(tempq90);
        pulse_length95.push_back(tempq95);
        pulse_length99.push_back(tempq99);

        temp_charge += SimpsIntegral(v, b, left, right) / resistance;
        if (i < 300)
          temp_ten_charge += SimpsIntegral(v, b, left, right) / resistance;
      } else {
        triggerHeight.push_back(max);
        triggerPosition.push_back(temp_peak);
        triggerWidth.push_back(right - left);
      }

      pulse_left_edge.push_back(left);
      pulse_right_edge.push_back(right);

    } // if statement
  }
  // getchar();
  if (!trig) {
    event_charge_ten.push_back(temp_ten_charge);
    event_charge.push_back(temp_charge);
  }
}
// Find the baseline
double baseline_rms(vector<float> &v, vector<float> &sample, double *irms) {
  double rms = 0;
  double temp_base = 0;
  // double baseline_samples = accumulate(v.begin(),v.end(),0);
  double baseline_samples = 0;
  for (int k = 0; k < (int)v.size(); k++) {
    baseline_samples += v[k];
  }
  // if (signal_start&&debug_mode)
  //    cout<<" baseline_samfloats is  : "<<baseline_samples<<endl;
  baseline_samples /= (float)v.size();
  for (int i = 0; i < (int)v.size(); i++) {
    //    if (signal_start&&debug_mode)
    //        cout<<" baseline sample is  : "<<v[i]<<endl;
    rms += pow(v[i] - baseline_samples, 2);
  }
  rms = sqrt(rms / v.size());
  if (signal_start && debug_mode) {
    cout << " rms is : " << rms << " baseline is : " << baseline_samples
         << endl;
    getchar();
  }
  if (smoothing)
    extract_event(sample, baseline_samples, rms, sample.size(), 0, false);

  //}
  irms[0] = rms;
  return baseline_samples;
}
// Trigger analysis
double Trigger_info(vector<float> waveform) {
  double rms_trigger = 0;
  double base_trigger = baseline_rms(
      trigbaselinev, waveform,
      &rms_trigger); // Calculate baseline and rms then pass to pulse finder
  extract_event(waveform, base_trigger, rms_trigger, number_of_samples, 0,
                true);
  double time;
  baselinev.clear();
  if (pulse_left_edge.size() == 0) {
    time = 0;
    cout << " Can not find the trigger pulse for this event ! " << endl;
  } else {
    time = (pulse_left_edge[0]);
  }
  pulse_left_edge.clear();
  pulse_right_edge.clear();
  trigbaselinev.clear();
  // waveform.clear();
  return time;
}

void getwaveform(vector<float> &v, int evt, PyArrayObject *arr, float mult = 1,
                 bool trig = false) {
  /*if(current_sweep%100000==0){
          cout<<"starting at "<<starti<<" in buffer"<<endl;
          cout<<"Evt "<<(current_sweep-numread)<<" of this buffer"<<endl;
  }*/
  double datum;
  for (int i = 0; i < number_of_samples; i++) {
    // cout << "reading sample " << i << endl;
    datum = *reinterpret_cast<int *>(PyArray_GETPTR2(arr, evt, i));
    datum *= mult / (double)adc_per_Volt;
    v.push_back(datum);
    if (i < baseline_samples_set) {
      if (!trig)
        baselinev.push_back(datum);
      else
        trigbaselinev.push_back(datum);
    }
  }
}
int calcnumchannels(int mask) {
  int numchans = 0;
  while (mask > 0) {
    numchans += mask % 2;
    mask /= 2;
  }
  return numchans;
}
vector<float> Smoothen(vector<float> &v) {
  vector<float> smoothv;

  double sum = 0.0;
  double movingAverage = 0.0;
  int WidthSample = 5;
  int rollingSum = 0.0;

  if (boxsmoothing) {                        // boxcar smoothing
    for (int it = 0; it < iteration; it++) { // how many passes
      if (it < 1) {
        for (int i = 0; i < (int)v.size(); i++) {
          sum = 0.0;
          movingAverage = 0.0;
          if (i < (MovingWindowSize - 1) / 2) {
            smoothv.push_back(v[i]);
          } else if (i >= (MovingWindowSize - 1) / 2 &&
                     i < ((int)v.size() - MovingWindowSize)) {
            for (int j = (i - (MovingWindowSize - 1) / 2);
                 j <= (i + (MovingWindowSize - 1) / 2); ++j) {
              sum += v[j];
            }
            movingAverage = sum / MovingWindowSize;
            smoothv.push_back(movingAverage);
          } else if (i >= ((int)v.size() - MovingWindowSize)) {
            smoothv.push_back(v[i]);
          }
        }
      } else {
        for (int i = 0; i < (int)smoothv.size(); i++) {
          sum = 0.0;
          movingAverage = 0.0;
          if (i >= (MovingWindowSize - 1) / 2 &&
              i < ((int)v.size() - MovingWindowSize)) {
            for (int j = (i - (MovingWindowSize - 1) / 2);
                 j <= (i + (MovingWindowSize - 1) / 2); ++j) {
              sum += v[j];
            }
            movingAverage = sum / MovingWindowSize;
            smoothv[i] = movingAverage;
          }
        }
      }
    }
    for (int i = 0; i < baseline_samples_set; i++) {
      smoothedBaseLine.push_back(smoothv[i]);
    }
  }
  if (triangle) { // triangular smoothing which preserves area under peaks ( 3
                  // points)
    for (int it = 0; it < iteration; it++) { // how many passes
      if (it < 1) {
        for (int i = 0; i < (int)v.size(); i++) {
          sum = 0.0;
          movingAverage = 0.0;
          if (i < (MovingWindowSize - 1) / 2) {
            smoothv.push_back(v[i]);
          } else if (i >= (MovingWindowSize - 1) / 2 &&
                     i < ((int)v.size() - MovingWindowSize)) {
            for (int j = (i - (MovingWindowSize - 1) / 2);
                 j <= (i + (MovingWindowSize - 1) / 2); ++j) {
              if (j == (i - (MovingWindowSize - 1) / 2) ||
                  j == (i + (MovingWindowSize - 1) / 2)) {
                sum += v[j];
              }
              if (j == i) {
                sum += 2 * v[j];
              }
            }
            movingAverage = sum / 4;
            smoothv.push_back(movingAverage);
          } else if (i >= ((int)v.size() - MovingWindowSize)) {
            smoothv.push_back(v[i]);
          }
        }
      } else {
        for (int i = 0; i < (int)smoothv.size(); i++) {
          sum = 0.0;
          movingAverage = 0.0;
          if (i >= (MovingWindowSize - 1) / 2 &&
              i < ((int)v.size() - MovingWindowSize)) {
            for (int j = (i - (MovingWindowSize - 1) / 2);
                 j <= (i + (MovingWindowSize - 1) / 2); ++j) {
              if (j == (i - (MovingWindowSize - 1) / 2) ||
                  j == (i + (MovingWindowSize - 1) / 2)) {
                sum += smoothv[j];
              }
              if (j == i) {
                sum += 2 * smoothv[j];
              }
            }
            movingAverage = sum / 4;
            smoothv[i] = movingAverage;
          }
        }
      }
    }
    for (int i = 0; i < baseline_samples_set; i++) {
      smoothedBaseLine.push_back(smoothv[i]);
    }
  }
  if (rolling) {
    smoothv.resize(8000, 0);
    for (int n = 0; n < WidthSample; n++) {
      rollingSum += v[n];
    }
    for (int n = 0; n < WidthSample; n++) {
      smoothv[n] = rollingSum / WidthSample;
    }
    for (int n = WidthSample; n < (int)v.size(); n++) {
      rollingSum += v[n];
      rollingSum -= v[n - WidthSample];
      smoothv[n - WidthSample / 2] = rollingSum / WidthSample;
    }
    for (int n = (int)v.size() - WidthSample / 2; n < (int)v.size(); n++) {
      smoothv[n] = smoothv[(int)v.size() - WidthSample / 2 - 1];
    }
    for (int i = 0; i < baseline_samples_set; i++) {
      smoothedBaseLine.push_back(smoothv[i]);
    }
  }
  if (debug_mode)
    cout << " Smooth vector has elements : " << smoothv.size() << endl;
  return smoothv;
}

static void show_usage(string name) {
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
int main(int argc, char *argv[]) {
  string filename;
  string outfilename;
  string working_dir;
  string baseline_file;
  std::string pydir;

  bool use_trigger = false;
  bool trigger_inversion = false;
  bool invert_waveform = false;

  int num_sams = 0;
  int trig_channel;
  int wform_channel = 1;

  if (argc < 4) {
    show_usage(argv[0]);
    return 1;
  }
  for (int i = 1; i < argc; ++i) {
    string arg = argv[i];
    if ((arg == "-h") || (arg == "--help")) {
      show_usage(argv[0]);
      return 0;
    } else if (arg == "-wd") {
      working_dir = argv[i + 1];
    } else if (arg == "-i") {
      filename = argv[i + 1];
    } else if (arg == "-o") {
      outfilename = argv[i + 1];
    } else if (arg == "-t") {
      trig_channel = atoi(argv[i + 1]);
      use_trigger = true;
    } else if (arg == "-wform") {
      wform_channel = atoi(argv[i + 1]);
    } else if (arg == "-bs") {
      baseline_samples_set = atoi(argv[i + 1]);
    } else if (arg == "-bf") {
      baseline_file = argv[i + 1];
      use_basefile = true;
    } else if (arg == "-pt") {
      pulseThresh = atof(argv[i + 1]);
    } else if (arg == "-win") {
      windowSize = atof(argv[i + 1]);
    } else if (arg == "-trigger") {
      trigger_inversion = true;
    } else if (arg == "-invert") {
      invert_waveform = true;
    } else if (arg == "-sams") {
      lim_sams = true;
      num_sams = atof(argv[i + 1]);
    } else if (arg == "-debug") {
      debug_mode = true;
    } else if (arg == "-s") {
      smoothing = true;
    } else if (arg == "-box") {
      boxsmoothing = true;
      smoothing = true;
    } else if (arg == "-mwin") {
      MovingWindowSize = atoi(argv[i + 1]);
    } else if (arg == "-roll") {
      rolling = true;
      smoothing = true;
    } else if (arg == "-tri") {
      triangle = true;
      smoothing = true;
    } else if (arg == "-sit") {
      iteration = atoi(argv[i + 1]);
    } else if (arg == "-pydir") {
      pydir = argv[i + 1];
      // pydir.append("/readlogs.py");
    }
  }

  double fixedbase;
  double fixedrms;
  if (use_basefile) {
    TFile *basefile = new TFile(baseline_file.c_str(), "READ");
    if (basefile == NULL) {
      cout << "Couldn't open baseline file | Using standard methods instead"
           << endl;
      use_basefile = false;
    } else {
      TVectorD *fbase = (TVectorD *)basefile->Get("baseline");
      TVectorD *frms = (TVectorD *)basefile->Get("rms");
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
  /*
    fin.open(open_filename, ios::binary | ios::in | ios::ate);

    if (fin.is_open()) {
      // memblock.resize(size);
      size = fin.tellg();
      fin.seekg(0, ios::beg);
      // cout<<"size: "<<size<<endl;
      // numevts = new char [5];
      fin.read((char *)&Nevts, sizeof(Nevts));
      cout << Nevts << " Events" << endl;

      fin.read((char *)&number_of_samples, sizeof(number_of_samples));
      cout << number_of_samples << " Samples" << endl;

      fin.read((char *)&mask, sizeof(mask));
      cout << "Mask: " << mask << endl;
      Nchannels = calcnumchannels(mask);
      cout << Nchannels << " channels" << endl;

      fin.read((char *)&dummy, sizeof(dummy));
    } else {
      cout << "Failed to open file" << endl;
      return -1;
    }
    int predsize =
        Nevts * Nchannels * (2 * 4 + 2 * number_of_samples + 4) + 4 * 4;
    if (size < predsize) {
      cout << "Warning::Size predicted from header is greater than actual size"
           << endl;
      return -1;
    }

    if (wform_channel >= Nchannels || wform_channel < 0 ||
        ((trig_channel >= Nchannels || trig_channel < 0) && use_trigger)) {
      cout << "Channel numbers given do not match file header, double check your
    " "inputs."
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
  */

  std::cout << "Initializing python" << std::endl;
  Py_Initialize();
  std::cout << "Updating syspath with pydir: " << pydir << std::endl;
  PyObject *sysPath = PySys_GetObject("path");
  PyList_Append(sysPath, PyUnicode_FromFormat("%s", pydir.c_str()));

  std::cout << "Loading module" << std::endl;
  // Load module
  PyObject *pName = PyUnicode_FromString("ReadNpyfiles");
  PyObject *pModule = PyImport_Import(pName);
  PyArrayObject *pData = NULL;
  Py_DECREF(pName);
  int nd1 = 0, nd2 = 0;
  if (pModule != NULL) {
    std::cout << "Py Module Found" << std::endl;

    // Get function from module
    PyObject *pFunc = PyObject_GetAttrString(pModule, "loadnpyfile");
    if (pFunc && PyCallable_Check(pFunc)) {
      PyObject *pargs = PyTuple_New(1);
      PyObject *pval = NULL;
      pval = PyUnicode_FromFormat("%s", open_filename);
      PyTuple_SetItem(pargs, 0, pval);

      pval = PyObject_CallObject(pFunc, pargs);
      if (pval != NULL && PyTuple_Check(pval)) {
        // PyObject *pTemp = PyTuple_GetItem(pval, 0);
        pData = reinterpret_cast<PyArrayObject *>(PyTuple_GetItem(pval, 0));
        nd1 = PyLong_AsLong(PyTuple_GetItem(pval, 1));
        nd2 = PyLong_AsLong(PyTuple_GetItem(pval, 2));
        npy_intp *pshape = PyArray_SHAPE(pData);

        std::cout << "Array has " << PyArray_NDIM(pData) << " dims and "
                  << pshape[0] << " x " << pshape[1] << " Elements"
                  << std::endl;
        // cout << "Array is int? " << PyArray_ISINTEGER(pData) << endl;
        // Py_XDECREF(pTemp);
      } else {
        PyErr_Print();
        std::cout << "Failed to get result" << std::endl;
      }
      Py_XINCREF(pData);
      Py_XDECREF(pval);
      Py_DECREF(pargs);
    } else {
      if (PyErr_Occurred())
        PyErr_Print();
      std::cout << "Couldn't find loadnpyfile " << std::endl;
    }
    Py_XDECREF(pFunc);
    Py_XDECREF(pModule);
  } else {
    PyErr_Print();
    std::cout << "Failed to load module" << std::endl;
  }
  // cout << "Typechecking pData" << endl;
  if (pData == NULL) {
    std::cout << "Failed to read data file. Exiting" << std::endl;
    Py_FinalizeEx();
    return -1;
  }

  number_of_samples = nd2;
  Nevts = nd1;
  // Plots for debugging pulse finding algorithm
  TGraph *t11;
  TGraph *t22;
  TGraph *t33;
  TGraph *t44;
  TGraph *t55;

  // char linea[200];// temp char for line in the file
  // TString* buff=new TString();
  char out_filename[200]; // full output filename
  double rms_value;

  sprintf(out_filename, "%s/%s", working_dir.c_str(), outfilename.c_str());
  cout << " Out put filename is : " << out_filename << endl;
  TFile *fout = new TFile(out_filename, "RECREATE");

  TH1D *h = new TH1D(("ADC_sum_waveform" + filename).c_str(),
                     ("#font[132]{WFD " + filename + " SumWaveForm}").c_str(),
                     10000, 0, 10000);
  h->SetXTitle("#font[132]{Sample (2ns)}");
  h->GetXaxis()->SetLabelFont(132);
  h->GetYaxis()->SetLabelFont(132);
  // Tetsing the dark hit counter
  TH1F *dark_hits = new TH1F("dark_hits", "dark_hits", 100, -0.5, 99.5);
  // dark_hits->SetBit(TH1::kCanRebin);

  // Create Ntuple to store properties of pulses found by pulse finder
  TNtuple *pulse;
  if (use_trigger)
    pulse = new TNtuple(
        "pulse", "pulse",
        "pulseHeight:pulseRightEdge:pulseLeftEdge:pulseCharge:pulsePeakTime:"
        "CalibratedTime:baselinerms:windowratio:sweep:bigstep:pulseLength5:"
        "pulseLength25:pulseLength50:pulseLength75:pulseLength80:pulseLength90:"
        "pulseLength95:pulseLength99:triggerpulseHeight:triggerpulseWidth:"
        "triggerpulsePeakTime");
  else
    pulse = new TNtuple(
        "pulse", "pulse",
        "pulseHeight:pulseRightEdge:pulseLeftEdge:pulseCharge:pulsePeakTime:"
        "CalibratedTime:baselinerms:windowratio:sweep:bigstep:pulseLength5:"
        "pulseLength25:pulseLength50:pulseLength75:pulseLength80:pulseLength90:"
        "pulseLength95:pulseLength99");

  TNtuple *event =
      new TNtuple("event", "event", "charge:charge_frac:baseline:rms:npulses");
  TTree *wforms_tree = new TTree("waveforms", "Waveform Tree");
  float waveforms[8192];
  float trigger_t;

  wforms_tree->Branch("pmt_waveforms", &waveforms[0],
                      TString::Format("waveforms[%i]/F", number_of_samples));
  if (use_trigger)
    wforms_tree->Branch("trigger_time", &trigger_t, "trigger_t/F");
  // Store the waveform plot for debugging
  TCanvas *waveplot[100];
  vector<float> baseline_sweep;
  vector<float> trigwaveform;

  for (int sweep = 0; sweep < Nevts; sweep++) {
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
    if (sweep % 100 == 0) {
      cout << " This is sweep : " << sweep << endl;
      cout << " Trigger time: " << trigger_t << endl;
    }
    current_sweep = sweep;

    if (use_trigger) {
      signal_start = false;
      getwaveform(trigwaveform, sweep, pData, (trigger_inversion ? -1.0 : 1.0),
                  true);
      trigger_t = Trigger_info(trigwaveform);
      trigwaveform.clear();
    }

    signal_start = true;
    // cout << "reading Waveform" << endl;
    getwaveform(raw_waveform, sweep, pData, (invert_waveform ? -1.0 : 1.0));
    // cout << "Waveform read" << endl;

    std::copy(raw_waveform.begin(), raw_waveform.end(), waveforms);
    wforms_tree->Fill();
    number_of_peaks = 0.0;
    double thisbase = (use_basefile ? fixedrms : 0);

    if (smoothing) {
      smoothingv.clear();
      smoothingv = Smoothen(raw_waveform);
      rms_value = (use_basefile
                       ? fixedbase
                       : baseline_rms(smoothedBaseLine, smoothingv, &thisbase));
    } else {
      rms_value =
          (use_basefile ? fixedbase
                        : baseline_rms(baselinev, raw_waveform, &thisbase));
      extract_event(raw_waveform, rms_value, thisbase, number_of_samples,
                    (use_trigger ? trigger_t : 0));
    }

    if (debug_mode) {
      cout << " basline is  : " << rms_value << " rms is : " << thisbase
           << endl;
      getchar();
    }

    event_baseline.push_back(rms_value);
    event_rms.push_back(thisbase);
    // event_time.push_back(number_of_samples);
    baseline_sweep.push_back(
        rms_value); // save baseline for checking baseline shifting
    dark_hits->Fill(number_of_peaks);
    npeaks.push_back(number_of_peaks);
    baselinev.clear();
    // fill canvases
    if (sweep < 100) {
      t11 = new TGraph();
      t22 = new TGraph();
      t33 = new TGraph();
      t55 = new TGraph();
      for (int j = 0; j < (int)pulse_left_edge.size(); j++) {
        t22->SetPoint(j, pulse_left_edge[j], startv[j]);
        t33->SetPoint(j, pulse_right_edge[j], endv[j]);
      }
      for (int sam = 0; sam < number_of_samples; sam++) {
        t11->SetPoint(sam, sam, raw_waveform[sam]);
        if (smoothing) {
          if (sam < (int)smoothingv.size())
            t55->SetPoint(sam, sam, smoothingv[sam]);
        }
      }
      char plotname[30];
      sprintf(plotname, "waveform%d", sweep);
      waveplot[sweep] = new TCanvas(plotname);
      TLine *line = new TLine(0, rms_value, number_of_samples, rms_value);
      TLine *line2 = new TLine(0, rms_value - thisbase, number_of_samples,
                               rms_value - thisbase);
      TLine *line3 = new TLine(0, rms_value + thisbase, number_of_samples,
                               rms_value + thisbase);
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
      t55->SetLineColor(2);
      t11->Draw("alp");
      t22->Draw("p");
      t33->Draw("p");
      t55->Draw("lp");
      line->Draw("");
      line2->Draw("");
      line3->Draw("");
      waveplot[sweep]->Write();
      cout << " plotting the waveform, this is sweep : " << sweep << endl;
    }
    pulse_left_edge.clear();
    pulse_right_edge.clear();
    startv.clear();
    endv.clear();
    raw_waveform.clear();
    // skip += ;
  }

  Py_FinalizeEx();

  // Fill Ntuple
  // pulseHeight:pulseRightEdge:pulseLeftEdge:pulseCharge:pulsePeakTime
  cout << " before tree fill !" << endl;
  for (int i = 0; i < (int)amplitude.size(); i++) {
    if (use_trigger) {
      float mtinput[] = {amplitude[i],
                         pr[i],
                         pl[i],
                         charge_v[i],
                         amplitude_position[i],
                         CalibratedTime[i],
                         pulsebaseline_rms[i],
                         windowratio[i],
                         (float)event_n[i],
                         biggeststep[i],
                         pulse_length5[i],
                         pulse_length25[i],
                         pulse_length50[i],
                         pulse_length75[i],
                         pulse_length80[i],
                         pulse_length90[i],
                         pulse_length95[i],
                         pulse_length99[i],
                         triggerHeight[i],
                         triggerWidth[i],
                         triggerPosition[i]};
      float *mtinputpoint = mtinput;
      pulse->Fill(mtinputpoint);
    } else {
      float mtinput[] = {amplitude[i],
                         pr[i],
                         pl[i],
                         charge_v[i],
                         amplitude_position[i],
                         CalibratedTime[i],
                         pulsebaseline_rms[i],
                         windowratio[i],
                         (float)event_n[i],
                         biggeststep[i],
                         pulse_length5[i],
                         pulse_length25[i],
                         pulse_length50[i],
                         pulse_length75[i],
                         pulse_length80[i],
                         pulse_length90[i],
                         pulse_length95[i],
                         pulse_length99[i]};
      float *mtinputpoint = mtinput;
      pulse->Fill(mtinputpoint);
    }
  }
  cout << " after tree fill ! " << endl;
  TGraph *baseline_plot = new TGraph();
  for (int i = 0; i < (int)baseline_sweep.size(); i++) {
    baseline_plot->SetPoint(i, i, baseline_sweep[i]);
  }
  for (int i = 0; i < (int)event_charge.size(); i++) {
    event->Fill(event_charge[i], event_charge_ten[i], event_baseline[i],
                event_rms[i], npeaks[i]);
  }

  // Baseline plot
  TCanvas *bplot = new TCanvas("bplot", "bplot");
  baseline_plot->SetMarkerStyle(22);
  baseline_plot->Draw("AP");

  TObject *nosinfo = new TObject();
  nosinfo->SetUniqueID(number_of_samples);

  cout << " Total sweeps is : " << Nevts << endl;

  nosinfo->Write("Nsamples");
  wforms_tree->Write();
  pulse->Write();
  event->Write();
  bplot->Write();
  fout->Write();
  fout->Close();

  return 0;
}
