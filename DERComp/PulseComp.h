//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jun 21 14:14:14 2019 by ROOT version 6.16/00
// from TTree Data/LZ Events
// found on file: lz_201906192219_000010_000000_raw.root
//////////////////////////////////////////////////////////

#ifndef PulseComp_h
#define PulseComp_h

#include <TChain.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

// Headers needed by this particular selector
#include <vector>

class PulseComp : public TSelector
{
public:
	TTreeReader fReader; //! the tree reader
	TTree* fChain = 0; //! pointer to the analyzed TTree or TChain

	// Readers to access the data (delete the ones you do not need).
	TTreeReaderValue<UShort_t> evt = { fReader, "evt" };
	TTreeReaderValue<UShort_t> channel = { fReader, "channel" };
	TTreeReaderValue<UShort_t> hit = { fReader, "hit" };
	// TTreeReaderValue<ULong64_t> startTime = { fReader, "startTime" };
	TTreeReaderValue<UShort_t> nSamples = { fReader, "nSamples" };
	TTreeReaderArray<short> zData = { fReader, "zData" };
	TH1F* h_avgphdlg = 0;
	TH1F* h_avgphdhg = 0;
	Int_t numpulseshg = 0;
	int mynumpulseshg = 0;
	Int_t numpulseslg = 0;
	int mynumpulseslg = 0;

	// int number_of_samples;
	// define pulse finding parameters
	double pulseThresh = 6;
	double trigPulseThresh = 8.0;
	double windowSize = 4.0;
	double edgeThresh = 3.0;
	double lookforward = 3.0;
	int current_sweep = 0;
	double number_of_peaks = 0.0;
	int adc_per_Volt = 8192;
	double resistance = 0.005;
	int baseline_samples_set = 15;
	int promptwindow = 300;
	double pth = 6;
	double rms;
	double baseline;
	int pulse_left_edge;
	int pulse_right_edge;
	int pulse_peak_time;
	int npulses;

	PulseComp(TTree* /*tree*/ = 0)
	{
	}
	virtual ~PulseComp()
	{
	}
	virtual Int_t Version() const
	{
		return 2;
	}
	virtual void Begin(TTree* tree);
	virtual void SlaveBegin(TTree* tree);
	virtual void Init(TTree* tree);
	virtual Bool_t Notify();
	virtual Bool_t Process(Long64_t entry);
	virtual Int_t GetEntry(Long64_t entry, Int_t getall = 0)
	{
		return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0;
	}
	virtual void SetOption(const char* option)
	{
		fOption = option;
	}
	virtual void SetObject(TObject* obj)
	{
		fObject = obj;
	}
	virtual void SetInputList(TList* input)
	{
		fInput = input;
	}
	virtual TList* GetOutputList() const
	{
		return fOutput;
	}
	virtual void SlaveTerminate();
	virtual void Terminate();

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
	void extract_event(vector<float>& v, double b, double rms)
	{

		double pThresh = pulseThresh * rms * windowSize;
		double eThresh = edgeThresh * rms * windowSize;

		double temp_charge = 0;
		double temp_ten_charge = 0;
		double pulse_height_thresh = pth * rms;
		// cout<<" vector size is : "<<v.size()<<endl;
		// getchar();
		// Let's looking for the Pulses
		for (int i = baseline_samples_set; i < v.size() - windowSize; i++)
		{
			// std::cout << "Sample " << i << std::endl;
			double integral = SimpsIntegral(v, b, i, i + windowSize);
			int left = 0;
			int right = 0;
			int temp_peak = 0;
			double temp_startv = 0;
			double temp_endv = 0;

			double temp_bigstep = 0;

			if (integral > pThresh)
			{
				left = i;
				integral = 1.0e9;
				while (integral > eThresh && left > windowSize)
				{
					left--;
					integral = SimpsIntegral(v, b, left, left + windowSize);
					temp_startv = v[left];
				}

				integral = 1.0e9;
				right = i + windowSize;
				double thischarge = 0;
				bool end = false;
				while (!end)
				{
					while (integral > eThresh && right < *nSamples - 1)
					{
						right++;
						integral = SimpsIntegral(v, b, right - windowSize, right);
						temp_endv = v[right];
					}

					end = true;
					int r = right;
					while (r < fmin(fmin((int)v.size() - 1, *nSamples - 1), right + lookforward))
					{
						r++;
						integral = SimpsIntegral(v, b, r - windowSize, r);
						if (integral > pThresh)
						{
							right = r;
							end = false;
							break;
						}
					}
				}
				double max = -1.0e9;
				double totalq = SimpsIntegral(v, b, left, right);
				double tempq5 = 0, tempq25 = 0, tempq50 = 0, tempq75 = 0, tempq90 = 0, tempq95 = 0, tempq99 = 0;

				if (right > *nSamples)
					continue;

				if (SimpsIntegral(v, b, left, right) <= eThresh)
				{
					i = right - 1;
					continue;
				}
				i = right;
				for (int j = left; j < right; j++)
				{
					double s = v[j] - b;
					if (s > max)
					{
						max = s;
						temp_peak = j;
					}
				}
				if (max < pulse_height_thresh)
					continue;

				pulse_left_edge = left;
				pulse_right_edge = right;
				pulse_peak_time = temp_peak;
				npulses++;

			} // if statement
		}
	}
	// Find the baseline
	double baseline_rms(vector<float>& v, double* irms)
	{
		double rms = 0;
		double temp_base = 0;
		// double baseline_samples = accumulate(v.begin(),v.end(),0);
		double baseline_samples = 0;
		for (int k = 0; k < baseline_samples_set; k++)
		{
			baseline_samples += v[k];
		}
		// if (signal_start&&debug_mode)
		//    cout<<" baseline_samples is  : "<<baseline_samples<<endl;
		baseline_samples /= (double)baseline_samples_set;
		for (int i = 0; i < baseline_samples_set; i++)
		{
			//    if (signal_start&&debug_mode)
			//        cout<<" baseline sample is  : "<<v[i]<<endl;
			rms += pow(v[i] - baseline_samples, 2);
		}
		rms = sqrt(rms / (double)baseline_samples_set);

		//}
		irms[0] = rms;
		return baseline_samples;
	}
	bool getwaveform(vector<float>& v)
	{
		for (int i = 0; i < *nSamples; i++)
		{
			double datum = -(double)zData[i] / (double)adc_per_Volt;
			v.push_back(datum);
		}
		baseline = baseline_rms(v, &rms);
		npulses = 0;
		extract_event(v, baseline, rms);
		// std::cout << baseline << " baseline, " << rms << std::endl;
		if (npulses > 1)
			std::cout << npulses << " pulses" << std::endl;
		if (npulses == 1 && pulse_left_edge > baseline_samples_set)
			return true;
		else
			return false;
	}

	ClassDef(PulseComp, 0);
};

#endif

#ifdef PulseComp_cxx
void PulseComp::Init(TTree* tree)
{
	// The Init() function is called when the selector needs to initialize
	// a new tree or chain. Typically here the reader is initialized.
	// It is normally not necessary to make changes to the generated
	// code, but the routine can be extended by the user if needed.
	// Init() will be called many times when running on PROOF
	// (once per file to be processed).

	fReader.SetTree(tree);
}

Bool_t PulseComp::Notify()
{
	// The Notify() function is called when a new file is opened. This
	// can be either for a new TTree in a TChain or when when a new TTree
	// is started when using PROOF. It is normally not necessary to make changes
	// to the generated code, but the routine can be extended by the
	// user if needed. The return value is currently not used.

	return kTRUE;
}

#endif // #ifdef PulseComp_cxx
