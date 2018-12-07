/////////////////////////////////////////////////////////////////////////////////////////
// To compile : g++ -I /usr/common/usg/software/ROOT/5.34.20/include/ `root-config --cflags --libs` -o DDC10_data_readout DDC10_data_readout.cc
// To execute (help infomation gives detail utility) : ./DDC10_data_readout -h
/* Revision log :
 *
	4/25/2018 RW : Code for reading scope's text file and do the pulse finding.
	7/23/2018 RW : Add option for user to identify each item easily.



*/
/////////////////////////////////////////////////////////////////////////////////////////

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
#include <TObject.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TVectorD.h>
#include <TMatrix.h>
#include <TH2Poly.h>
using namespace std;
ifstream fin;
bool signal_start = false;
bool debug_mode = false;
//define pulse finding parameters
double pulseThresh = 8.0 ;
double tpulseThresh = 8.0 ;
double windowSize = 3.0;
double edgeThresh = 3.0;
double lookforward = 3.0;
int current_sweep = 0;
double number_of_peaks = 0.0;
int adc_per_Volt = 8192;
double resistance = 0.005;
int baseline_samples_set = 160;

//vector<int> run_info;
vector<float> raw_waveform;
vector<double> amplitude;
vector<double> charge_v;
vector<double> startv;
vector<double> endv;
vector<double> start_pointv;
vector<double> end_pointv;
vector<double> baseline;
vector<double> baselinev;
vector<double> trigbaselinev;
vector<double> pulse_left_edge;
vector<double> pulse_right_edge;
vector<double> amplitude_position;
vector<double> pl;
vector<double> pr;
vector<float> npeaks;
vector<int> cand_pulses; //indexes of windows which pass first cut
//vector<double> trigger_time;
vector<double> CalibratedTime;
vector<double> windowratio;
vector<double> pulsebaseline_rms;
vector<short int> event_n;

vector<double> event_charge;
vector<double> event_charge_ten;
vector<double> event_baseline;
vector<double> event_rms;

int number_of_samples;
int Nchannels;
int Nevts;
short int *buff;

// Simpson Integral
double SimpsIntegral(const vector<float>& samples, double baseline,int start, int end){
	int len;
	double qsum = 0.0;
	if ((end - start) % 2 == 0){
	/* If there are an even number of samples, then there are an odd
	number of intervals; but Simpson's rule works only on an even
	number of intervals. Therefore we use Simpson's method on the
	all but the final sample, and integrate the last interval
	using the trapezoidal rule */
		len = end - start - 1;
		qsum += (samples[end-1] + samples[end-2] - 2*baseline)/2.0;
	}
	else
		len = end - start;

	double qsimps;
	qsimps = samples[start]-baseline;
	for (int i=start; i<start+len; i+=2)
		qsimps += (samples[i]-baseline)*4;
	for (int i=start+1; i<len+start-1; i+=2)
		qsimps += (samples[i]-baseline)*2;
	qsimps += samples[start+len-1]-baseline;
	qsimps /= 3.0;

	qsum += qsimps;
	return qsum;
}


// Pulse Finding
void extract_event(vector<float> &v, double b ,double rms,int nos,int trigger,bool trig=false){
	double pThresh = (trig ? tpulseThresh : pulseThresh) * rms * windowSize;
	double eThresh = edgeThresh * rms * windowSize;
	double temp_charge = 0;
	double temp_ten_charge = 0;
	double window1 = 0;
	double window2 = 0;
	double window3 = 0;
	//cout<<" vector size is : "<<v.size()<<endl;
	//getchar();
	//Let's looking for the Pulses
	for (int i=0;i<nos;i++){
		double integral = SimpsIntegral(v,b,i,i+windowSize);
		int left = 0;
		int right = 0;
		int temp_peak = 0;
		double temp_startv = 0;
		double temp_endv = 0;

		window1 = window2;
		window2 = window3;
		window3 = accumulate(v.begin()+i,v.begin()+i+windowSize,0.0) - b*windowSize;

		if (integral>pThresh && window1<window2 && window3<window2){
			double iw1= 1.0e9;
			//cout<<" This is sample : "<<i<<" integral value is : "<<integral<<" pThresh is : "<<pThresh<<endl;
			left = i;
			integral = 1.0e9;
			while (integral > eThresh && left > windowSize && iw1>=integral){
				left --;
				iw1 = integral;
				integral = SimpsIntegral(v,b,left,left+windowSize);
				temp_startv = v[left];

			}

			integral = 1.0e9;
			right = i + windowSize;
			double thischarge = 0;
			bool end=false;
			iw1 = 1.0e9;
			while (!end){
				while (integral > eThresh && right<number_of_samples-1 && iw1>=integral){
					right ++;
					iw1 = integral;
					integral = SimpsIntegral(v,b,right-windowSize,right);
					temp_endv = v[right];

				}
				end = true ;
				int r = right;
				while (r < fmin((int)v.size()-1, right+lookforward)){
					r++ ;
					integral = SimpsIntegral(v,b,r-windowSize,r);
					if (integral > pThresh){
						right = r;
						end = false;
						break;
					}
				}
			}

			double max = -1.0e9;
			for (int j=left;j<right;j++){
				double s = v[j] - b;
				if (s > max){
					max = s;
					temp_peak = j;
				}
			}
			if (right >nos)
				continue;

			// noise veto
			float width = (right - left);
			float nwidth = 2*width;
			int nright = right + nwidth;
			int nleft = left - nwidth;
			if(nright>nos){
				nright = nos;
				nleft -= (nleft>nwidth ? nwidth : 0);
			}
			if(nleft<0){
				nleft = 0;
				nright += (nright<(nos-nwidth) ? nwidth : nos);
			}
			double dratio = SimpsIntegral(v,b,nleft,nright)/( ( (nwidth+1)/width )*SimpsIntegral(v,b,left,right) );

			//cout<<" Peak is : "<<temp_peak<<" max is : "<<max<<endl;
			if (signal_start){
				amplitude.push_back(max);
				amplitude_position.push_back(temp_peak);
				pl.push_back(left);
				pr.push_back(right);
				charge_v.push_back(SimpsIntegral(v,b,left,right)/resistance);
				startv.push_back(temp_startv);
				endv.push_back(temp_endv);
				CalibratedTime.push_back(temp_peak-trigger);
				windowratio.push_back(dratio);
				pulsebaseline_rms.push_back(rms);
				event_n.push_back(current_sweep);
				temp_charge +=SimpsIntegral(v,b,left,right)/resistance;
				if (i<300)
					temp_ten_charge += SimpsIntegral(v,b,left,right)/resistance;
			}

			pulse_left_edge.push_back(left);
			pulse_right_edge.push_back(right);

			//if (temp_peak>0 &&temp_peak<8000){
			//if (thischarge<1.0)
			if(trig) break;
			number_of_peaks ++;
			//}
			//cout<<" This is sample : "<<i<<" Charge integral is : "<<SimpsIntegral(v,b,left,right)<<endl;

			if (SimpsIntegral(v,b,left,right) <= eThresh){
				i = right -1;
				continue;
			}
			i = right;
		}//if statement

	}
	//getchar();
	if (!trig){
		event_charge_ten.push_back(temp_ten_charge);
		event_charge.push_back(temp_charge);
	}
}

// Find the baseline
double baseline_rms(vector<double> &v, vector<float> &sample,double rms){
	rms = 0;
	double temp_base = 0;
	//double baseline_samples = accumulate(v.begin(),v.end(),0);
	double baseline_samples=0;
	for (int k=0;k<v.size();k++){
		   baseline_samples += v[k];
	}
	if (signal_start&&debug_mode)
		cout<<" baseline_samples is  : "<<baseline_samples<<endl;
	baseline_samples /= v.size();
	for (int i=0;i<v.size();i++){
		if (signal_start&&debug_mode)
			cout<<" baseline sample is  : "<<v[i]<<endl;
		rms += pow(v[i] - baseline_samples,2);
	}
	rms = sqrt(rms / v.size());
	if (signal_start&&debug_mode){
		cout<<" rms is : "<<rms<<" baseline is : "<<baseline_samples<<endl;
		getchar();
	}

	return baseline_samples;
}
//Trigger analysis
double Trigger_info(vector<float> waveform){
	double rms_trigger;
	double base_trigger = baseline_rms(trigbaselinev,waveform,rms_trigger);// Calculate baseline and rms then pass to pulse finder
	extract_event(waveform,base_trigger,rms_trigger,number_of_samples,0,true);
	double time;
	baselinev.clear();
	if (pulse_left_edge.size() == 0){
		time = 0;
		cout<<" Can not find the trigger pulse for this event ! "<<endl;
	}
	else{
		time = (pulse_left_edge[0]);
	}
	pulse_left_edge.clear();
	pulse_right_edge.clear();
	trigbaselinev.clear();
	//waveform.clear();
	return time;
}

void getwaveform(vector<float> &v,int channel,int numread,float mult=1,bool trig=false){
	int starti = ((current_sweep-numread)*Nchannels + channel)*(4 + 2 + number_of_samples) + 4;
	/*if(current_sweep%100000==0){
		cout<<"starting at "<<starti<<" in buffer"<<endl;
		cout<<"Evt "<<(current_sweep-numread)<<" of this buffer"<<endl;
	}*/
	double datum;
	for(int i=0;i<number_of_samples;i++){
		datum = (double)buff[i+starti]*mult/(double)adc_per_Volt;
		v.push_back(datum);
		if (i<baseline_samples_set){
			if(!trig) baselinev.push_back(datum);
			else trigbaselinev.push_back(datum);
		}
	}
}
int calcnumchannels(int mask){
	int numchans=0;
	while (mask>0){
		numchans += mask%2;
		mask /= 2;
	}
	return numchans;
}


static void show_usage(string name){
	cout<<" Usage : ./DDC10_data_readout [-co] file1 "<<name<<" Options:\n"
	<<" -o : Name of output file.\n"
	<<" -i : Name of input file.\n"
	<<" -wd : Working directory\n"
	<<" -wform : Waveform channel\n"
	<<" -t : trigger channel\n"
	<<" -invert: invert waveform\n"
	<<" -trigger : invert trigger pulse \n"
	<<" -pt : pulse threshold\n"
	<<" -e : write event tree\n"
	<<" -debug : Get in the debugging mode.\n"
	<<" -h or --help : Show the usage\n"
	<<" Enjoy ! -Ryan Wang"<<endl;
}
int main(int argc, char *argv[]){
	string filename;
	string outfilename;
	string working_dir;
	string baseline_file;

	bool use_trigger = false;
	bool trigger_inversion = false;
	bool invert_waveform = false;
	bool use_basefile = false;

	int trig_channel;
	int wform_channel=1;

	if (argc<4){
		show_usage(argv[0]);
		return 1;
	}
	for (int i=1;i<argc;++i){
		string arg = argv[i];
		if ((arg=="-h") || (arg=="--help")){
			show_usage(argv[0]);
			return 0;
		}
		else if (arg=="-wd"){
			working_dir = argv[i+1] ;
		}
		else if (arg=="-i") {
			filename = argv[i+1];
		}
		else if (arg=="-o"){
			outfilename = argv[i+1];
		}
		else if (arg=="-t"){
			trig_channel = atoi(argv[i+1]);
			use_trigger = true;
		}
		else if (arg=="-wform"){
			wform_channel = atoi(argv[i+1]);
		}
		else if (arg=="-bs"){
			baseline_samples_set = atoi(argv[i+1]);
		}
		else if (arg=="-bf"){
			baseline_file = argv[i+1];
			use_basefile = true;
		}
		else if (arg=="-pt"){
				pulseThresh = atof(argv[i+1]);
		}
		else if (arg=="-win"){
				windowSize = atof(argv[i+1]);
		}
		else if (arg=="-trigger"){
			trigger_inversion = true;
		}
		else if (arg=="-invert"){
			invert_waveform = true;
		}
		else if (arg=="-debug"){
			debug_mode = true;
		}
	}

	double fixedbase;
	double fixedrms;
	if(use_basefile){
		TFile *basefile = new TFile(baseline_file.c_str(),"READ");
		if(basefile == NULL){
			cout<<"Couldn't open baseline file | Using standard methods instead"<<endl;
			use_basefile = false;
		}
		else{
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
	int mask,size;
	char open_filename[200];// full input filename
	sprintf(open_filename,"%s/%s",working_dir.c_str(),filename.c_str());
	cout<<" We are opening "<<open_filename<<endl;
	fin.open(open_filename,ios::binary|ios::in|ios::ate);

	if (fin.is_open()){
		//memblock.resize(size);
		size = fin.tellg();
		fin.seekg(0,ios::beg);
		//cout<<"size: "<<size<<endl;
		//numevts = new char [5];
		fin.read((char*)&Nevts,sizeof(Nevts));
		cout<<Nevts<<" Events"<<endl;

		fin.read((char*)&number_of_samples,sizeof(number_of_samples));
		cout<<number_of_samples<<" Samples"<<endl;

		fin.read((char*)&mask,sizeof(mask));
		cout<<"Mask: "<<mask<<endl;
		Nchannels = calcnumchannels(mask);
		cout<<Nchannels<<" channels"<<endl;

		fin.read((char*)&dummy,sizeof(dummy));
	}
	else{
		cout<<"Failed to open file"<<endl;
		return -1;
	}
	int predsize = Nevts*Nchannels*(2*4 + 2*number_of_samples + 4) + 4*4;
	if(size<predsize){
		cout<<"Warning::Size predicted from header is greater than actual size"<<endl;
		return -1;
	}

	if(wform_channel>=Nchannels || wform_channel<0 || ((trig_channel>=Nchannels || trig_channel<0) && use_trigger)){
		cout<<"Channel numbers given do not match file header, double check your inputs."<<endl;
		return -1;
	}
	int evtsize = Nchannels*(2*4 + 2*number_of_samples + 4);
	int buffsize;
	if((Nevts*evtsize)>20971520) buffsize = 20971520/evtsize;
	else buffsize = Nevts;
	cout<<"Using buffer of "<<buffsize<<" events"<<endl;

	// Plots for debugging pulse finding algorithm
	TGraph* t11;
	TGraph* t22;
	TGraph* t33;
	TGraph* t44;

	//char linea[200];// temp char for line in the file
	//TString* buff=new TString();
	char out_filename[200]; // full output filename
	int iana = 0; // TGraph counter
	int pcount=0; // Pulse counter
	double temp_sum=0;
	double rms_value;


	sprintf(out_filename,"%s/%s",working_dir.c_str(),outfilename.c_str());
	cout<<" Out put filename is : "<<out_filename<<endl;
	TFile* fout = new TFile(out_filename,"RECREATE");

	TH1D* h = new TH1D(("ADC_sum_waveform" + filename).c_str(),
					 ("#font[132]{WFD "+filename+" SumWaveForm}").c_str(),
						   10000, 0, 10000);
	h->SetXTitle("#font[132]{Sample (2ns)}");
	h->GetXaxis()->SetLabelFont(132);
	h->GetYaxis()->SetLabelFont(132);
	//Tetsing the dark hit counter
	TH1F* dark_hits = new TH1F("dark_hits","dark_hits",100,0,10);
	//dark_hits->SetBit(TH1::kCanRebin);

	//Create Ntuple to store properties of pulses found by pulse finder
	TNtuple *pulse = new TNtuple("pulse","pulse","pulseHeight:pulseRightEdge:pulseLeftEdge:pulseCharge:pulsePeakTime:CalibratedTime:baselinerms:windowratio:sweep");
	TNtuple *event = new TNtuple("event","event","charge:charge_frac:baseline:rms:npulses");
	TTree *wforms_tree = new TTree("waveforms","Waveform Tree");
	float waveforms[8192];
	float trigger_t;

	wforms_tree->Branch("pmt_waveforms",&waveforms[0],TString::Format("waveforms[%i]/F",number_of_samples));
	if (use_trigger)
		wforms_tree->Branch("trigger_time",&trigger_t,"trigger_t/F");
	// Store the waveform plot for debugging
	TCanvas *waveplot[100];
	vector<double> baseline_sweep;
	vector<float> trigwaveform;

	int skip=0;
	int readin=0;
	int lastadd=0;
	for(int sweep=0;sweep<Nevts;sweep++){
		if((readin)<(sweep+1)){
			int arrsize = buffsize*evtsize/2;
			if((Nevts-readin)<buffsize){
				arrsize = (Nevts-readin)*evtsize/2;
			}
			buff = new short int[arrsize];
			fin.read((char*)&buff[0],arrsize*sizeof(buff[0]));
			lastadd = readin;
			readin += arrsize*2/evtsize;

			cout<<"Read in "<<readin<<" events so far"<<endl;
		}
		if (sweep%1000==0){
			cout<<" This is sweep : "<<sweep<<endl;
		}
		current_sweep = sweep;

		if(use_trigger){
			signal_start = false;
			getwaveform(trigwaveform,trig_channel,lastadd,(trigger_inversion ? -1.0 : 1.0),true);
			trigger_t = Trigger_info(trigwaveform);
			trigwaveform.clear();
		}

		signal_start = true;
		getwaveform(raw_waveform,wform_channel,lastadd,(invert_waveform ? -1.0 : 1.0));

		std::copy(raw_waveform.begin(),raw_waveform.end(),waveforms);
		wforms_tree->Fill();
		number_of_peaks = 0.0;
		double thisbase = (use_basefile ? fixedrms : 0);
		rms_value = (use_basefile ? fixedbase : baseline_rms(baselinev,raw_waveform,thisbase));
		event_baseline.push_back(rms_value);
		event_rms.push_back(thisbase);
		baseline_sweep.push_back(rms_value);// save baseline for checking baseline shifting

		extract_event(raw_waveform,rms_value,thisbase,number_of_samples,(use_trigger ? trigger_t : 0));
		//second pass goes here

		dark_hits->Fill(number_of_peaks);
		npeaks.push_back(number_of_peaks);


		baselinev.clear();
		//fill canvases
		if(sweep<100){
			t11 = new TGraph();
			t22 = new TGraph();
			t33 = new TGraph();
			for (int j=0;j<pulse_left_edge.size();j++){
				t22->SetPoint(j,pulse_left_edge[j],startv[j]);
				t33->SetPoint(j,pulse_right_edge[j],endv[j]);
			}
			for (int sam=0;sam<number_of_samples;sam++){
				t11->SetPoint(sam,sam,raw_waveform[sam]);
			}
			char plotname [30];
			sprintf(plotname,"waveform%d",sweep);
			waveplot[sweep] = new TCanvas(plotname);
			TLine* line = new TLine(0,rms_value,number_of_samples,rms_value);
			TLine* line2 = new TLine(0,rms_value-thisbase,number_of_samples,rms_value-thisbase);
			TLine* line3 = new TLine(0,rms_value+thisbase,number_of_samples,rms_value+thisbase);
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
			waveplot[sweep]->Write();
			cout<<" plotting the waveform, this is sweep : "<<sweep<<endl;
		}
		for (int i=0;i<amplitude.size();i++){
		  pulse->Fill(amplitude[i],pr[i],pl[i],charge_v[i],amplitude_position[i],CalibratedTime[i],pulsebaseline_rms[i],windowratio[i],(float)event_n[i]);
		}
		event->Fill(event_charge[sweep],event_charge_ten[sweep],event_baseline[sweep],event_rms[sweep],npeaks[sweep]);

		amplitude.clear();
		pr.clear();
		pl.clear();
		charge_v.clear();
		amplitude_position.clear();
		CalibratedTime.clear();
		pulsebaseline_rms.clear();
		windowratio.clear();
		event_n.clear();

		pulse_left_edge.clear();
		pulse_right_edge.clear();
		startv.clear();
		endv.clear();
		raw_waveform.clear();
		//skip += ;
	}

	//Fill Ntuple
	//pulseHeight:pulseRightEdge:pulseLeftEdge:pulseCharge:pulsePeakTime

	TGraph* baseline_plot = new TGraph();
	for (int i=0;i<baseline_sweep.size();i++){
		baseline_plot->SetPoint(i,i,baseline_sweep[i]);
	}

	//Baseline plot
	TCanvas* bplot = new TCanvas("bplot","bplot");
	baseline_plot->SetMarkerStyle(22);
	baseline_plot->Draw("AP");

	TObject *nosinfo = new TObject();
	nosinfo->SetUniqueID(number_of_samples);

	cout<<" Total sweeps is : "<<Nevts<<endl;

	nosinfo->Write("Nsamples");
	wforms_tree->Write();
	pulse->Write();
	event->Write();
	bplot->Write();
	fout->Write();
	fout->Close();

	return 0;
}
