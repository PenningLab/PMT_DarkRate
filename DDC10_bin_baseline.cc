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
bool debug_mode = false;

int current_sweep = 0;
int adc_per_Volt = 8192;

vector<float> raw_waveform;
vector<double> event_baseline;
vector<double> event_rms;

int number_of_samples;
int Nchannels;
int Nevts;
short int *buff;

// Find the baseline
double baseline_rms(vector<float> &v){
    double rms=0;
    double temp_base = 0;
    //double baseline_samples = accumulate(v.begin(),v.end(),0);
    double baseline_samples=0;
    for (int k=0;k<v.size();k++){
	       baseline_samples += v[k];
    }
    if (debug_mode)
        cout<<" baseline_samples is  : "<<baseline_samples<<endl;
    baseline_samples /= v.size();
    for (int i=0;i<v.size();i++){
        if (debug_mode)
            cout<<" baseline sample is  : "<<v[i]<<endl;
        rms += pow(v[i] - baseline_samples,2);
    }
    rms = sqrt(rms / v.size());

    event_baseline.push_back(baseline_samples);
    event_rms.push_back(rms);
    return baseline_samples;
}

void getwaveform(vector<float> &v,int channel,int numread,float mult=1){
	int starti = ((current_sweep-numread)*Nchannels + channel)*(4 + 2 + number_of_samples) + 4;
	/*if(current_sweep%100000==0){
		cout<<"starting at "<<starti<<" in buffer"<<endl;
		cout<<"Evt "<<(current_sweep-numread)<<" of this buffer"<<endl;
	}*/
	double datum;
	for(int i=0;i<number_of_samples;i++){
		datum = (double)buff[i+starti]*mult/(double)adc_per_Volt;
		v.push_back(datum);
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
	<<" -bo : Name of output text file.\n"
    <<" -i : Name of input file.\n"
    <<" -wd : Working directory\n"
	<<" -wform : Waveform channel\n"
	<<" -invert: invert waveform\n"
    <<" -debug : Get in the debugging mode.\n"
    <<" -h or --help : Show the usage\n"
    <<" Enjoy ! -Ryan Wang"<<endl;
}
int main(int argc, char *argv[]){
	string filename;
	string outfilename;
	string baseoutfilename;
	string working_dir;

	bool invert_waveform = false;

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
		else if (arg=="-bo"){
            baseoutfilename = argv[i+1];
        }
		else if (arg=="-wform"){
            wform_channel = atoi(argv[i+1]);
        }
		else if (arg=="-invert"){
			invert_waveform = true;
		}
        else if (arg=="-debug"){
            debug_mode = true;
        }
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

	if(wform_channel>=Nchannels || wform_channel<0){
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

    TNtuple *event = new TNtuple("event","event","baseline:rms");
	TTree *wforms_tree = new TTree("waveforms","Waveform Tree");
	float waveforms[8192];
	float trigger_t;

	wforms_tree->Branch("pmt_waveforms",&waveforms[0],TString::Format("waveforms[%i]/F",number_of_samples));

    // Store the waveform plot for debugging
    TCanvas *waveplot[100];

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

		getwaveform(raw_waveform,wform_channel,lastadd,(invert_waveform ? -1.0 : 1.0));

		std::copy(raw_waveform.begin(),raw_waveform.end(),waveforms);
		wforms_tree->Fill();
		rms_value = baseline_rms(raw_waveform);
		//fill canvases
		if(sweep<100){
			t11 = new TGraph();
			for (int sam=0;sam<number_of_samples;sam++){
				t11->SetPoint(sam,sam,raw_waveform[sam]);
			}
			char plotname [30];
			sprintf(plotname,"waveform%d",sweep);
			waveplot[sweep] = new TCanvas(plotname);
            TLine* line = new TLine(0,rms_value,number_of_samples,rms_value);
            line->SetLineColor(3);
            line->SetLineStyle(3);
            line->SetLineWidth(3);
            t11->Draw("alp");
            line->Draw("");
            waveplot[sweep]->Write();
            cout<<" plotting the waveform, this is sweep : "<<sweep<<endl;
		}
		raw_waveform.clear();
		//skip += ;
	}

    TGraph* baseline_plot = new TGraph();

	double mean_baseline = 0;
	double mean_rms = 0;

    for (int i=0;i<event_baseline.size();i++){
		event->Fill(event_baseline[i],event_rms[i]);
		baseline_plot->SetPoint(i,i,event_baseline[i]);
		mean_baseline += event_baseline[i];
		mean_rms += event_rms[i];
    }
	mean_baseline /= (double)event_baseline.size();
	mean_rms /= (double)event_rms.size();

    //Baseline plot
    TCanvas* bplot = new TCanvas("bplot","bplot");
    baseline_plot->SetMarkerStyle(22);
    baseline_plot->Draw("AP");

	TObject *nosinfo = new TObject();
	nosinfo->SetUniqueID(number_of_samples);

    cout<<" Total sweeps is : "<<Nevts<<endl;

	nosinfo->Write("Nsamples");
	wforms_tree->Write();
	event->Write();
    bplot->Write();
    fout->Write();
    fout->Close();

	ofstream baseout;
	char base_filename[200];// full input filename
	sprintf(base_filename,"%s/%s",working_dir.c_str(),baseoutfilename.c_str());
	baseout.open(base_filename);
	baseout << mean_baseline << " " << mean_rms << endl;
	baseout.close();
    return 0;
}
