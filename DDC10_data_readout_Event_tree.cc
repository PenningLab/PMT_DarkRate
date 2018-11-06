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
#include <TMatrixDSymEigen.h>
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
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TVectorD.h>
#include <TMatrix.h>
#include <TH2Poly.h>
//#include <boost/date_time.hpp>
using namespace std;
ifstream fin;
bool signal_start = false;
bool debug_mode = false;
//define pulse finding parameters
double pulseThresh =8.0 ;
double windowSize = 3.0;
double edgeThresh = 3.0;
double fix_threshold = 0.0028;
double lookforward = 3.0;
double number_of_peaks = 0.0;
int adc_per_Volt = 8192;
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
vector<double> pulse_left_edge;
vector<double> pulse_right_edge;
vector<double> amplitude_position;
vector<double> pl;
vector<double> pr;
vector<double> trigger_time;
vector<double> CalibratedTime;
vector<double> windowratio;

vector<double> event_charge;
vector<double> event_charge_ten;
vector<double> event_baseline;
vector<double> event_rms;
int number_of_samples;
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
void extract_event(vector<float> &v, double b ,double rms,int nos,int trigger){
    //double pThresh = pulseThresh * rms * sqrt(windowSize);
    double pThresh = fix_threshold;
    double eThresh = edgeThresh * rms * sqrt(windowSize);
    double temp_charge = 0;
    double temp_ten_charge = 0;
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

        if (integral > pThresh){
            //cout<<" This is sample : "<<i<<" integral value is : "<<integral<<" pThresh is : "<<pThresh<<endl;
            left = i;
            integral = 1.0e9;
            while (integral > eThresh && left > windowSize){
                left --;
                integral = SimpsIntegral(v,b,left,left+windowSize);
                temp_startv = v[left];

            }

            integral = 1.0e9;
            right = i + windowSize;
            bool end=false;
            while (!end){
                while (integral > eThresh && right<number_of_samples-1){
                    right ++;
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
			int width = right - left;
			int nright = right + width;
			int nleft = left - width;
			if(nright>nos){
				nright = nos;
				nleft -= (nleft>width ? width : 0);
			}
			if(nleft<0){
				nleft = 0;
				nright += (nright<(nos-width) ? width : nos);
			}
			double dratio = SimpsIntegral(v,b,nleft,nright)/SimpsIntegral(v,b,left,right);

            //cout<<" Peak is : "<<temp_peak<<" max is : "<<max<<endl;
            if (signal_start){
                amplitude.push_back(max);
                amplitude_position.push_back(temp_peak);
                pl.push_back(left);
                pr.push_back(right);
                charge_v.push_back(SimpsIntegral(v,b,left,right));
                startv.push_back(temp_startv);
                endv.push_back(temp_endv);
                CalibratedTime.push_back(temp_peak-trigger);
				windowratio.push_back(dratio);
                temp_charge +=SimpsIntegral(v,b,left,right);
                if (i<300)
                    temp_ten_charge += SimpsIntegral(v,b,left,right);
            }

            pulse_left_edge.push_back(left);
            pulse_right_edge.push_back(right);

            //if (temp_peak>0 &&temp_peak<8000){
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
    event_charge_ten.push_back(temp_ten_charge);
    event_charge.push_back(temp_charge);
}
// Find the baseline
double baseline_rms(vector<double> &v, vector<float> &sample, int nosamples,int ttime = 0){
    double rms=0;
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
    event_baseline.push_back(baseline_samples);
    event_rms.push_back(rms);
    extract_event(sample,baseline_samples,rms,nosamples,ttime);
    return baseline_samples;
}
/*void Trigger_info(string tfile,string wd){
    char trigger_path [120];
    sprintf(trigger_path,"%s/%s",wd.c_str(),tfile.c_str());
    cout<<" We are opening trugger file : "<<trigger_path<<endl;
    char linet[200];// temp char for line in the file
    int pc =0, sw =0;
    double rms_value_trigger = 0;
    TString* bufft=new TString();
    ifstream fint;
    fint.open(trigger_path,ios::in);
    cout<<" try to open file here"<<endl;
    if (!fint.is_open()){
        cout<<" Can't open the file "<<endl;
    }
    else{
        while(!fint.rdstate() ){
        //while(sweep < 682){
            fint.getline(linet,200);//>>datum;
            //cout<<" each line is : "<<linet<<endl;;
    		*bufft=linet;


            if (bufft->Contains("Event")){
                //cout<<" What we are seeing in the loop ? >> "<<linet<<endl;
                if (sw>0){
                    //cout<<" size baseline : "<<baselinev.size()<<endl;
                    rms_value_trigger = baseline_rms(baselinev,raw_waveform,number_of_samples,0);// Calculate baseline and rms then pass to pulse finder

                    baselinev.clear();
		            if (sw%1000==0){
                    	cout<<" This is sweep : "<<sw<<endl;
                    }
                    //if (pulse_left_edge[0] > 1000)
                    //cout<<" size is : "<<pulse_left_edge.size()<<endl;
                    //cout<<" pulse finding passed !, Trigger time is  : "<<pulse_left_edge[0]<<endl;
                    if (pulse_left_edge.size() == 0){
                        trigger_time.push_back(0);
                        cout<<" Can not find the trigger pulse for this event ! "<<endl;
                    }
                    else{
                        trigger_time.push_back(pulse_left_edge[0]);
                    }


                    pulse_left_edge.clear();
                    pulse_right_edge.clear();


                    raw_waveform.clear();
                }
                pc = 0;
                sw ++;
            }
            //cout<<"pcount is : "<<pcount<<endl;
            else{
                double datum = bufft->Atof();
                if (pc<number_of_samples+1){
                    //if (isinf(datum) || isnan(datum) || datum > 200000 || datum<-1){
                    //    datum=0;
                    //}
                    //cout<<" adc count : "<<datum<<endl;
                    datum /= adc_per_Volt;
                    //cout<<" volt : "<<datum<<endl;
                    raw_waveform.push_back(datum);
                    if (pc<baseline_samples_set){
                        baselinev.push_back(datum);
                    }

                }
                pc++;

            }




        }
    }
    fint.close();
}*/
static void show_usage(string name){
    cout<<" Usage : ./DDC10_data_readout [-co] file1 "<<name<<" Options:\n"
    <<" -o : Name of output file.\n"
    <<" -i : Name of input file.\n"
    <<" -n : Number of Samples in one event\n"
    <<" -wd : Working directory\n"
    <<" -t filename for trigger pulse\n"
    <<" -debug : Get in the debugging mode.\n"
    <<" -h or --help : Show the usage\n"
    <<" Enjoy ! -Ryan Wang"<<endl;
}
int main(int argc, char *argv[]){

    string filename;
    string outfilename;
    string triggerfilename;
    string working_dir;

    bool trigger_inversion = false;

    if (argc<2){
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
        else if (arg=="-n"){
            number_of_samples = atoi(argv[i+1]);
        }
        else if (arg=="-i") {
            filename = argv[i+1];
        }
        else if (arg=="-o"){
            outfilename = argv[i+1];
        }
        else if (arg=="-t"){
            triggerfilename = argv[i+1];
        }
        else if (arg=="-trigger"){
            trigger_inversion = true;
        }
        else if (arg=="-debug"){
            debug_mode = true;
        }
    }
    //cout<<" Before getting trigger info, trigger filename is : "<<triggerfilename<<endl;
    //Trigger_info(triggerfilename,working_dir);


    // Plots for debugging pulse finding algorithm
    TGraph* t11;
    TGraph* t22;
    TGraph* t33;
    TGraph* t44;

    char linea[200];// temp char for line in the file
    TString* buff=new TString();
    char open_filename[200];// full input filename
    char out_filename[200]; // full output filename
    int iana = 0; // TGraph counter
    int pcount=0; // Pulse counter
    int sweep=0; //event counter
    double temp_sum=0;
    double rms_value;

    sprintf(open_filename,"%s/%s",working_dir.c_str(),filename.c_str());

    cout<<" We are opening "<<open_filename<<endl;

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
    TNtuple *pulse = new TNtuple("pulse","pulse","pulseHeight:pulseRightEdge:pulseLeftEdge:pulseCharge:pulsePeakTime:CalibratedTime:windowratio");
    TNtuple *event = new TNtuple("event","event","charge:charge_frac:baseline:rms");
	TTree *wforms_tree = new TTree("waveforms","Waveform Tree");
	float waveforms[8192];
	wforms_tree->Branch("pmt_waveforms",waveforms[0],TString::Format("waveforms[%i]/F",number_samples));
    // Store the waveform plot for debugging
    TCanvas *waveplot[100];
    vector<double> baseline_sweep;
    signal_start = true;
    baselinev.clear();
    //reading data
    fin.open(open_filename,ios::in);
    if (!fin.is_open()){
        cout<<" Can't open the file "<<endl;
    }
    else{
        while(!fin.rdstate() ){
        //while(sweep < 682){
            fin.getline(linea,200);//>>datum;
            //cout<<" each line is : "<<linea;
    		*buff=linea;


            if (buff->Contains("Event")){

                if (sweep>0){
                    // For counting dark rate
                    dark_hits->Fill(number_of_peaks);
                    number_of_peaks = 0.0;
					std::copy(raw_waveform.begin(),raw_waveform.end(),waveforms);
					wforms_tree->Fill();
                    rms_value = baseline_rms(baselinev,raw_waveform,number_of_samples);// Calculate baseline and rms then pass to pulse finder
		            baseline_sweep.push_back(rms_value);// save baseline for checking baseline shifting
                    //cout<<"This is sweep : "<<sweep<<" baseline is : "<<rms_value<<endl;
                    //getchar();
		            baselinev.clear();
                    // store first 100 waveform for debugging pulse finder
                    if (sweep < 100){
                        for (int j=0;j<pulse_left_edge.size();j++){
                            t22->SetPoint(j,pulse_left_edge[j],startv[j]);
                            t33->SetPoint(j,pulse_right_edge[j],endv[j]);
                        }

                        char plotname [30];
                        sprintf(plotname,"waveform%d",iana);
                        waveplot[iana] = new TCanvas(plotname);
                        TLine* line = new TLine(0,rms_value,number_of_samples,rms_value);
                        line->SetLineColor(3);
                        line->SetLineStyle(3);
                        line->SetLineWidth(3);
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
                        waveplot[iana]->Write();
                        cout<<" plotting the waveform, this is sweep : "<<sweep<<endl;
                    }
		            if (sweep%1000==0){
                    	cout<<" This is sweep : "<<sweep<<endl;
                    }
		            iana++;
                    pulse_left_edge.clear();
                    pulse_right_edge.clear();
                    startv.clear();
                    endv.clear();
                }

                raw_waveform.clear();
                if (sweep < 100){
                    //cout<<" Create new TGraph ! "<<endl;
                    t11 = new TGraph();
                    t22 = new TGraph();
                    t33 = new TGraph();
                }
                pcount = 0;
                sweep ++;
            }
            //cout<<"pcount is : "<<pcount<<endl;



            else{
                double datum = buff->Atof();
                if (pcount<number_of_samples+1){
                    datum*=-1;

                    datum /= adc_per_Volt;
                    //cout<<" This is sample "<<pcount<<" with value : "<<datum<<endl;
                    //if (isinf(datum) || isnan(datum) || datum > 200000 || datum<-1){
                    //    datum=0;
                    //}
                    if (pcount<100){
                        //if (debug_mode)
                        //    cout<<" Raw data is for baseline : "<<datum<<endl;
                        baselinev.push_back(datum);
                    }

                    raw_waveform.push_back(datum);
                    t11->SetPoint(pcount,pcount,datum);
                    temp_sum = h->GetBinContent(pcount+1);
                    h->SetBinContent(pcount+1,temp_sum+datum);

                }
                pcount++;

            }

        }
    }

    //Fill Ntuple
    //pulseHeight:pulseRightEdge:pulseLeftEdge:pulseCharge:pulsePeakTime
    for (int i=0;i<amplitude.size();i++){
        pulse->Fill(amplitude[i],pr[i],pl[i],charge_v[i],amplitude_position[i],CalibratedTime[i],windowratio[i]);
    }
    TGraph* baseline_plot = new TGraph();
    for (int i=0;i<baseline_sweep.size();i++){
        baseline_plot->SetPoint(i,i,baseline_sweep[i]);
    }
    for (int i=0;i<event_charge.size();i++){
            event->Fill(event_charge[i],event_charge_ten[i],event_baseline[i],event_rms[i]);
    }
    //Baseline plot
    TCanvas* bplot = new TCanvas("bplot","bplot");
    baseline_plot->SetMarkerStyle(22);
    baseline_plot->Draw("AP");

    cout<<" Total sweeps is : "<<sweep<<endl;

	wforms_tree->Write();
    pulse->Write();
    bplot->Write();
    fout->Write();
    fout->Close();

    return 0;
}
