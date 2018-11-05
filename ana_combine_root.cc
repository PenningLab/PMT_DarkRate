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
static void show_usage(string name){
    cout<<" Usage : ./ana_combine_root [-co] file1 "<<name<<" Options:\n"
    <<" -o : Name of output file.\n"
    <<" -i : Name of input file.\n"
    <<" -n : Number of root files\n"
    <<" -wd : Working directory\n"
    <<" -od : Output directory\n"
    <<" -ct : Charge threshold for dark rate\n"
    <<" -debug : Get in the debugging mode.\n"
    <<" -h or --help : Show the usage\n"
    <<" Enjoy ! -Ryan Wang"<<endl;
}
int main(int argc, char *argv[]){
    string filename;
    string outfilename;
    string working_dir;
    string out_dir = "";
    int number_of_files = 0;
    bool event_tree_enable = false;
		bool use_temp = false;
    float charge_threshold = -1;
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
        else if (arg=="-od"){
            out_dir = argv[i+1] ;
        }
        else if (arg=="-n"){
            number_of_files = atoi(argv[i+1]);
        }
        else if (arg=="-i") {
            filename = argv[i+1];
        }
        else if (arg=="-o"){
            outfilename = argv[i+1];
        }
        else if (arg=="-ct"){
            charge_threshold = atof(argv[i+1]);
        }
				else if (arg=="-t"){
            use_temp = true;
        }
	else if (arg=="-e"){
		event_tree_enable = true;
	}
    }
    if (out_dir=="") out_dir = working_dir;
    ifstream temp_file;
		if (use_temp){
			char in_temp[320];
			sprintf(in_temp,"%s/%s",working_dir.c_str(),"temptimes.txt");
			temp_file.open(in_temp,std::ifstream::in);
		}
    char out_path [320];
    sprintf(out_path,"%s/%s",out_dir.c_str(),outfilename.c_str());
    TFile *fout = new TFile(out_path,"RECREATE");
    TNtuple *pulse = new TNtuple("pulse","pulse","pulseHeight:pulseRightEdge:pulseLeftEdge:pulseCharge:pulsePeakTime:CalibratedTime:Run");
    //if (event_tree_enable){
    	TNtuple *event = new TNtuple("event","event","charge:charge_frac:baseline:rms");
    //}
    // variables
    float pulseHeight=0,pulseRightEdge=0,pulseLeftEdge=0,pulseCharge=0,pulsePeakTime=0,CalibratedTime=0;
    float charge=0,charge_frac=0,baseline=0,rms=0;
    //int run=0;
    vector<double> dark_count;
    vector<double> dark_count_error;
    vector<double> dark_countcut;
    vector<double> dark_countcut_error;
		vector<double> rtd1;
		vector<double> rtd2;
		vector<double> rtd3;
		vector<double> rtd4;
    TTree* tree;
    TTree* event_tree;
    cout<<"start looping"<<endl;
    for (int i=0;i<number_of_files;i++){
        char root_file_name [320];
        sprintf(root_file_name,"%s/%u_%s",working_dir.c_str(),i,filename.c_str());
        TFile *fin = new TFile(root_file_name,"READ");
        if (fin == NULL){
            cout<<" File is corrupted ! "<<endl;
            continue;
        }
        else{
            tree = (TTree*) fin->Get("pulse");
	    if (event_tree_enable){
	    	event_tree = (TTree*) fin->Get("event");
		cout<<" event tree enabled "<<endl;
	    }
        }
        TH1F* dark_hit = (TH1F*) fin->Get("dark_hits");
        dark_count.push_back(dark_hit->GetMean());
        dark_count_error.push_back(dark_hit->GetMeanError());

        TH2F* dark_hit2d = (TH2F*) fin->Get("dark_hits2d");
	if(charge_threshold!=-1){
		dark_hit2d->GetYaxis()->SetRange(0,dark_hit2d->GetYaxis()->FindBin(charge_threshold));
		dark_countcut.push_back(dark_hit2d->GetMean(1));
		dark_countcut_error.push_back(dark_hit2d->GetMean(11));
	}
			if(use_temp){
				double n,irtd1,irtd2,irtd3,irtd4;
				temp_file >> n >> irtd1 >> irtd2 >> irtd3 >> irtd4;
				rtd1.push_back(irtd1);
				rtd2.push_back(irtd2);
				rtd3.push_back(irtd3);
				rtd4.push_back(irtd4);
			}

        tree->SetBranchAddress("pulseHeight",&pulseHeight);
        tree->SetBranchAddress("pulseRightEdge",&pulseRightEdge);
        tree->SetBranchAddress("pulseLeftEdge",&pulseLeftEdge);
        tree->SetBranchAddress("pulseCharge",&pulseCharge);
        tree->SetBranchAddress("pulsePeakTime",&pulsePeakTime);
        tree->SetBranchAddress("CalibratedTime",&CalibratedTime);
        cout<<" processing root file No. "<<i<<endl;

	if (event_tree_enable){
        	event_tree->SetBranchAddress("charge",&charge);
		event_tree->SetBranchAddress("charge_frac",&charge_frac);
		event_tree->SetBranchAddress("baseline",&baseline);
		event_tree->SetBranchAddress("rms",&rms);
        }
        int nument = tree->GetEntries();

        for (int j=0;j<nument;j++){
            tree->GetEntry(j);
            pulse->Fill(pulseHeight,pulseRightEdge,pulseLeftEdge,pulseCharge,pulsePeakTime,CalibratedTime,i);
            //cout<<" This is root file : "<<i<<" we are reading entry : "<<j<<" with pulseHeight : "<<pulseHeight<<endl;
        }
	if (event_tree_enable){
		for (int j =0;j<event_tree->GetEntries();j++){
			event_tree->GetEntry(j);
			event->Fill(charge,charge_frac,baseline,rms);
		}
	}
        fin->Close();
    }//main for loop
    fout->cd();

    TGraphErrors* dark_plot = new TGraphErrors();
		
    for (int h=0;h<dark_count.size();h++){
        double temp_dark_rate = dark_count[h]*1E6/80;
        double temp_dark_rate_error = dark_count_error[h]*1E6/80;
        dark_plot->SetPoint(h,h,temp_dark_rate);
        dark_plot->SetPointError(h,0,temp_dark_rate_error);
    }

		if(use_temp){

			TGraphErrors* prtd1 = new TGraphErrors();
			TGraphErrors* prtd2 = new TGraphErrors();
			TGraphErrors* prtd3 = new TGraphErrors();
			TGraphErrors* prtd4 = new TGraphErrors();
				for(int h=0;h<dark_count.size();h++){
					double temp_dark_rate = dark_count[h]*1E6/80;
					double temp_dark_rate_error = dark_count_error[h]*1E6/80;
					prtd1->SetPoint(h,rtd1[h],temp_dark_rate);
					prtd2->SetPoint(h,rtd2[h],temp_dark_rate);
					prtd3->SetPoint(h,rtd3[h],temp_dark_rate);
					prtd4->SetPoint(h,rtd4[h],temp_dark_rate);

					prtd1->SetPointError(h,0,temp_dark_rate_error);
					prtd2->SetPointError(h,0,temp_dark_rate_error);
					prtd3->SetPointError(h,0,temp_dark_rate_error);
					prtd4->SetPointError(h,0,temp_dark_rate_error);
				}
			
			prtd1->Sort();
			prtd2->Sort();
			prtd3->Sort();
			prtd4->Sort();

			prtd1->SetName("RTD1");
			prtd2->SetName("RTD2");
			prtd3->SetName("RTD3");
			prtd4->SetName("RTD4");

			prtd1->SetTitle("RTD1;Temp [C];DarkRate (Hz)");
			prtd2->SetTitle("RTD2;Temp [C];DarkRate (Hz)");
			prtd3->SetTitle("RTD3;Temp [C];DarkRate (Hz)");
			prtd4->SetTitle("RTD4;Temp [C];DarkRate (Hz)");
			
			prtd1->Write(); prtd2->Write(); prtd3->Write(); prtd4->Write();
		}
    dark_plot->SetName("dark_plotraw");
    dark_plot->SetTitle(";Run (100s);Dark Rate (Hz)");
    TCanvas* cdark = new TCanvas("cdark","cdark");
    dark_plot->SetMarkerStyle(24);
    dark_plot->SetMarkerColor(2);
    dark_plot->Draw("AP");



    if(charge_threshold!=-1){
    TGraphErrors* dark_plotcut = new TGraphErrors();
    for (int h=0;h<dark_countcut.size();h++){
        double temp_dark_rate = dark_countcut[h]*1E6/80;
        double temp_dark_rate_error = dark_countcut_error[h]*1E6/80;
        dark_plotcut->SetPoint(h,h,temp_dark_rate);
        dark_plotcut->SetPointError(h,0,temp_dark_rate_error);
    }
    dark_plotcut->SetName("dark_plotcut");
    dark_plotcut->Write();
    }

    dark_plot->Write();
    cdark->Write();
    pulse->Write();
    fout->Write();
    fout->Close();
    return 0;
}
