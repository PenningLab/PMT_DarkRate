/////////////////////////////////////////////////////////////////////////////////////////
// To compile : g++ ana_combine_root.cc -o ana_combine_root -I${ROOTSYS}/include/ `root-config --cflags --libs` `python3.6m-config --cflags --ldflags`
// Requires an install of python 3.6, point to your location for python3.6m-config
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
	<<" -t : Use temperature measurements\n"
	<<" -s : Number of sweeps\n"
	<<" -debug : Get in the debugging mode.\n"
	<<" -h or --help : Show the usage\n"
	<<" Enjoy ! -Ryan Wang"<<endl;
}
int main(int argc, char *argv[]){
	int number_files=0;
	string outfilename;
	string infiledir;
	if (argc<2){
		show_usage(argv[0]);
		return 1;
	}
	cout<<"loading in arguments"<<endl;
	for (int i=1;i<argc;++i){
		string arg = argv[i];
		if ((arg=="-h") || (arg=="--help")){
			show_usage(argv[0]);
			return 0;
		}
		else if (arg=="-n"){
			number_files = atoi(argv[i+1]);
		}
		else if (arg=="-o"){
			outfilename = argv[i+1];
		}
		else if (arg=="-i"){
			infiledir = argv[i+1];
		}
	}
	string measurement[2] = {"BNL_test_50ns_181_5cm_2_28_2019_1000_samples_10000_events","BNL_test_50ns_144_7cm_2_28_2019_1000_samples_10000_events"};
	for (int i=0;i<2;i++){

		for (int j=0;j<number_files;j++){
			char filename[200];
			sprintf(filename,"%s/%s/%u_PMT_Trigger.root",infiledir.c_str(),measurement[i].c_str(),j);
			TFile* f = new TFile(filename,"READ");
			TTree* event_tree;
			f->GetObject("event",event_tree);

			vector<double> *charge = 0;
			vector<double> *charge_frac = 0;
			vector<double> *QPE = 0;
			vector<double> *Height = 0;
			vector<double> *start = 0;
			vector<double> *end = 0;
			vector<double> *PeakTime = 0;

			TBranch *bcharge = 0;
			TBranch *bcharge_frac = 0;
			TBranch *bQPE = 0;
			TBranch *bHeight = 0;
			TBranch *bstart = 0;
			TBranch *bend = 0;
			TBranch *bPeakTime = 0;

			event_tree->SetBranchAddress("charge",&charge,&bcharge);
			event_tree->SetBranchAddress("charge_frac",&charge_frac,&bcharge_frac);
			event_tree->SetBranchAddress("QPE",&QPE,&bQPE);
            event_tree->SetBranchAddress("Height",&Height,&bHeight);
            event_tree->SetBranchAddress("stime",&start,&bstart);
            event_tree->SetBranchAddress("etime",&end,&bend);
            event_tree->SetBranchAddress("ptime",&PeakTime,&bPeakTime);

			for (int ie=0;ie<event_tree->GetEntries();ie++){
				Long64_t tentry = event_tree->LoadTree(ie);
      			bcharge->GetEntry(tentry);
				bcharge_frac->GetEntry(tentry);
				bQPE->GetEntry(tentry);
				bHeight->GetEntry(tentry);
				bstart->GetEntry(tentry);
				bend->GetEntry(tentry);
				bPeakTime->GetEntry(tentry);

				cout<<" This is file : "<<filename<<" event : "<<ie<<" has pulses : "<<QPE->size()<<endl;
			}
			getchar();
			//f.Close();
		}



	}




    return 0;
}