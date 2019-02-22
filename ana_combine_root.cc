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

#include <Python.h>

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
	std::string filename;
	std::string outfilename;
	std::string working_dir;
	std::string out_dir = "";
	int number_of_files = 0;
	bool event_tree_enable = false;
	bool use_temp = false;
	bool trig_pmt = false;
	bool use_frac = false;
	bool calcrate = false;
	int frac_time = 0;
	int frac_start = 0;
	float charge_threshold = -1;
	double num_sweeps=-1;
	std::string pydir;
	int initial_run = 0;
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
		else if (arg=="-s"){
			num_sweeps = atof(argv[i+1]);
		}
		else if (arg=="-t"){
			use_temp = true;
		}
		else if (arg=="-init"){
			initial_run = atoi(argv[i+1]);
		}
		else if (arg=="-e"){
			event_tree_enable = true;
		}
		else if (arg=="-pmt"){
		  trig_pmt = true;
        }
		else if (arg=="-frac"){
			use_frac = true;
			frac_time = atof(argv[i+1]);
			frac_start = atof(argv[i+2]);
		}
		else if (arg=="-rate"){
			calcrate = true;
			pydir = argv[i+1];
			//pydir.append("/readlogs.py");
		}

	}
	if (out_dir=="") out_dir = working_dir;
	std::ifstream temp_file;
	if (use_temp){
		char in_temp[320];
		sprintf(in_temp,"%s/%s",working_dir.c_str(),"temptimes.txt");
		temp_file.open(in_temp,std::ifstream::in);
	}

	double myrate = 0;
	if(calcrate){
		std::cout<<"Initializing python"<<std::endl;
		Py_Initialize();
		std::cout<<"Updating syspath with pydir: "<<pydir<<std::endl;
		PyObject* sysPath = PySys_GetObject("path");
		PyList_Append(sysPath, PyUnicode_FromFormat("%s",pydir.c_str()));

		std::cout<<"Loading module"<<std::endl;

		//Load module
		PyObject *pName = PyUnicode_FromString("readlogs");
		PyObject *pModule = PyImport_Import(pName);
		Py_DECREF(pName);
		if (pModule != NULL){
			std::cout << "Py Module Found" << std::endl;

			// Get function from module
			PyObject *pFunc = PyObject_GetAttrString(pModule,"calcrates");
			if(pFunc && PyCallable_Check(pFunc)){
				PyObject* pargs = PyTuple_New(2);
				PyObject* pval;
				pval =  PyUnicode_FromFormat("%s",working_dir.c_str());
				PyTuple_SetItem(pargs,0,pval);
				pval = PyLong_FromLong(number_of_files);
				PyTuple_SetItem(pargs,1,pval);

				PyObject* myresult = PyObject_CallObject(pFunc,pargs);
				Py_DECREF(pargs);
				Py_DECREF(pval);
				if(myresult != NULL){
					myrate = PyFloat_AsDouble(myresult);
					std::cout<<"Rate is "<<myrate<<std::endl;
					Py_DECREF(myresult);
				}
				else{
					Py_DECREF(pFunc);
					Py_DECREF(pModule);
					PyErr_Print();
					std::cout<<"Failed to get result"<<std::endl;
				}

			}
			else{
				if (PyErr_Occurred())
                	PyErr_Print();
				std::cout<<"Couldn't find calcrates"<<std::endl;
			}
			Py_XDECREF(pFunc);
        	Py_DECREF(pModule);
		}
		else{
			PyErr_Print();
			std::cout<<"Failed to load module"<<std::endl;
		}
		Py_FinalizeEx();
	}

	std::cout<<"\nCreating output file"<<std::endl;
	char out_path [320];
	sprintf(out_path,"%s/%s",out_dir.c_str(),outfilename.c_str());

	TFile *fout = new TFile(out_path,"RECREATE");
	TNtuple *pulse = new TNtuple("pulse","pulse","pulseHeight:pulseRightEdge:pulseLeftEdge:pulseCharge:pulsePeakTime:CalibratedTime:baselinerms:windowratio:Run:sweep:bigstep:pulseLength5:pulseLength25:pulseLength50:pulseLength75:pulseLength80:pulseLength90:pulseLength95:pulseLength99:triggerpulseHeight:triggerpulseWidth:triggerpulsePeakTime");
	//TNtuple *event = new TNtuple("event","event","charge:charge_frac:baseline:rms:npulses");

	//}
	std::cout<<"creating output ntuple"<<std::endl;
	// variables
	float pulseHeight=0,pulseRightEdge=0,pulseLeftEdge=0,pulseCharge=0,pulsePeakTime=0,CalibratedTime=0,windowRatio=0,baselinerms=0,sweep=0,bigstep=0,pl5=0,pl25=0,pl50=0,pl75=0,pl80=0,pl90=0,pl95=0,pl99=0;
	float charge=0,charge_frac=0,baseline=0,rms=0,npeaks=0,firstTime=0,triggerpulseHeight=0,triggerpulseWidth=0,triggerpulsePeakTime=0,mycharge_fracj=0;
	TTree *event = new TTree("event","event");
	event->Branch("charge",&charge,"charge/F");
	event->Branch("charge_frac",&charge_frac,"charge_frac/F");
	event->Branch("baseline",&baseline,"baseline/F");
	event->Branch("rms",&rms,"rms/F");
	event->Branch("npulses",&npeaks,"npulse/F");
	if(use_frac) event->Branch("mycharge_frac",&mycharge_fracj,"mycharge_fracj/F");
	//    short int sweep=0;
	//int run=0;
	vector<double> dark_count;
	vector<double> dark_count_error;

	vector<double> rtd1;
	vector<double> rtd2;
	vector<double> rtd3;
	vector<double> rtd4;

	cout<<"start looping"<<endl;
	for (int i=initial_run;i<number_of_files;i++){
		char root_file_name [320];
		sprintf(root_file_name,"%s/%u_%s",working_dir.c_str(),i,filename.c_str());
		cout<<"Reading in file"<<endl;
		TTree* tree;
		TTree* event_tree;
		TFile *fin = new TFile(root_file_name,"READ");
		if (fin == NULL || fin->IsZombie()){
			cout<<" File is corrupted ! "<<endl;
			dark_count.push_back(-1);
			dark_count_error.push_back(-1);
			continue;
		}
		else{
			tree = (TTree*) fin->Get("pulse");
			if (event_tree_enable){
				event_tree = (TTree*) fin->Get("event");
				cout<<" event tree enabled "<<endl;
			}
		}
		double nos = fin->Get("Nsamples")->GetUniqueID();
		TH1F* dark_hit = (TH1F*) fin->Get("dark_hits");
		dark_count.push_back(dark_hit->GetMean()/nos);
		dark_count_error.push_back(dark_hit->GetMeanError()/nos);


		if(use_temp){
			double n,irtd1,irtd2,irtd3,irtd4;
			temp_file >> n >> irtd1 >> irtd2 >> irtd3 >> irtd4;
			rtd1.push_back(irtd1);
			rtd2.push_back(irtd2);
			rtd3.push_back(irtd3);
			rtd4.push_back(irtd4);
		}
		cout<<"Loading tree branches"<<endl;
		tree->SetBranchAddress("pulseHeight",&pulseHeight);
		tree->SetBranchAddress("pulseRightEdge",&pulseRightEdge);
		tree->SetBranchAddress("pulseLeftEdge",&pulseLeftEdge);
		tree->SetBranchAddress("pulseCharge",&pulseCharge);
		tree->SetBranchAddress("pulsePeakTime",&pulsePeakTime);
		tree->SetBranchAddress("CalibratedTime",&CalibratedTime);
		tree->SetBranchAddress("windowratio",&windowRatio);
		tree->SetBranchAddress("baselinerms",&baselinerms);
		tree->SetBranchAddress("sweep",&sweep);
		tree->SetBranchAddress("bigstep",&bigstep);
		tree->SetBranchAddress("pulseLength5",&pl5);
		tree->SetBranchAddress("pulseLength25",&pl25);
		tree->SetBranchAddress("pulseLength50",&pl50);
		tree->SetBranchAddress("pulseLength75",&pl75);
		tree->SetBranchAddress("pulseLength80",&pl80);
		tree->SetBranchAddress("pulseLength90",&pl90);
		tree->SetBranchAddress("pulseLength95",&pl95);
		tree->SetBranchAddress("pulseLength99",&pl99);
		//tree->SetBranchAddress("triggerpulseHeight",&triggerpulseHeight);
		//tree->SetBranchAddress("triggerpulseWidth",&triggerpulseWidth);
		//tree->SetBranchAddress("triggerpulsePeakTime",&triggerpulsePeakTime);
		cout<<" processing root file No. "<<i<<endl;

		int nument = tree->GetEntries();
		std::vector<double> mycharge_frac(1,0);
		if(event_tree_enable){
			int numevts = event_tree->GetEntries();
			if(use_frac) mycharge_frac.resize(numevts,0);
		}
		for (int j=0;j<nument;j++){
			tree->GetEntry(j);
			float inmount[] = {pulseHeight,pulseRightEdge,pulseLeftEdge,pulseCharge,pulsePeakTime,CalibratedTime,baselinerms,windowRatio,i,sweep,bigstep,pl5,pl25,pl50,pl75,pl80,pl90,pl95};
			float* inmountpoint = inmount;
			pulse->Fill(inmountpoint);
			if(use_frac && pulsePeakTime>frac_start && pulsePeakTime<(frac_start+frac_time)) mycharge_frac[sweep]+=pulseCharge;
			//cout<<" This is root file : "<<i<<" we are reading entry : "<<j<<" with pulseHeight : "<<pulseHeight<<endl;
		}
		if (event_tree_enable){
			event_tree->SetBranchAddress("charge",&charge);
			event_tree->SetBranchAddress("charge_frac",&charge_frac);
			event_tree->SetBranchAddress("baseline",&baseline);
			event_tree->SetBranchAddress("rms",&rms);
			event_tree->SetBranchAddress("npulses",&npeaks);

			int numevts = event_tree->GetEntries();
			for (int j =0;j<event_tree->GetEntries();j++){
				event_tree->GetEntry(j);
				if(use_frac) mycharge_fracj = mycharge_frac[j];
				event->Fill();
			}
			delete event_tree;
		}
		delete tree;
		fin->Close();
	}//main for loop
	fout->cd();

	TGraphErrors* dark_plot = new TGraphErrors();
	int backcount =0;
	for (int h=0;h<dark_count.size();h++){
	  double temp_dark_rate = dark_count[h]*1E8;
		double temp_dark_rate_error = dark_count_error[h]*1E8;
		if(temp_dark_rate<0){
		  backcount++;
		  continue;
		}
		dark_plot->SetPoint(h-backcount,h,temp_dark_rate);
		dark_plot->SetPointError(h-backcount,0,temp_dark_rate_error);
	}

	if(use_temp){
		TGraphErrors* prtd1 = new TGraphErrors();
		TGraphErrors* prtd2 = new TGraphErrors();
		TGraphErrors* prtd3 = new TGraphErrors();
		TGraphErrors* prtd4 = new TGraphErrors();
		for(int h=0;h<dark_count.size();h++){
			double temp_dark_rate = dark_count[h]*1E8;
			double temp_dark_rate_error = dark_count_error[h]*1E8;
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

	dark_plot->SetName("dark_plotfromhist");
	dark_plot->SetTitle(";Run (100s);Dark Rate (Hz)");
	TCanvas* cdark = new TCanvas("cdark","cdark");
	dark_plot->SetMarkerStyle(24);
	dark_plot->SetMarkerColor(2);
	dark_plot->Draw("AP");

	TVectorD mrate(1);
	mrate[0] = myrate;
	mrate.Write("Rate");

	dark_plot->Write();
	cdark->Write();
	pulse->Write();
	fout->Write();
	fout->Close();
	return 0;
}
