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
	TFile* fout = new TFile(outfilename.c_str(),"RECREATE");
	//string measurement[6] = {"BNL_test_50ns_0cm_2_28_2019_1000_samples_10000_events","BNL_test_50ns_38_0cm_2_28_2019_1000_samples_10000_events","BNL_test_50ns_73_0cm_2_28_2019_1000_samples_10000_events","BNL_test_50ns_108_5cm_2_28_2019_1000_samples_10000_events","BNL_test_50ns_144_7cm_2_28_2019_1000_samples_10000_events","BNL_test_50ns_181_5cm_2_28_2019_1000_samples_10000_events"};
	//string measurement[6] = {"BNL_test_50ns_0cm_2_28_2019_1000_samples_10000_events","BNL_test_50ns_NewWater_35_0cm_3_1_2019_1000_samples_10000_events","BNL_test_50ns_NewWater_74_7cm_3_1_2019_1000_samples_10000_events","BNL_test_50ns_NewWater_102_4cm_3_1_2019_1000_samples_10000_events","BNL_test_50ns_NewWater_147_1cm_3_1_2019_1000_samples_10000_events","BNL_test_50ns_NewWater_180_8cm_3_1_2019_1000_samples_10000_events"};
	string measurement[6] = {"BNL_test_50ns_0cm_2_28_2019_1000_samples_10000_events","BNL_new_LED435nm_full_water_stability_500_samples_10000_events","BNL_test_50ns_NewWater_74_7cm_3_1_2019_1000_samples_10000_events","BNL_test_50ns_NewWater_102_4cm_3_1_2019_1000_samples_10000_events","BNL_test_50ns_NewWater_147_1cm_3_1_2019_1000_samples_10000_events","BNL_test_50ns_NewWater_180_8cm_3_1_2019_1000_samples_10000_events"};
	int run_number[6] = {60,1500,60,60,60,60};
	vector<double> charge_ratio;
	vector<double> height_ratio;
	vector<double> charge_ratio_std;
	vector<double> height_ratio_std;

	vector<TH1F*> ratio_q;
	vector<TH1F*> ratio_h;
	vector<TH2F*> charge_dis;
	vector<TH2F*> height_dis;

	TCanvas* ratio_time[number_files];
	TCanvas* top_charge_time[number_files];
	TCanvas* bottom_charge_time[number_files];
	TCanvas* top_height_time[number_files];
	TCanvas* bottom_height_time[number_files];

	TGraphErrors* RPlot[number_files];
	TGraphErrors* TOPlot_charge[number_files];
	TGraphErrors* BOTTOMlot_charge[number_files];
	TGraphErrors* TOPlot_height[number_files];
	TGraphErrors* BOTTOMlot_height[number_files];



	for (int i=0;i<number_files;i++){
		char runname_q[100] ;
		char runname_h[100] ;
		sprintf(runname_q,"%s_chargeRatio",measurement[i].c_str());
		sprintf(runname_h,"%s_heightRatio",measurement[i].c_str());
		char runname_q_2d[100] ;
		char runname_h_2d[100] ;
		sprintf(runname_q_2d,"%s_chargeRatio_2d",measurement[i].c_str());
		sprintf(runname_h_2d,"%s_heightRatio_2d",measurement[i].c_str());
		//cout<<" This is corruption ? "<<endl;
		TH1F* h = new TH1F(runname_q,"",100,0,10);
		TH1F* h2 = new TH1F(runname_h,"",100,0,10);

		TH2F* h3 = new TH2F(runname_q_2d,"",2000,0,500,2000,0,500);
		TH2F* h4 = new TH2F(runname_h_2d,"",1000,0,1,1000,0,1);

		ratio_q.push_back(h);
		ratio_h.push_back(h2);
		charge_dis.push_back(h3);
		height_dis.push_back(h4);

	}

	//double depth[6] = {0,38,73,108.5,144.7,181.5};
	double depth[6] = {0,35,74.7,102.4,147.1,180.8};
	for (int i=0;i<number_files;i++){
		double total_ratio = 0;
		double total_ratio_std = 0;
		double total_ratio_height = 0;
		double total_ratio_height_std = 0;
		double total_ratio_counter = 0;
		//double plotcounter = 0;
		RPlot[i] = new TGraphErrors;
		TOPlot_charge[i] = new TGraphErrors;
		BOTTOMlot_charge[i] = new TGraphErrors;
		TOPlot_height[i] = new TGraphErrors;
		BOTTOMlot_height[i] = new TGraphErrors;
		//bool total_ratio_fire = false;
		for (int j=0;j<run_number[i];j++){
			//if (i==1 && j>9)
			//	break;
			char filename[200];
			sprintf(filename,"%s/%s/%u_PMT_Trigger.root",infiledir.c_str(),measurement[i].c_str(),j);
			TFile* f = new TFile(filename,"READ");
			TTree* event_tree;
			f->GetObject("event",event_tree);


			//ratio_q[i] = new TH1F(runname_q,"",100,0,10);
			//ratio_h[i] = new TH1F(runname_h,"",100,0,10);

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

			double ratio = 0;
			double ratio_std = 0;
			double ratio_height = 0;
			double ratio_height_std = 0;
			double ratio_counter = 0;
			double top_charge_avg = 0;
			double bottom_charge_avg = 0;
			double top_height_avg = 0;
			double bottom_height_avg = 0;
			double top_charge_std = 0;
			double bottom_charge_std = 0;
			double top_height_std = 0;
			double bottom_height_std = 0;

			for (int ie=0;ie<event_tree->GetEntries();ie++){
				Long64_t tentry = event_tree->LoadTree(ie);
      			bcharge->GetEntry(tentry);
				bcharge_frac->GetEntry(tentry);
				bQPE->GetEntry(tentry);
				bHeight->GetEntry(tentry);
				bstart->GetEntry(tentry);
				bend->GetEntry(tentry);
				bPeakTime->GetEntry(tentry);

				//cout<<" This is file : "<<filename<<" event : "<<ie<<" has pulses : "<<QPE->size()<<endl;
				double top_charge = 0;
				double bottom_charge = 0;
				double temp_ratio = 0;
				double top_height = 0;
				double bottom_height = 0;
				double temp_ratio_height = 0;
				bool top_fire = false;
				bool bottom_fire = false;
				//cout<<" Ratio : "<<ratio<<endl;
				for (unsigned int ipulse=0;ipulse<Height->size();ipulse++){
					if (PeakTime->at(ipulse)>195&&PeakTime->at(ipulse)<203){
						top_charge = QPE->at(ipulse);
						top_height = Height->at(ipulse);
						//cout<<" Ratio pulse width : "<<Height->at(ipulse)/(end->at(ipulse)-start->at(ipulse))<<endl;
						//top_fire = true;
						if ((end->at(ipulse)-start->at(ipulse))<15 && (end->at(ipulse)-start->at(ipulse))>10){
							top_fire = true;
						}
						else{
							top_fire = false;
						}
						//cout<<" top charge : "<<top_charge<<endl;
					}
					if (PeakTime->at(ipulse)>210&&PeakTime->at(ipulse)<222){
						bottom_charge = QPE->at(ipulse);
						bottom_height = Height->at(ipulse);
						//bottom_fire = true;
						if ((end->at(ipulse)-start->at(ipulse))<15 && (end->at(ipulse)-start->at(ipulse))>10){
							bottom_fire = true;
						}
						else{
							bottom_fire = false;
						}
						//cout<<" bottom charge : "<<bottom_charge<<endl;
						break;
					}
				}
				if (top_fire && bottom_fire){
					temp_ratio = top_charge/bottom_charge;
					temp_ratio_height = top_height/bottom_height;
					ratio += temp_ratio;
					ratio_std += temp_ratio*temp_ratio;

					charge_dis[i]->Fill(top_charge,bottom_charge);

					ratio_height += temp_ratio_height;
					ratio_height_std += temp_ratio_height*temp_ratio_height;
					ratio_counter ++;
					//cout<<" temp_ratio : "<<temp_ratio<<endl;
					height_dis[i]->Fill(top_height,bottom_height);
					ratio_q[i]->Fill(temp_ratio);
					ratio_h[i]->Fill(temp_ratio_height);

					top_charge_avg += top_charge;
					bottom_charge_avg += bottom_charge;
					top_height_avg += top_height;
					bottom_height_avg += bottom_height;

					top_charge_std += top_charge*top_charge;
					bottom_charge_std += bottom_charge*bottom_charge;
					top_height_std += top_height*top_height;
					bottom_height_std += bottom_height*bottom_height;
				}


			}
			//cout<<" ratio_counter : "<<ratio_counter<<endl;
			ratio /= ratio_counter;
			ratio_std -= ratio*ratio*ratio_counter;
			ratio_std = sqrt(ratio_std/(ratio_counter-1));

			ratio_height /= ratio_counter;
			ratio_height_std -= ratio_height*ratio_height*ratio_counter;
			ratio_height_std = sqrt(ratio_height_std/(ratio_counter-1));

			top_charge_avg /= ratio_counter;
			bottom_charge_avg /= ratio_counter;
			top_height_avg /= ratio_counter;
			bottom_height_avg /= ratio_counter;

			top_charge_std -= top_charge_avg*top_charge_avg*ratio_counter;
			top_charge_std = sqrt(top_charge_std/(ratio_counter-1));

			bottom_charge_std -= bottom_charge_avg*bottom_charge_avg*ratio_counter;
			bottom_charge_std = sqrt(bottom_charge_std/(ratio_counter-1));

			top_height_std -= top_height_avg*top_height_avg*ratio_counter;
			top_height_std = sqrt(top_height_std/(ratio_counter-1));

			bottom_height_std -= bottom_height_avg*bottom_height_avg*ratio_counter;
			bottom_height_std = sqrt(bottom_height_std/(ratio_counter-1));

			TOPlot_charge[i]->SetPoint(j,j,top_charge_avg);
			TOPlot_charge[i]->SetPointError(j,0,top_charge_std);

			BOTTOMlot_charge[i]->SetPoint(j,j,top_charge_avg);
			BOTTOMlot_charge[i]->SetPointError(j,0,bottom_charge_std);

			TOPlot_height[i]->SetPoint(j,j,top_height_avg);
			TOPlot_height[i]->SetPointError(j,0,top_height_std);

			BOTTOMlot_height[i]->SetPoint(j,j,bottom_height_avg);
			BOTTOMlot_height[i]->SetPointError(j,0,bottom_height_std);
			//cout<<"File : "<<measurement[i]<<" This is file : "<<j<<" it has ratio : "<<ratio<<" with std : "<<ratio_std<<endl;
			//cout<<"File : "<<measurement[i]<<" This is file : "<<j<<" it has height ratio : "<<ratio_height<<" with std : "<<ratio_height_std<<endl;
			total_ratio += ratio;
			total_ratio_std += ratio*ratio;
			total_ratio_counter++;

			total_ratio_height += ratio_height;
			total_ratio_height_std += ratio_height*ratio_height;
			//total_ratio_counter++;

			RPlot[i]->SetPoint(j,j,ratio);
			RPlot[i]->SetPointError(j,0,ratio_std);
			//getchar();
			//plotcounter ++;
			f->Close();
		}

		total_ratio /= total_ratio_counter;
		total_ratio_std -= total_ratio*total_ratio*total_ratio_counter;
		total_ratio_std = sqrt(total_ratio_std/(total_ratio_counter-1));

		total_ratio_height /= total_ratio_counter;
		total_ratio_height_std -= total_ratio_height*total_ratio_height*total_ratio_counter;
		total_ratio_height_std = sqrt(total_ratio_height_std/(total_ratio_counter-1));

		cout<<"File : "<<measurement[i]<<" it has ratio : "<<total_ratio<<" with std : "<<total_ratio_std<<endl;
		cout<<"File : "<<measurement[i]<<" it has height ratio : "<<total_ratio_height<<" with height std : "<<total_ratio_height_std<<endl;

		charge_ratio.push_back(total_ratio);
		charge_ratio_std.push_back(total_ratio_std);

		height_ratio.push_back(total_ratio_height);
		height_ratio_std.push_back(total_ratio_height_std);





		//getchar();


	}

	fout->cd();

	TGraphErrors *charge_pl = new TGraphErrors();
	TGraphErrors *height_pl = new TGraphErrors();

	for (int i=0;i<6;i++){
		charge_pl->SetPoint(i,depth[i],charge_ratio[i]);
		charge_pl->SetPointError(i,0,charge_ratio_std[i]);

		height_pl->SetPoint(i,depth[i],height_ratio[i]);
		height_pl->SetPointError(i,0,height_ratio_std[i]);

		char plotname[30];
		sprintf(plotname,"Ratio_plot_%u",i);
		ratio_time[i] = new TCanvas(plotname);
		RPlot[i]->SetMarkerStyle(3);
		RPlot[i]->SetMarkerSize(3);
		RPlot[i]->Draw("AP");
		ratio_time[i]->Write();

		char plotname2[30];
		sprintf(plotname2,"Top_charge_plot_%u",i);
		top_charge_time[i] = new TCanvas(plotname2);
		TOPlot_charge[i]->SetMarkerStyle(3);
		TOPlot_charge[i]->SetMarkerSize(3);
		TOPlot_charge[i]->Draw("AP");
		top_charge_time[i]->Write();

		char plotname3[30];
		sprintf(plotname3,"Bottom_charge_plot_%u",i);
		bottom_charge_time[i] = new TCanvas(plotname3);
		BOTTOMlot_charge[i]->SetMarkerStyle(3);
		BOTTOMlot_charge[i]->SetMarkerSize(3);
		BOTTOMlot_charge[i]->Draw("AP");
		bottom_charge_time[i]->Write();

		char plotname4[30];
		sprintf(plotname4,"Top_height_plot_%u",i);
		top_height_time[i] = new TCanvas(plotname4);
		TOPlot_height[i]->SetMarkerStyle(3);
		TOPlot_height[i]->SetMarkerSize(3);
		TOPlot_height[i]->Draw("AP");
		top_height_time[i]->Write();

		char plotname5[30];
		sprintf(plotname5,"Bottom_height_plot_%u",i);
		bottom_height_time[i] = new TCanvas(plotname5);
		BOTTOMlot_height[i]->SetMarkerStyle(3);
		BOTTOMlot_height[i]->SetMarkerSize(3);
		BOTTOMlot_height[i]->Draw("AP");
		bottom_height_time[i]->Write();
	}


	TCanvas* charge_c = new TCanvas("charge_c","charge_c");
	charge_pl->SetMarkerStyle(22);
	charge_pl->SetMarkerColor(2);
	charge_pl->Draw("AP");

	TCanvas* height_c = new TCanvas("height_c","height_c");
	height_pl->SetMarkerStyle(23);
	height_pl->SetMarkerColor(4);
	height_pl->Draw("AP");

	TCanvas* charge_height_c = new TCanvas("charge_height_c","charge_height_c");
	charge_pl->Draw("AP");
	height_pl->Draw("P");

	charge_c->Write();
	height_c->Write();
	charge_height_c->Write();
	fout->Write();
	fout->Close();



    return 0;
}
