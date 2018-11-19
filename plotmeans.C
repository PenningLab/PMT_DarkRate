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
using namespace std;

void plotmeans(TString filename){
 TFile *infile = new TFile(filename,"READ");
 TGraphErrors *rtd4 = (TGraphErrors*)infile->Get("RTD4");
 int npoints = rtd4->GetN();
 vector<double> temps;
 vector<double> temperr;
 vector<double> darkrate;
 vector<double> darkerr;

 TFile* outfile = new TFile("outplots.root","RECREATE");
 for(int i=0;i<npoints;i+=5){
  double temprate=0;
  double temprate_err=0;
  double temptemp=0;
  double temptemperr=0;
  temps.push_back(0);
  darkrate.push_back(0);
  darkerr.push_back(0);
  temperr.push_back(0);

  rtd4->GetPoint(i,temptemp,temprate);
  temps.back()+=temptemp;
  darkrate.back()+=temprate;
  temprate_err+=rtd4->GetErrorY(i) * rtd4->GetErrorY(i);
  temptemperr+=rtd4->GetErrorX(i) * rtd4->GetErrorX(i);

  rtd4->GetPoint(i+1,temptemp,temprate);
  temps.back()+=temptemp;
  darkrate.back()+=temprate;
  temprate_err+=rtd4->GetErrorY(i+1) * rtd4->GetErrorY(i+1);
  temptemperr+=rtd4->GetErrorX(i+1) * rtd4->GetErrorX(i+1);

  rtd4->GetPoint(i+2,temptemp,temprate);
  temps.back()+=temptemp;
  darkrate.back()+=temprate;
  temprate_err+=rtd4->GetErrorY(i+2) * rtd4->GetErrorY(i+2);
  temptemperr+=rtd4->GetErrorX(i+2) * rtd4->GetErrorX(i+2);

  rtd4->GetPoint(i+3,temptemp,temprate);
  temps.back()+=temptemp;
  darkrate.back()+=temprate;
  temprate_err+=rtd4->GetErrorY(i+3) * rtd4->GetErrorY(i+3);
  temptemperr+=rtd4->GetErrorX(i+3) * rtd4->GetErrorX(i+3);

  rtd4->GetPoint(i+4,temptemp,temprate);
  temps.back()+=temptemp;
  darkrate.back()+=temprate;
  temprate_err+=rtd4->GetErrorY(i+4) * rtd4->GetErrorY(i+4);
  temptemperr+=rtd4->GetErrorX(i+4) * rtd4->GetErrorX(i+4);

  temps.back()*=0.2;
  darkrate.back()*=0.2;
  darkerr.back()=sqrt(temprate_err);
  temperr.back()=sqrt(temptemperr);
 }

 TGraphErrors* rtd4means = new TGraphErrors();
 for (int h=0;h<temps.size();h++){
  rtd4means->SetPoint(h,temps[h],darkrate[h]);
  rtd4means->SetPointError(h,0,darkerr[h]);
 }
 rtd4means->SetName("rtd4");
 TCanvas* cdark = new TCanvas("cdark","cdark");
 rtd4means->SetMarkerStyle(24);
 rtd4means->SetMarkerColor(2);
 rtd4means->Draw("AP");

 cdark->Draw();
 float pulseCharge=0,Run=0;
 TTree* tree = (TTree*)infile->Get("pulse");
 tree->SetBranchAddress("pulseCharge",&pulseCharge);
 tree->SetBranchAddress("Run",&Run);
 int nument = tree->GetEntries();
 vector<double> multiplierup(temps.size(),0);
 vector<double> multiplierdown(temps.size(),0);
 vector<double> count(temps.size(),0);
 for (int i=0; i<nument;i++){
  tree->GetEntry(i);
  count[(int)Run/5]++;
  if(pulseCharge/2.0 <1.5) multiplierdown[(int)Run/5]++;
  else if(pulseCharge/2.0 >2) multiplierup[(int)Run/5]++;
 }
 TGraphErrors* chargelow,*chargeup;
 chargeup=new TGraphErrors();
 chargelow = new TGraphErrors();
 for(int j=0;j<temps.size();j++){
  double uprate = darkrate[j]*multiplierup[j]/count[j];
  double downrate = darkrate[j]*multiplierdown[j]/count[j];
  chargeup->SetPoint(j,temps[j],uprate);
  chargelow->SetPoint(j,temps[j],downrate);

  double uperr = darkerr[j]*multiplierup[j]/count[j];
  double downerr = darkerr[j]*multiplierdown[j]/count[j];
  chargeup->SetPointError(j,0,uperr);
  chargelow->SetPointError(j,0,downerr);
 }

 chargeup->SetName("highcharge");
 chargelow->SetName("lowcharge");
 rtd4means->Write();
 chargeup->Write();
 chargelow->Write();
 outfile->Write();

}
