#define PulseComp_cxx
// The class definition in PulseComp.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("PulseComp.C")
// root> T->Process("PulseComp.C","some options")
// root> T->Process("PulseComp.C+")
//

#include "PulseComp.h"
#include <TH2.h>
#include <TStyle.h>

#include <cstdio>
#include <cstdlib>
#include <dirent.h>
#include <fstream> // std::ifstream
#include <iostream> // std::cout
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <sys/stat.h> // for check directory
#include <sys/types.h>
#include <time.h>
#include <vector>

void PulseComp::Begin(TTree* /*tree*/)
{
	// The Begin() function is called at the start of the query.
	// When running with PROOF Begin() is only called on the client.
	// The tree argument is deprecated (on PROOF 0 is passed).

	TString option = GetOption();
}

void PulseComp::SlaveBegin(TTree* /*tree*/)
{
	// The SlaveBegin() function is called after the Begin() function.
	// When running with PROOF SlaveBegin() is called on each slave server.
	// The tree argument is deprecated (on PROOF 0 is passed).

	TString option = GetOption();
	h_avgphdlg = new TH1F("h_avgphdlg", "Average LG Pulse;Sample (10ns);V", 200, -100, 100);
	h_avgphdlg->Sumw2();
	fOutput->Add(h_avgphdlg);
	h_avgphdhg = new TH1F("h_avgphdhg", "Average HG Pulse;Sample (10ns);V", 200, -100, 100);
	h_avgphdhg->Sumw2();
	fOutput->Add(h_avgphdhg);
	// fOutput->Add(numpulses);
}

Bool_t PulseComp::Process(Long64_t entry)
{
	// The Process() function is called for each entry in the tree (or possibly
	// keyed object in the case of PROOF) to be processed. The entry argument
	// specifies which entry in the currently loaded tree is to be processed.
	// When processing keyed objects with PROOF, the object is already loaded
	// and is available via the fObject pointer.
	//
	// This function should contain the \"body\" of the analysis. It can contain
	// simple or elaborate selection criteria, run algorithms on the data
	// of the event and typically fill histograms.
	//
	// The processing can be stopped by calling Abort().
	//
	// Use fStatus to set the return value of TTree::Process().
	//
	// The return value is currently not used.

	fReader.SetLocalEntry(entry);
	std::vector<float> rawform;
	if (getwaveform(rawform))
	{
		// std::cout << "blastin off again" << std::endl;
		int initials = pulse_left_edge - 10;
		int finals = pulse_right_edge + 10;
		for (int k = initials; k < finals; k++)
		{
			if (*channel >= 1000)
				h_avgphdlg->Fill(k - pulse_peak_time, rawform[k] - baseline);
			else
				h_avgphdhg->Fill(k - pulse_peak_time, rawform[k] - baseline);
		}
		if (*channel >= 1000)
			mynumpulseslg++;
		else
			mynumpulseshg++;
	}
	return kTRUE;
}

void PulseComp::SlaveTerminate()
{
	// The SlaveTerminate() function is called after all entries or objects
	// have been processed. When running with PROOF SlaveTerminate() is called
	// on each slave server.
	numpulseshg += mynumpulseshg;
	numpulseslg += mynumpulseslg;
}

void PulseComp::Terminate()
{
	// The Terminate() function is the last function to be called during
	// a query. It always runs on the client, it can be used to present
	// the results graphically or save the results to file.
	TFile* outfile = new TFile("MeanODPulse.root", "recreate");
	outfile->cd();

	h_avgphdhg = dynamic_cast<TH1F*>(fOutput->FindObject("h_avgphdhg"));
	h_avgphdhg->Scale(1.0 / (double)numpulseshg);
	h_avgphdhg->Write();

	h_avgphdlg = dynamic_cast<TH1F*>(fOutput->FindObject("h_avgphdlg"));
	h_avgphdlg->Scale(1.0 / (double)numpulseslg);
	h_avgphdlg->Write();
	outfile->Close();
}
