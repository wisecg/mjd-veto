// auto-veto.cc
// Runs during production, creates ROOT files of MJD veto data.
// C. Wiseman, A. Lopez, 10/24/16.
//
// NOTE: The scans are split across a few different loops over the events in the run.
// This is done to increase the flexibility of the code, since it checks many different
// quantities.  The performance hit should be minimal, since the size of the veto trees
// is relatively small.

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"
#include "TFile.h"
#include "MJVetoEvent.hh"
#include "GATDataSet.hh"

using namespace std;

vector<int> VetoThreshFinder(TChain *vetoChain, string outputDir, bool makePlots=false);
void ProcessVetoData(TChain *vetoChain, vector<int> thresholds, string outputDir, bool errorCheckOnly=false);

void SetCardNumbers(int runNum, int &card1, int &card2);
int FindPanelThreshold(TH1D *qdcHist, int threshVal, int panel, int runNum);
double InterpTime(int entry, vector<double> times, vector<bool> badScaler);
int PlaneMap(int qdcChan, int runNum=0);
bool CheckEventErrors(MJVetoEvent veto, MJVetoEvent prev, MJVetoEvent first, long prevGoodEntry, vector<int> &ErrorVec);
bool CheckEventErrors(MJVetoEvent veto, MJVetoEvent prev, MJVetoEvent first, long prevGoodEntry);

int main(int argc, char** argv)
{
	// get command line args
	if (argc < 2) {
		cout << "Usage: ./auto-veto [run number]\n"
			  << "                   [-d (optional: draws QDC & multiplicity plots)]\n"
			  << "                   [-e (optional: error check only)]\n";
		return 1;
	}
	int run = stoi(argv[1]);
	if (run > 60000000 && run < 70000000) {
		cout << "Veto data not present in Module 2 runs.  Exiting ...\n";
		return 1;
	}
	bool draw = false, errorCheckOnly = false;
	string opt1 = "", opt2 = "";
	if (argc > 2) opt1 = argv[2];
	if (argc > 3) opt2 = argv[3];
	if (argc > 2 && (opt1 == "-d" || opt2 == "-d"))
		draw = true;
	if (argc > 2 && (opt1 == "-e" || opt2 == "-e"))
		errorCheckOnly = true;

	// output file directory
	string outputDir = "./avout";

	// Only get the run path (so we can run with veto-only runs too)
	// GATDataSet ds;
	// string runPath = ds.GetPathToRun(run,GATDataSet::kBuilt);

	string runPath = "./stage/OR_run"+std::to_string(run)+".root";

	TChain *vetoChain = new TChain("VetoTree");

	// Verify that the veto data exists in the given run
	if (!vetoChain->Add(runPath.c_str())){
		cout << "File doesn't exist.  Exiting ...\n";
		return 1;
	}

	printf("\n========= Processing run %i ... %lli entries. =========\n",run,vetoChain->GetEntries());
	cout << "Path: " << runPath << endl;

	// Find the QDC pedestal location in each channel.
	// Set a software threshold value above this location,
	// and optionally output plots that confirm this choice.
	vector<int> thresholds = VetoThreshFinder(vetoChain, outputDir, draw);

	// Check for data quality errors,
	// tag muon and LED events in veto data,
	// and output a ROOT file for further analysis.
	ProcessVetoData(vetoChain,thresholds,outputDir,errorCheckOnly);

	printf("=================== Done processing. ====================\n\n");
	return 0;
}

vector<int> VetoThreshFinder(TChain *vetoChain, string outputDir, bool makePlots)
{
	// format: (panel 1) (threshold 1) (panel 2) (threshold 2) ...
	vector<int> thresholds;
	int threshVal = 35;	// how many QDC above the pedestal we set the threshold at

	long vEntries = vetoChain->GetEntries();
	TTreeReader reader(vetoChain);
	TTreeReaderValue<unsigned int> vMult(reader, "mVeto");
	TTreeReaderValue<uint32_t> vBits(reader, "vetoBits");
	TTreeReaderValue<MGTBasicEvent> vEvt(reader,"vetoEvent");
	TTreeReaderValue<MJTRun> vRun(reader,"run");
	reader.Next();
	int runNum = vRun->GetRunNumber();
	reader.SetTree(vetoChain);  // resets the reader

	gStyle->SetOptStat(0);
	int bins=500, lower=0, upper=500;
	TH1D *hLowQDC[32];
	TH1D *hFullQDC[32];
	char hname[50];
	for (int i = 0; i < 32; i++) {
		sprintf(hname,"hLowQDC%d",i);
		hLowQDC[i] = new TH1D(hname,hname,bins,lower,upper);
		sprintf(hname,"hFullQDC%d",i);
		hFullQDC[i] = new TH1D(hname,hname,420,0,4200);
	}
	sprintf(hname,"Run %i Hit Multiplicity",runNum);
	TH1D *hMultip = new TH1D("hMultip",hname,32,0,32);

	// Set all thresholds to 1, causing all entries to have a multiplicity of 32
	int def[32] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
	long skippedEvents = 0, prevGoodEntry=0;

	// MJVetoEvent variables, with run-based card numbers
	int card1=0, card2=0;
	SetCardNumbers(runNum,card1,card2);
	MJVetoEvent veto(card1,card2);
	MJVetoEvent prev, first;
	cout << "QDC 1 in slot " << card1 << ", QDC 2 in slot " << card2 << endl;

	while (reader.Next())
	{
		long i = reader.GetCurrentEntry();
		int run = vRun->GetRunNumber();

		veto.SetSWThresh(def);
    	veto.WriteEvent(i,&*vRun,&*vEvt,*vBits,run,true);
		if (CheckEventErrors(veto,prev,first,prevGoodEntry))
		{
    		skippedEvents++;
			// do the end of event resets before continuing
			prev = veto;
			prevGoodEntry = i;
			veto.Clear();
    		continue;
    	}
    	for (int q = 0; q < 32; q++) {
    		hLowQDC[q]->Fill(veto.GetQDC(q));
    		hFullQDC[q]->Fill(veto.GetQDC(q));
		}
		// save previous entries for the event error check
		prev = veto;
		prevGoodEntry = i;
		veto.Clear();
	}
	if (skippedEvents > 0) printf("VetoThreshFinder skipped %li of %li entries.\n",skippedEvents,vEntries);

	int thresh[32] = {9999};
	for (int i = 0; i < 32; i++)
	{
		thresh[i] = FindPanelThreshold(hLowQDC[i],threshVal,i,runNum);
		thresholds.push_back(i);
		thresholds.push_back(thresh[i]);
	}
	// cout << "vtf thresholds: " << endl;
	// for (int i = 0; i < 32; i++)
		// cout << i << " " << thresh[i] << endl;

	// re-scan with the found thresholds to make a multiplicity plot
	reader.SetTree(vetoChain);  // resets the reader
	while (reader.Next())
	{
		long i = reader.GetCurrentEntry();
		int run = vRun->GetRunNumber();
		veto.SetSWThresh(thresh);
    	veto.WriteEvent(i,&*vRun,&*vEvt,*vBits,run,true);
		if (CheckEventErrors(veto,prev,first,prevGoodEntry))
		{
    		skippedEvents++;
			// do the end of event resets before continuing
			prev = veto;
			prevGoodEntry = i;
			veto.Clear();
    		continue;
    	}
		hMultip->Fill(veto.GetMultip());
		// save previous entries for the event error check
		prev = veto;
		prevGoodEntry = i;
		veto.Clear();
	}
	if (makePlots)
	{
		TCanvas *c1 = new TCanvas("c1","full QDC",1600,1200);
		c1->Divide(8,4,0,0);
		for (int i=0; i<32; i++)
		{
			c1->cd(i+1);
			TVirtualPad *vpad1 = c1->cd(i+1); vpad1->SetLogy();
			hFullQDC[i]->Draw();
		}
		TCanvas *c2 = new TCanvas("c2","QDC thresholds",1600,1200);
		c2->Divide(8,4,0,0);
		for (int i=0; i<32; i++)
		{
			c2->cd(i+1);
			TVirtualPad *vpad2 = c2->cd(i+1); vpad2->SetLogy();
			hLowQDC[i]->GetXaxis()->SetRange(lower,upper);
			hLowQDC[i]->Draw();
			double ymax = hLowQDC[i]->GetMaximum();
			TLine *line = new TLine(thresh[i],0,thresh[i],ymax+10);
			line->SetLineColor(kRed);
			line->SetLineWidth(2.0);
			line->Draw();
		}
		TCanvas *c3 = new TCanvas("c3","multiplicity",800,600);
		c3->cd();
		c3->SetLogy();
		hMultip->Draw();

		char plotName[200];
		sprintf(plotName,"%s/veto-%i-qdc.pdf",outputDir.c_str(),runNum);
		c1->Print(plotName);

		sprintf(plotName,"%s/veto-%i-qdcThresh.pdf",outputDir.c_str(),runNum);
		c2->Print(plotName);

		sprintf(plotName,"%s/veto-%i-multip.pdf",outputDir.c_str(),runNum);
		c3->Print(plotName);

	}
	return thresholds;
}

void ProcessVetoData(TChain *vetoChain, vector<int> thresholds, string outputDir, bool errorCheckOnly)
{
	// QDC software threshold (obtained from VetoThreshFinder)
	int swThresh[32] = {0};
	for (int i = 0; i < (int)thresholds.size(); i+=2)
		swThresh[thresholds[i]] = thresholds[i+1];
	// cout << "ProcessVetoData is using these QDC thresholds: " << endl;
	// for (int i = 0; i < 32; i++)
	// 	cout << i << " " << swThresh[i] << endl;

	// LED variables
	int LEDMultipThreshold = 5;  // "multipThreshold" = "highestMultip" - "LEDMultipThreshold"
	int LEDSimpleThreshold = 10;  // used when LED frequency measurement is bad.
	int highestMultip=0, multipThreshold=0;
	double LEDWindow = 0.1;
	double LEDfreq=0, LEDperiod=0, LEDrms=0;
	bool badLEDFreq=false;
	bool IsLED = false, IsLEDPrev = false;
	bool LEDTurnedOff = false;
	int simpleLEDCount=0;
	bool useSimpleThreshold=false;
	bool firstLED=false;

	// Error check variables
	const int nErrs = 29; // error 0 is unused
	int SeriousErrorCount = 0;
	int TotalErrorCount = 0;
	vector<double> EntryTime;
	vector<bool> BadScalers;
	vector<int> EventError(nErrs), ErrorCount(nErrs);
	bool badEvent = false;
	// Specify which error types to print during the loop over events
	vector<int> SeriousErrors = {1, 13, 14, 18, 19, 20, 21, 22, 23, 24, 26};

	// muon ID variables
	bool LEDCut = true;	// actually an LED cut
	bool EnergyCut = false;
	vector<int> CoinType(32), CutType(32), PlaneTrue(32);

	// time variables
	double SBCOffset=0;
	bool ApproxTime = false;
	bool foundScalerJump = false;	// this is related to ApproxTime but should be kept separate.
	double deltaTadj = 0;
	double timeSBC=0, prevtimeSBC=0;
	long skippedEvents=0;
	long start=0, stop=0;
	double duration=0, livetime=0, xTime=0, x_deltaT=0, x_LEDDeltaT=0;
	double xTimePrev=0, x_deltaTPrev=0, xTimePrevLED=0, xTimePrevLEDSimple=0, TSdifference=0;
	double prevGoodTime=0, firstGoodScaler=0;
	int prevGoodEntry=0;
	bool foundFirst=false, foundFirstScaler=false;

	// initialize input data
	long vEntries = vetoChain->GetEntries();
	TTreeReader reader(vetoChain);
	TTreeReaderValue<unsigned int> vMult(reader, "mVeto");
	TTreeReaderValue<uint32_t> vBits(reader, "vetoBits");
	TTreeReaderValue<MGTBasicEvent> vEvt(reader,"vetoEvent");
	TTreeReaderValue<MJTRun> vRun(reader,"run");
	reader.SetEntry(0);
	int runNum = vRun->GetRunNumber();
	start = vRun->GetStartTime();
	stop = vRun->GetStopTime();
	duration = (double)(stop - start);
	reader.SetTree(vetoChain);  // resets the reader

	// MJVetoEvent variables, with run-based card numbers
	int card1=0, card2=0;
	SetCardNumbers(runNum,card1,card2);
	MJVetoEvent veto(card1,card2);
	MJVetoEvent first, prev, last, prevLED, out;

	// initialize output file
	char outputFile[200];
	sprintf(outputFile,"%s/veto_run%i.root",outputDir.c_str(),runNum);
	TFile *RootFile = new TFile(outputFile, "RECREATE");
	TTree *vetoTree = new TTree("vetoTree","MJD Veto Events");
	// event info
	vetoTree->Branch("run",&runNum);
	vetoTree->Branch("vetoEvent","MJVetoEvent",&out,32000,1);
	// LED variables
	vetoTree->Branch("LEDfreq",&LEDfreq);
	vetoTree->Branch("LEDrms",&LEDrms);
	vetoTree->Branch("multipThreshold",&multipThreshold);
	vetoTree->Branch("highestMultip",&highestMultip);
	vetoTree->Branch("LEDWindow",&LEDWindow);
	vetoTree->Branch("LEDMultipThreshold",&LEDMultipThreshold);
	vetoTree->Branch("LEDSimpleThreshold",&LEDSimpleThreshold);
	vetoTree->Branch("useSimpleThreshold",&useSimpleThreshold);
	// time variables
	vetoTree->Branch("start",&start,"start/L");
	vetoTree->Branch("stop",&stop,"stop/L");
	vetoTree->Branch("duration",&duration);
	vetoTree->Branch("livetime",&livetime);
	vetoTree->Branch("xTime",&xTime);
	vetoTree->Branch("timeSBC",&timeSBC);
	vetoTree->Branch("x_deltaT",&x_deltaT);
	vetoTree->Branch("x_LEDDeltaT",&x_LEDDeltaT);
	// muon ID variables
	vetoTree->Branch("CoinType",&CoinType);
	vetoTree->Branch("CutType",&CutType);
	vetoTree->Branch("PlaneTrue",&PlaneTrue);
	// error variables
	vetoTree->Branch("badEvent",&badEvent);
	vetoTree->Branch("EventErrors",&EventError);
	vetoTree->Branch("ErrorCount",&ErrorCount);

	// Error tree
	TTree *skipTree = new TTree("skipTree","skipped veto events");
	skipTree->Branch("badEvent",&badEvent);
	skipTree->Branch("EventErrors",&EventError);
	skipTree->Branch("ErrorCount",&ErrorCount);
	skipTree->Branch("run",&runNum);
	skipTree->Branch("vetoEvent","MJVetoEvent",&out,32000,1);


	// ================= 1st loop over veto entries  ==============
	// Measure the LED frequency, find the first good entry,
	// highest multiplicity, and SBC offset.

	TH1D *LEDDeltaT = new TH1D("LEDDeltaT","LEDDeltaT",100000,0,100); // 0.001 sec/bin
	while(reader.Next())
	{
		long i = reader.GetCurrentEntry();
		int run = vRun->GetRunNumber();

		veto.SetSWThresh(swThresh);
		veto.WriteEvent(i,&*vRun,&*vEvt,*vBits,run,true);
		if (!veto.GetBadScaler()) {
			BadScalers.push_back(0);
			xTime = veto.GetTimeSec();
		}
		else {
			BadScalers.push_back(1);
			xTime = ((double)i / vEntries) * duration; // this breaks if we have corrupted duration
		}
		EntryTime.push_back(xTime);

		if (foundFirst && veto.GetError(1)) {
			foundFirst = false;
		}
		if (!foundFirstScaler && !veto.GetError(4)){
			foundFirstScaler = true;
			firstGoodScaler = veto.GetTimeSec();
		}
		if (CheckEventErrors(veto,prev,first,prevGoodEntry)){
			skippedEvents++;
			// do end of loop resets before continuing
			prev = veto;
			prevGoodTime = xTime;
			prevGoodEntry = i;
			veto.Clear();
			continue;
		}
		if (!foundFirst && veto.GetTimeSBC() > 0 && veto.GetTimeSec() > 0 && !veto.GetError(4)) {
			first = veto;
			foundFirst = true;
		}
		if (veto.GetMultip() > highestMultip) {
			highestMultip = veto.GetMultip();
		}
		// cout << veto.GetMultip() << endl;
		if (veto.GetMultip() > LEDSimpleThreshold) {
			LEDDeltaT->Fill(veto.GetTimeSec()-prev.GetTimeSec());
			simpleLEDCount++;
		}
		// end of loop
		prev = veto;
		prevGoodTime = xTime;
		prevGoodEntry = i;
		veto.Clear();
	}

	// =====================  Run-level checks =====================

	SBCOffset = first.GetTimeSBC() - first.GetTimeSec();

	// 27. QDC threshold not found
	// 28. No events above QDC threshold
	for (int i=0; i < 32; i++) if (swThresh[i] == 9999) {
		EventError[27]=true;
		ErrorCount[27]++;
		cout << "Warning: Couldn't find QDC threshold for panel " << i << ". Set to 9999\n";
	}

	// duration and live time
	// TODO: implement an alternate way of finding the run stop packet "fStopTime".
	if (duration <= 0 || duration > 9999) {
		duration = prevGoodTime-firstGoodScaler;
		cout << "Warning: bad run duration.  start: " << start << " stop: " << stop
			  << "  Using last good timestamp: " << prevGoodTime-firstGoodScaler << "\n";
	}
	livetime = duration - (first.GetTimeSec() - firstGoodScaler);
	cout << "Veto livetime: " << livetime << " seconds\n";
	if (livetime < 0) {
		cout << "Warning: corrupted livetime.  Setting it to 0.\n";
		livetime=0;
	}

	// set LED multiplicity threshold
	multipThreshold = highestMultip - LEDMultipThreshold;
	if (multipThreshold < 0) multipThreshold = 0;

	// find LED frequency, and try to adjust if we have a short run (only a few LED events)
	int dtEntries = LEDDeltaT->GetEntries();
	if (dtEntries > 0) {
		int maxbin = LEDDeltaT->GetMaximumBin();
		LEDDeltaT->GetXaxis()->SetRange(maxbin-100,maxbin+100); // looks at +/- 0.1 seconds of max bin.
		LEDrms = LEDDeltaT->GetRMS();
		LEDfreq = 1/LEDDeltaT->GetMean();
	}
	else {
		cout << "Warning! No multiplicity > " << LEDSimpleThreshold << " events.  LED may be off.\n";
		LEDrms = 9999;
		LEDfreq = 9999;
		badLEDFreq = true;
	}
	LEDperiod = 1/LEDfreq;
	delete LEDDeltaT;
	if (LEDperiod > 9 || vEntries < 100)
	{
		cout << "Warning: Short run.\n";
		if (simpleLEDCount > 3) {
			cout << "   From histo method, LED freq is " << LEDfreq
				 << "  Using approximate rate: " << simpleLEDCount/duration << endl;
			LEDperiod = duration/simpleLEDCount;
			useSimpleThreshold=true;
		}
		else {
			LEDperiod = 9999;
			badLEDFreq = true;
		}
	}
	// 26. LED frequency very low/high, corrupted, or LED's off.
	if (LEDperiod > 20 || LEDperiod < 0 || badLEDFreq) {
		ErrorCount[26]++;
		EventError[26] = true;
	}

	// ========= 2nd loop over entries - Error checks =========
	// this time we don't skip anything, and count up each type of error.

	reader.SetTree(vetoChain); // reset the reader
	std::fill(EventError.begin(), EventError.end(), 0); // reset error bools
	while(reader.Next())
	{
		long i = reader.GetCurrentEntry();
		int run = vRun->GetRunNumber();

		veto.SetSWThresh(swThresh);
		veto.WriteEvent(i,&*vRun,&*vEvt,*vBits,run,true);
		CheckEventErrors(veto,prev,first,prevGoodEntry,EventError);
		for (int j=0; j<nErrs; j++) if (EventError[j]==1) ErrorCount[j]++;

		// xTime = FindEventTime(veto,SBCOffset,EntryTime,BadScalers,timeSBC);

		// Find event time.
		if (!veto.GetBadScaler())
		{
			xTime = veto.GetTimeSec();
			if(run > 8557 && veto.GetTimeSBC() < 2000000000)
				timeSBC = veto.GetTimeSBC() - SBCOffset;
		}
		else if (run > 8557 && run < 45000000 && veto.GetTimeSBC() < 2000000000)
			xTime = veto.GetTimeSBC() - SBCOffset;
		else {
			xTime = InterpTime(i,EntryTime,BadScalers);
			ErrorCount[25]++;
		}
		EntryTime[i] = xTime;	// replace entry with the most accurate one

		// Print errors to screen
		bool PrintError = false;
		for (auto i : SeriousErrors){
			if (EventError[i]) {
				PrintError = true;
				break;
			}
		}
		if (PrintError)
		// if (i==717)
		{
			cout << "\nSerious error found in entry " << i << ":\n";

			// debug block (don't delete!) used to compare with original vetoCheck code
			// cout << veto.GetTimeSec() << "\n" 	// STime
			// 	  << prev.GetTimeSec() << "\n" 	// STimePrev
			// 	  << timeSBC << "\n"					// SBCTime
			// 	  << SBCOffset << "\n"				// SBCOffset
			// 	  << prevtimeSBC << "\n"			// SBCTimePrev
			// 	  << veto.GetError(1) << "\n"
			// 	  << i  << "\n"						// entry
			// 	  << first.GetEntry() << "\n"		// firstGoodEntry
			// 	  << "error 18: " << EventError[18] << "\n"			// Error[18]
			// 	  << veto.GetSEC() << "\n"
			// 	  << prev.GetSEC() << "\n"
			// 	  << prevGoodEntry << "\n"			//  EventNumPrev_good
			// 	  << "error 20: " << EventError[20] << "\n";
			// long qdc1=veto.GetQEC(), qdc2=veto.GetQEC2(), scaler=veto.GetSEC();
			// cout << i << "  " << scaler << "  " << qdc1 << "  " << qdc2
				//   << " q1-s " << qdc1-scaler << " q2-s: " << qdc2-scaler << "  "
				//   << veto.GetError(11) << "  " << veto.GetError(12) << endl;

			if (EventError[1]) {
				cout << "EventError[1] Missing Packet."
					 << "  Scaler index " << veto.GetScalerIndex()
					 << "  Scaler Time " << veto.GetTimeSec()
					 << "  SBC Time " << veto.GetTimeSBC() << "\n";
			}
			if (EventError[13]) {
				cout << "EventError[13] ORCA packet indexes of QDC1 and Scaler differ by more than 2."
					 << "\n    Scaler Index " << veto.GetScalerIndex()
					 << "  QDC1 Index " << veto.GetQDC1Index()
					 << "\n    Previous scaler Index " << prev.GetScalerIndex()
					 << "  Previous QDC1 Index " << prev.GetQDC1Index() << endl;
			}
			if (EventError[14]) {
				cout << "EventError[14] ORCA packet indexes of QDC2 and Scaler differ by more than 2."
					 << "\n    Scaler Index " << veto.GetScalerIndex()
					 << "  QDC2 Index " << veto.GetQDC2Index()
					 << "\n    Previous scaler Index " << prev.GetScalerIndex()
					 << "  Previous QDC2 Index " << prev.GetQDC2Index() << endl;
			}
			if (EventError[18]) {
				cout << "EventError[18] Scaler/SBC Desynch."
					 << "\n    DeltaT (adjusted) " << veto.GetTimeSec() - timeSBC - TSdifference
					 << "  DeltaT " << veto.GetTimeSec() - timeSBC
					 << "\n    Prev TSdifference " << TSdifference
					 << "  Scaler DeltaT " << veto.GetTimeSec()-prev.GetTimeSec()
					 << "\n    Scaler Index " << veto.GetScalerIndex()
					 << "  Previous Scaler Index " << prev.GetScalerIndex()
					 << "\n    Scaler Time " << veto.GetTimeSec()
					 << "  SBC Time " << timeSBC << "\n";
				TSdifference = veto.GetTimeSec() - timeSBC;
			}
			if (EventError[19]) {
				cout << "EventError[19] Scaler Event Count Reset. "
					 << "\n    Scaler Index " << veto.GetScalerIndex()
					 << "  SEC " << veto.GetSEC()
					 << "  Previous SEC " << prev.GetSEC() << "\n";
			}
			if (EventError[20]) {
				cout << "EventError[20] Scaler Event Count Jump."
					 << "    xTime " << xTime
					 << "  Scaler Index " << veto.GetScalerIndex()
					 << "\n    SEC " << veto.GetSEC()
					 << "  Previous SEC " << prev.GetSEC() << "\n";
			}
			if (EventError[21]) {
				cout << "EventError[21] QDC1 Event Count Reset."
					 << "\n    Scaler Index " << veto.GetScalerIndex()
					 << "  QEC1 " << veto.GetQEC()
					 << "  Previous QEC1 " << prev.GetQEC() << "\n";
			}
			if(EventError[22]) {
				cout << "EventError[22] QDC 1 Event Count Jump."
					 << "\n    xTime " << xTime
					 << "  QDC 1 Index " << veto.GetQDC1Index()
					 << "  QEC 1 " << veto.GetQEC()
					 << "  Previous QEC 1 " << prev.GetQEC() << "\n";
			}
			if (EventError[23]) {
				cout << "EventError[23] QDC2 Event Count Reset."
					 << "\n    Scaler Index " << veto.GetScalerIndex()
					 << "  QEC2 " << veto.GetQEC2()
					 << "  Previous QEC2 " << prev.GetQEC2() << "\n";
			}
			if(EventError[24]) {
				cout << "EventError[24] QDC 2 Event Count Jump."
					 << "\n    xTime " << xTime
					 << "  QDC 2 Index " << veto.GetQDC2Index()
					 << "  QEC 2 " << veto.GetQEC2()
					 << "  Previous QEC 2 " << prev.GetQEC2() << "\n";
			}
		}
		// end of event resets
		TSdifference = veto.GetTimeSec() - timeSBC;
		prevtimeSBC = timeSBC;
		timeSBC = 0;
		prev = veto;
		last = prev;
		prevGoodEntry = i;
		veto.Clear();
		std::fill(EventError.begin(), EventError.end(), 0);
	}

	// Calculate total errors and total serious errors
	for (int i = 1; i < nErrs; i++)
	{
			// Ignore 10 and 11, these will always be present
			// as long as the veto counters are not reset at the beginning of runs.
			if (i != 10 && i != 11) TotalErrorCount += ErrorCount[i];

			// count up the serious errors in the vector
			for (auto j : SeriousErrors)
				if (i == j)
					SeriousErrorCount += ErrorCount[i];
	}
	cout << "=================== Veto Error Report ===================\n";
	cout << "Serious errors found :: " << SeriousErrorCount << endl;
	if (SeriousErrorCount > 0)
	{
		cout << "Total Errors : " << TotalErrorCount << endl;
		if (duration != livetime)
			cout << "Run duration (" << duration << " sec) doesn't match live time: " << livetime << endl;

		for (int i = 1; i < nErrs; i++)
		{
			if (ErrorCount[i] > 0 && (i!=7 && i!=10 && i!=11))
			{
				if (i != 26)
					cout << "  Error[" << i <<"]: " << ErrorCount[i] << " events ("
						 << 100*(double)ErrorCount[i]/vEntries << " %)\n";
		 		else if (i == 26)
				{
					cout << "  Error[26]: Bad LED rate: " << LEDfreq << "  Period: " << LEDperiod << endl;
					if (LEDperiod > 0.1 && (abs(duration/LEDperiod) - simpleLEDCount) > 5)
					{
						cout << "   Simple LED count: " << simpleLEDCount
							 << "  Expected: " << (int)(duration/LEDperiod) << endl;
					}
				}
			}
		}
		cout << "For reference, \"serious\" error types are: ";
		for (auto i : SeriousErrors) cout << i << " ";
		cout << "\nPlease report these to the veto group.\n";
	}
	if (errorCheckOnly) return;

	// ========= 3rd loop over veto entries - Find muons! Write ROOT output! =========

	cout << "================= Scanning for muons ... ================\n";
	reader.SetTree(vetoChain); // reset the reader
	prev.Clear();
	skippedEvents = 0;
	printf("Highest multiplicity found: %i.  Using LED threshold: %i\n",highestMultip,multipThreshold);

	while(reader.Next())
	{
		long i = reader.GetCurrentEntry();
		int run = vRun->GetRunNumber();

		veto.SetSWThresh(swThresh);
		veto.WriteEvent(i,&*vRun,&*vEvt,*vBits,run,true);
		CheckEventErrors(veto,prev,first,prevGoodEntry,EventError);
		for (int j=0; j<nErrs; j++) if (EventError[j]==1) ErrorCount[j]++;

		// Calculate time of event.
		// TODO: implement an estimate of the error when alternate methods are used.
		ApproxTime = false;
		if (!veto.GetBadScaler())
		{
			xTime = veto.GetTimeSec();
			if(run > 8557 && veto.GetTimeSBC() < 2000000000)
				timeSBC = veto.GetTimeSBC() - SBCOffset;
		}
		else if (run > 8557 && veto.GetTimeSBC() < 2000000000) {
			timeSBC = veto.GetTimeSBC() - SBCOffset;
			xTime = timeSBC;
		}
		else {
			xTime = InterpTime(i,EntryTime,BadScalers);
			ApproxTime = true;
			EventError[25]=true;
		}
		// Scaler jump handling:
		// When the error is initally found, we adjust xTime by "DeltaT (adjusted)" from EventError[18].
		// If we find a scaler jump in one event, we force the scaler and SBC times to match
		// for the rest of the run, or until they come back in sync on their own.
		deltaTadj = veto.GetTimeSec() - timeSBC - TSdifference;
		if (EventError[18]) {
			foundScalerJump=true;
			xTime = xTime - deltaTadj;
		}
		if (foundScalerJump) {
			ApproxTime = true;
			xTime = xTime - TSdifference;
		}
		if (fabs(TSdifference) < 0.001)	// SBC is accurate to microseconds
			foundScalerJump = false;
		// debug block for a run w/ jump at i==248 (don't delete!)
		// if (i==247||i==248||i==249||i==250)
			// printf("%i  e18: %i  %-6.3f  %-6.3f  %-6.3f  %-6.3f  %-6.3f\n", 		i,EventError[18],veto.GetTimeSec(),timeSBC,xTime,TSdifference,deltaTadj);

		// Skip events after the event time is calculated,
		// so we still properly reset for the next event.
		// These events also go into the skipTree.
		badEvent = CheckEventErrors(veto,prev,first,prevGoodEntry,EventError);
		if (badEvent)
		{
			skippedEvents++;
			// do the end-of-event reset
			TSdifference = veto.GetTimeSec() - timeSBC;
			prevtimeSBC = timeSBC;
			timeSBC = 0;
			prev = veto;
			last = prev;
			prevGoodEntry = i;
			skipTree->Fill();
			veto.Clear();
			continue;
		}

		// LED Cut
		LEDCut = false;
		LEDTurnedOff = ErrorCount[26];
		x_deltaT = xTime - xTimePrevLED;

		// TODO: right now this enforces a hard multiplicity cut.
		// Need to evaluate if the frequency-based multiplicity cut will
		// recover high-multiplicity muon events without mis-identifying LED's as muons.
		if (veto.GetMultip() < multipThreshold && !LEDTurnedOff) {
			// IsLED = false;
			LEDCut = true;
		}
		else if (LEDTurnedOff)
			LEDCut = true;

		/*
		// if (!LEDTurnedOff && !badLEDFreq && fabs(LEDperiod - x_deltaT) < LEDWindow && veto.GetMultip() > multipThreshold)
		// {
		// 	LEDCut = false;
		// 	IsLED = true;
		// }
		// // almost missed a high-multiplicity event somehow ...
		// // often due to skipping previous events.
		// else if (!LEDTurnedOff && !badLEDFreq && fabs(LEDperiod - x_deltaT) >= (LEDperiod - LEDWindow)
		// 				&& veto.GetMultip() > multipThreshold && i > 1)
		// {
		// 	LEDCut = false;
		// 	IsLED = true;
		// 	cout << "Almost missed LED:\n";
		// 	printf("Current: %-3li  m %-3i LED? %i t %-6.2f LEDP %-5.2f  XDT %-6.2f LEDP-XDT %-6.2f\n"
		// 		,i,veto.GetMultip(),IsLED,xTime,LEDperiod,x_deltaT,LEDperiod-x_deltaT);
		// }
		// else LEDCut = true;
		// // Grab first LED
		// if (!LEDTurnedOff && !firstLED && veto.GetMultip() > multipThreshold) {
		// 	printf("Found first LED.  i %-2li m %-2i t %-5.2f thresh:%i  LEDoff:%i\n\n",i,veto.GetMultip(),xTime,multipThreshold,LEDTurnedOff);
		// 	IsLED=true;
		// 	firstLED=true;
		// 	LEDCut=false;
		// 	x_deltaT = -1;
		// }
		// // If frequency measurement is bad, revert to standard multiplicity cut
		// if (badLEDFreq && veto.GetMultip() >= LEDSimpleThreshold){
		// 	IsLED = true;
		// 	LEDCut = false;
		// }
		// // Simple x_LEDDeltaT uses the multiplicity-only threshold, veto.GetMultip() > multipThreshold.
		// x_LEDDeltaT = xTime - xTimePrevLEDSimple;
		//
		// // If LED is off, all events pass time cut.
		// if (LEDTurnedOff) {
		// 	IsLED = false;
		// 	LEDCut = true;
		// }
		// // Check output
		// printf("%-3li  m %-3i LED? %i t %-6.2f LEDP %-5.2f  XDT %-6.2f LEDP-XDT %-6.2f\n"
		// 	,i,veto.GetMultip(),IsLED,xTime,LEDperiod,x_deltaT,LEDperiod-x_deltaT);
		*/

    	// Energy (Gamma) Cut
    	// The measured muon energy threshold is QDC = 500.
    	// Set TRUE if at least TWO panels are over 500.
		EnergyCut = false;
    	int over500Count = 0;
    	for (int q = 0; q < 32; q++) {
    		if (veto.GetQDC(q) > 500)
    			over500Count++;
    	}
    	if (over500Count >= 2) EnergyCut = true;

		// debug block (don't delete): check entries
		printf("Entry %li  Time %-6.2f  QDC %-5i  Mult %i  Ov500 %i  Loff? %i  LEDCut %i  ECut %i  \n",i,xTime,veto.GetTotE(),veto.GetMultip(),over500Count,LEDTurnedOff,LEDCut,EnergyCut);


		// Muon Identification:
		// Use EnergyCut, LEDCut, and the Hit Pattern to identify them sumbitches.
		for (int k = 0; k < 12; k++) PlaneTrue[k] = 0;
		for (int k = 0; k < 32; k++) {
			if (veto.GetQDC(k) > veto.GetSWThresh(k)) {
				if (PlaneMap(k,run)==0)       PlaneTrue[0]=1;  // 0: Lower Bottom
				else if (PlaneMap(k,run)==1)  PlaneTrue[1]=1;  // 1: Upper Bottom
				else if (PlaneMap(k,run)==2)  PlaneTrue[2]=1;  // 3: Top Inner
				else if (PlaneMap(k,run)==3)  PlaneTrue[3]=1;  // 4: Top Outer
				else if (PlaneMap(k,run)==4)  PlaneTrue[4]=1;  // 5: North Inner
				else if (PlaneMap(k,run)==5)  PlaneTrue[5]=1;  // 6: North Outer
				else if (PlaneMap(k,run)==6)  PlaneTrue[6]=1;  // 7: South Inner
				else if (PlaneMap(k,run)==7)  PlaneTrue[7]=1;  // 8: South Outer
				else if (PlaneMap(k,run)==8)  PlaneTrue[8]=1;  // 9: West Inner
				else if (PlaneMap(k,run)==9)  PlaneTrue[9]=1;  // 10: West Outer
				else if (PlaneMap(k,run)==10) PlaneTrue[10]=1; // 11: East Inner
				else if (PlaneMap(k,run)==11) PlaneTrue[11]=1; // 12: East Outer
				else if (PlaneMap(k,run)==-1)
					cout << "Error: Panel " << k << " was not installed for this run and should not be giving counts above threshold.\n";
			}
		}
		for (int r = 0; r < 32; r++) {CoinType[r]=0; CutType[r]=0;}
		if (LEDCut && EnergyCut)
		{
			CoinType[0] = true;
			int type = 0;
			bool a=0,b=0,c=0;

			// debug block (don't delete!)
			// cout << "\nQDC-panel-plane: ";
			// for (int i=0; i<32; i++) {
			// 	int qdc=0;
			// 	if (veto.GetQDC(i) > veto.GetSWThresh(i)) {
			// 		qdc = veto.GetQDC(i);
			// 		cout << i << "-" << i+1 << "-p" << PlaneMap(i,run) << ":" << qdc << "  ";}
			// }
			// cout << endl;
			// cout << "Planes:\n";
			// cout << "0: Lower Bottom  1: Upper Bottom\n"
			// 	  << "2: Top Inner     3: Top Outer\n"
			// 	  << "4: North Inner   5: North Outer\n"
			// 	  << "6: South Inner   7: South Outer\n"
			// 	  << "8: West Inner    9: West Outer\n"
			// 	  << "10: East Inner   11: East Outer\n";
			// for (int i=0; i<12; i++) cout << "p" << i << ":" << PlaneTrue[i] << "  ";
			// cout << endl;

			// Type 1: vertical muons
			if (PlaneTrue[0] && PlaneTrue[1] && PlaneTrue[2] && PlaneTrue[3]) {
				CoinType[1]=true;
				a=true;
				type=1;
			}
			// Type 2: (both planes of a side) + both bottom planes
			if ((PlaneTrue[0] && PlaneTrue[1]) &&
				((PlaneTrue[4] && PlaneTrue[5]) || (PlaneTrue[6] && PlaneTrue[7])
				  || (PlaneTrue[8] && PlaneTrue[9])
				  || (PlaneTrue[10] && PlaneTrue[11]) ) )
			{
				CoinType[2] = true;
				b=true;
				type = 2;
			}
			// Type 3: (both top planes) + (both planes of a side)
			if ((PlaneTrue[2] && PlaneTrue[3]) && ((PlaneTrue[4] && PlaneTrue[5]) || (PlaneTrue[6] && PlaneTrue[7])
					|| (PlaneTrue[8] && PlaneTrue[9]) || (PlaneTrue[10] && PlaneTrue[11]))) {
				CoinType[3] = true;
				c=true;
				type = 3;
			}
			// Type 4: compound hit (combination of types 1-3)
			if ((a && b)||(a && c)||(b && c)) type = 4;
			cout << "a:" << a << "  b:" << b << "  c:" << c << endl;

			char hitType[200];
			if (type==0) sprintf(hitType,"2+ panels");
			if (type==1) sprintf(hitType,"vertical");
			if (type==2) sprintf(hitType,"side+bottom");
			if (type==3) sprintf(hitType,"top+sides");
			if (type==4) sprintf(hitType,"compound");

			// print the details of the hit
			printf("Hit: %-12s Entry %-4li Time %-6.2f  QDC %-5i  Mult %i  Ov500 %i  LEDoff %i  ApxT %i\n", hitType,i,xTime,veto.GetTotE(),veto.GetMultip(),over500Count,LEDTurnedOff,ApproxTime);
		}

		// Output

		CutType[0] = LEDTurnedOff;
		CutType[1] = EnergyCut;
		CutType[2] = ApproxTime;
		CutType[3] = LEDCut;
		CutType[4] = IsLED;
		CutType[5] = firstLED;
		CutType[6] = badLEDFreq;

		out = veto;
		vetoTree->Fill();

		// Reset for next entry

		// resets used by error check (don't delete!)
		TSdifference = veto.GetTimeSec() - timeSBC;
		prevtimeSBC = timeSBC;
		timeSBC = 0;
		prev = veto;
		last = prev;
		prevGoodEntry = i;

		// resets used by muon finder (may need revision)
		if (IsLED) {
			prevLED = veto;
			xTimePrevLED = xTime;
		}
		if (veto.GetMultip() > multipThreshold) {
			xTimePrevLEDSimple = xTime;
		}
		IsLEDPrev = IsLED;
		xTimePrev = xTime;
		x_deltaTPrev = x_deltaT;

		veto.Clear();
	}
	if (skippedEvents > 0) printf("ProcessVetoData skipped %li of %li entries.\n",skippedEvents,vEntries);

	vetoTree->Write("",TObject::kOverwrite);
	skipTree->Write("",TObject::kOverwrite);
	RootFile->Close();
	cout << "Wrote ROOT file: " << outputFile << endl;
}
// ====================================================================================
// ====================================================================================

void SetCardNumbers(int runNum, int &card1, int &card2)
{
	if (runNum > 45000000){
		card1 = 11;
		card2 = 18;
	}
	else {
		card1 = 13;
		card2 = 18;
	}
}

int FindPanelThreshold(TH1D *qdcHist, int threshVal, int panel, int runNum)
{
	if (runNum > 45000000 && panel > 23)
		return 9999;

	int firstNonzeroBin = qdcHist->FindFirstBinAbove(1,1);
	qdcHist->GetXaxis()->SetRange(firstNonzeroBin-10,firstNonzeroBin+50);
	//qdcHist->GetXaxis()->SetRangeUser(0,500); //alternate method of finding pedestal
	int bin = qdcHist->GetMaximumBin();
	if (firstNonzeroBin == -1) return 9999;
	double xval = qdcHist->GetXaxis()->GetBinCenter(bin);
	return xval+threshVal;
}

double InterpTime(int entry, vector<double> times, vector<bool> badScaler)
{
	if (times.size() != badScaler.size())
	{
		cout << "Vectors are different sizes!\n";
		if (entry >= (int)times.size())
			cout << "Entry is larger than number of entries in vector!\n";
		return -1;
	}
	double iTime = 0, lower=0, upper=0;
	if (!badScaler[entry]) iTime = times[entry];
	else
	{
		for (int i = entry; i < (int)times.size(); i++)
			if (badScaler[i] == 0) { upper = times[i]; break; }
		for (int j = entry; j > 0; j--)
			if (badScaler[j] == 0) { lower = times[j]; break; }
		iTime = (upper + lower)/2.0;
	}
	return iTime;
}

int PlaneMap(int qdcChan, int runNum)
{
	// For tagging plane-based coincidences.
	// This uses a zero-indexed map: {panel number, plane number}
	map<int,int> planes;

	// Geometric planes:
	// 0: Lower Bottom 	1: Upper Bottom
	// 2: Top Inner		3: Top Outer
	// 4: North Inner		5: North Outer
	// 6: South Inner 	7: South Outer
	// 8: West Inner		9: West Outer
	// 10: East Inner 	11: East Outer

	// 32-panel config (default) - began 7/10/15
	if (runNum < 45000000 && runNum >= 3057) {
		map<int,int> tempMap =
		{{0,0},  {1,0},  {2,0},  {3,0}, {4,0}, {5,0},
		 {6,1},  {7,1},  {8,1},  {9,1}, {10,1}, {11,1},
		 {12,8}, {13,8}, {14,9}, {15,5},
		 {16,5}, {17,3}, {18,3}, {19,4},
		 {20,2}, {21,2}, {22,9}, {23,4},
		 {24,6}, {25,7}, {26,6}, {27,7},
		 {28,10},{29,11},{30,10},{31,11}};
		planes = tempMap;
	}
	// 1st prototype config (24 panels)
 	else if (runNum >= 45000509 && runNum <= 45004116) {
		map<int,int> tempMap =
		{{0,0},  {1,0},  {2,0},  {3,0}, {4,0}, {5,0},
		 {6,1},  {7,1},  {8,1},  {9,1}, {10,1}, {11,1},
		 {12,2}, {13,2}, {14,9}, {15,5},
		 {16,5}, {17,3}, {18,3}, {19,4},
		 {20,8}, {21,8}, {22,9}, {23,4},
		 {24,-1},{25,-1},{26,-1},{27,-1},
		 {28,-1},{29,-1},{30,-1},{31,-1}};
		planes = tempMap;
	}
	// 2nd prototype config (24 panels)
 	else if (runNum >= 45004117 && runNum <= 45008659) {
		map<int,int> tempMap =
		{{0,0},  {1,0},  {2,0},  {3,0}, {4,0}, {5,0},
		 {6,1},  {7,1},  {8,1},  {9,1}, {10,1}, {11,1},
		 {12,8}, {13,8}, {14,9}, {15,5},
		 {16,5}, {17,3}, {18,3}, {19,4},
		 {20,2}, {21,2}, {22,9}, {23,4},
		 {24,-1},{25,-1},{26,-1},{27,-1},
		 {28,-1},{29,-1},{30,-1},{31,-1}};
		planes = tempMap;
	}
	// 1st module 1 (P3JDY) config, 6/24/15 - 7/7/15
	else if (runNum > 0 && runNum <=3056){
		map<int,int> tempMap =
		{{0,0},  {1,0},  {2,0},  {3,0}, {4,0}, {5,0},
		 {6,1},  {7,1},  {8,1},  {9,1}, {10,1}, {11,1},
		 {12,8}, {13,8}, {14,9}, {15,5},
		 {16,5}, {17,3}, {18,3}, {19,4},
		 {20,2}, {21,2}, {22,9}, {23,4},
		 {24,-1},{25,-1},{26,-1},{27,-1},
		 {28,-1},{29,-1},{30,-1},{31,-1}};
		planes = tempMap;
	}
	// If the panel is not installed for this run, return -1.
	else return -1;

	// Otherwise, return the plane index for this panel.
	int plane = -1;
	auto search = planes.find(qdcChan);
	if(search != planes.end()) plane=search->second;
	return plane;
}

// This is overloaded so we don't have to use the error vector if we don't need it
bool CheckEventErrors(MJVetoEvent veto, MJVetoEvent prev, MJVetoEvent first, long prevGoodEntry)
{
	vector<int> ErrorVec(29);
	std::fill(ErrorVec.begin(), ErrorVec.end(), 0);
	return CheckEventErrors(veto,prev,first,prevGoodEntry,ErrorVec);
}
bool CheckEventErrors(MJVetoEvent veto, MJVetoEvent prev, MJVetoEvent first, long prevGoodEntry, vector<int> &ErrorVec)
{
	// Returns skip = false if the event is analyzable (either clean, or a workaround exists)
	bool skip = false;

	/*
	Recoverable:
	Actionable:
		LED Off -> Output string Dave can grep for, send veto experts an email
		No QDC events -> Need to have people check if a panel went down
		High error rate (desynchs, bad scalers, etc) -> Need to try and restart data taking
	Diagnostic:

	// Event-level error checks ('s' denotes setting skip=true)
	// s 1. Missing channels (< 32 veto datas in event)
	// s 2. Extra Channels (> 32 veto datas in event)
	// s 3. Scaler only (no QDC data)
	//   4. Bad Timestamp: FFFF FFFF FFFF FFFF
	// s 5. QDCIndex - ScalerIndex != 1 or 2
	// s 6. Duplicate channels (channel shows up multiple times)
	//   7. HW Count Mismatch (SEC - QEC != 1 or 2)
	//   8. MJTRun run number doesn't match input file
	// s 9. MJTVetoData cast failed (missing QDC data)
	//   10. Scaler EventCount doesn't match ROOT entry
	//   11. Scaler EventCount doesn't match QDC1 EventCount
	//   12. QDC1 EventCount doesn't match QDC2 EventCount
	// s 13. Indexes of QDC1 and Scaler differ by more than 2
	// s 14. Indexes of QDC2 and Scaler differ by more than 2
	//   15. Indexes of either QDC1 or QDC2 PRECEDE the scaler index
	//   16. Indexes of either QDC1 or QDC2 EQUAL the scaler index
	//   17. Unknown Card is present.
	// s 18. Scaler & SBC Timestamp Desynch.
	// s 19. Scaler Event Count reset.
	// s 20. Scaler Event Count increment by > +1.
	// s 21. QDC1 Event Count reset.
	// s 22. QDC1 Event Count increment by > +1.
	// s 23. QDC2 Event Count reset.
	// s 24. QDC2 Event Count increment > +1.
	//   25. Used interpolated time

	// Run-level error checks (not checked in this function)
	// 26. LED frequency very low/high, corrupted, or LED's off.
	// 27. QDC threshold not found
	// 28. No events above QDC threshold
	*/

	vector<int> EventErrors(29);
	std::fill(EventErrors.begin(), EventErrors.end(), 0);

	// Errors 1-18 are checked automatically when we call MJVetoEvent::WriteEvent
	for (int i=0; i < 18; i++) EventErrors[i] = veto.GetError(i);

	for (int q=0; q<18; q++) {
		if (EventErrors[q]==1 && (q==1||q==2||q==3||q==5||q==6||q==9||q==13||q==14))
			skip = true;
	}

	// condition to use interpolated time
	if (veto.GetBadScaler() && ((veto.GetRun() < 8557 || veto.GetRun() < 4500000) && veto.GetTimeSBC() > 2000000000))
		EventErrors[25] = true;

	int entry = veto.GetEntry();
	int firstGoodEntry = first.GetEntry();

	// assign EventErrors to the output vetor
	ErrorVec = EventErrors;

	// debug: print error vector
	// for (int i=0; i < 28; i++) cout << i << ":" << EventErrors[i] << "  ";
	// cout << endl;

	// return if we haven't found the first good entry yet
	if (firstGoodEntry == -1) return skip;

	double SBCOffset = first.GetTimeSBC() - first.GetTimeSec();
	double timeSBC = veto.GetTimeSBC() - SBCOffset;
	double prevtimeSBC = prev.GetTimeSBC() - SBCOffset;

	if (veto.GetTimeSec() > 0 && timeSBC > 0 && SBCOffset !=0
			&& !veto.GetError(1) && entry > firstGoodEntry
			&& fabs((veto.GetTimeSec() - prev.GetTimeSec())-(timeSBC - prevtimeSBC)) > 2)
				EventErrors[18] = true;

	if (veto.GetSEC() == 0 && entry != 0 && entry > firstGoodEntry)
		EventErrors[19] = true;

	if (abs(veto.GetSEC() - prev.GetSEC()) > entry-prevGoodEntry
			&& entry > firstGoodEntry && veto.GetSEC()!=0)
		EventErrors[20] = true;

	if (veto.GetQEC() == 0 && entry != 0 && entry > firstGoodEntry && !veto.GetError(1))
		EventErrors[21] = true;

	if (abs(veto.GetQEC() - prev.GetQEC()) > entry-prevGoodEntry
			&& entry > firstGoodEntry && veto.GetQEC() != 0)
		EventErrors[22] = true;

	if (veto.GetQEC2() == 0 && entry != 0 && entry > firstGoodEntry && !veto.GetError(1))
		EventErrors[23] = true;

	if (abs(veto.GetQEC2() - prev.GetQEC2()) > entry-prevGoodEntry
			&& entry > firstGoodEntry && veto.GetQEC2() != 0)
		EventErrors[24] = true;

	// Check errors 18-24 (don't accidentally change the run-level errors 25-28)
	for (int q=18; q<25; q++) {
		if (EventErrors[q]==1 && (q==18||q==19||q==20||q==21||q==22||q==23||q==24))
			skip = true;
	}
	// assign EventErrors to the output vetor
	ErrorVec = EventErrors;

	return skip;
}
