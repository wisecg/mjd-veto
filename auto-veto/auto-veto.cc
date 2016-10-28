// auto-veto.cc
// Runs during production, creates ROOT files of MJD veto data.
// C. Wiseman, A. Lopez, 10/24/16.

#include <iostream>
#include <fstream>
#include "TTreeReader.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"
#include "TFile.h"
#include "MJVetoEvent.hh"
#include "GATDataSet.hh"

using namespace std;

bool vetoFileCheck(GATDataSet *ds);
vector<int> vetoThreshFinder(GATDataSet *ds, bool makeQDCPlot=false);
void processVetoData(GATDataSet *ds, vector<int> thresholds, bool errorCheckOnly=false);

int findPanelThreshold(TH1F *qdcHist);
bool CheckEventErrors(MJVetoEvent veto, MJVetoEvent prev, MJVetoEvent first, long prevGoodEntry, int errorCode, int *ErrorArray=NULL);
double InterpTime(int entry, vector<double> times, vector<double> entries, vector<bool> badScaler);
int PanelMap(int i);

int main(int argc, char** argv)
{
	// NOTE: running over multiple runs is not recommended at this point due to lack of testing.
	if (argc < 1) {
		cout << "Usage: ./auto-veto [run number] [optional: upper run number]\n";
		return 0;
	}

	// TODO: clean up the output!
	// TODO: remove unused variables!

	// TODO: Change this to a TChain object so that it will work with veto-only files as well.
	// Create a GATDataSet object
	int run = stoi(argv[1]);
	int hiRun=0;
	if (argc > 2) {
		hiRun = stoi(argv[2]);
		printf("Processing runs %i through %i ...\n",run,hiRun);
	}
	else printf("Processing run %i ... \n",run);
	GATDataSet ds;
	if (hiRun == 0) ds.AddRunNumber(run);
	else if (hiRun != 0) ds.AddRunRange(run,hiRun);

	// TODO: options from vetoCheck - need to merge them in
	// int run = atoi(argv[1]);
	// if (run > 60000000){
	// 	cout << "Veto data not present in Module 2 runs.\n";
	// 	return 1;
	// }
	// bool draw = false;
	// string opt = "";
	// if (argc > 2) opt = argv[2];
	// if (argc > 2 && opt == "-d")
	// 	draw = true;
	// vetoCheck(run,draw);

	// Verify that the veto data exists in the given run(s)
	bool goodFile = vetoFileCheck(&ds);
	if (!goodFile) {
		cout << "Found corrupted files!  Exiting ...\n";
		return 0;
	}

	// Find the QDC pedestal location in each channel.
	// Give a threshold that is 35 QDC above this location,
	// and optionally output a plot that confirms this choice.
	gStyle->SetOptStat(0);
	vector<int> thresholds = vetoThreshFinder(&ds,false);	// true-makes qdc plots (should make multiplicity too)

	// Check for errors (using a special tag - loop 1),
	// tag muon and LED events in veto data,
	// and output a ROOT file for further analysis.
	bool errorCheckOnly = false;
	processVetoData(&ds,thresholds,errorCheckOnly);
	printf("Done processing.\n");
}

bool vetoFileCheck(GATDataSet *ds)
{
	double duration = ds->GetRunTime();
	TChain *v = ds->GetVetoChain();
	long vEntries = v->GetEntries();
	if (duration == 0){
		cout << "File not found.\n";
		return 0;
	}
	else if (vEntries == 0){
		cout << "No veto entries found.\n";
		return 0;
	}
	return 1;
}

vector<int> vetoThreshFinder(GATDataSet *ds, bool makeQDCPlot)
{
	// format: (panel 1) (threshold 1) (panel 2) (threshold 2) ...
	vector<int> thresholds;

	int loRun = ds->GetRunNumber(0), hiRun = 0;
	if (ds->GetNRuns() > 1) hiRun = ds->GetRunNumber(ds->GetNRuns()-1);

	int bins=500, lower=0, upper=500;
	TH1F *hLowQDC[32];
	TH1F *hFullQDC[32];
	char hname[50];
	for (int i = 0; i < 32; i++) {
		sprintf(hname,"hLowQDC%d",i);
		hLowQDC[i] = new TH1F(hname,hname,bins,lower,upper);
		sprintf(hname,"hFullQDC%d",i);
		hFullQDC[i] = new TH1F(hname,hname,420,0,4200);
	}

	TChain *v = ds->GetVetoChain();
	long vEntries = v->GetEntries();
	TTreeReader reader(v);
	TTreeReaderValue<unsigned int> vMult(reader, "mVeto");
	TTreeReaderValue<uint32_t> vBits(reader, "vetoBits");
	TTreeReaderValue<MGTBasicEvent> vEvt(reader,"vetoEvent");
	TTreeReaderValue<MJTRun> vRun(reader,"run");

	// Set all thresholds to 1, causing all entries to have a multiplicity of 32
	int def[32] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
	int errorCode = 0;
	long skippedEvents = 0, prevGoodEntry=0;
	MJVetoEvent veto, prev, first;

	while (reader.Next())
	{
		long i = reader.GetCurrentEntry();
		int run = vRun->GetRunNumber();
		MJVetoEvent veto;
		veto.SetSWThresh(def);
    	errorCode = veto.WriteEvent(i,&*vRun,&*vEvt,*vBits,run,true);
		if (CheckEventErrors(veto,prev,first,prevGoodEntry,errorCode)) {
    		skippedEvents++;
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
	if (skippedEvents > 0) printf("vetoThreshFinder skipped %li of %li entries.\n",skippedEvents,vEntries);

	int thresh[32] = {9999};
	for (int i = 0; i < 32; i++) {
		thresh[i] = findPanelThreshold(hLowQDC[i]);
		thresholds.push_back(i);
		thresholds.push_back(thresh[i]);
	}
	if (makeQDCPlot)
	{
		TCanvas *c1 = new TCanvas("c1","full QDC",0,0,1600,1200);
		c1->Divide(8,4,0,0);
		for (int i=0; i<32; i++)
		{
			c1->cd(i+1);
			TVirtualPad *vpad1 = c1->cd(i+1); vpad1->SetLogy();
			hFullQDC[i]->Draw();
		}
		TCanvas *c2 = new TCanvas("c2","QDC thresholds",0,0,1600,1200);
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
		char plotName[200], runRange[200];
		sprintf(runRange,"%i",loRun);
		if (hiRun != 0)
			sprintf(runRange,"%i-%i",loRun,hiRun);
		sprintf(plotName,"./output/qdc-%s-Full.pdf",runRange);
		c1->Print(plotName);
		sprintf(plotName,"./output/qdc-%s-Thresh.pdf",runRange);
		c2->Print(plotName);
	}
	return thresholds;
}

void processVetoData(GATDataSet *ds, vector<int> thresholds, bool errorCheckOnly)
{
	// TODO: organize all variables according to scope and use.
	// choose which ones go into the ROOT output file and which ones are just "internal"

	// Error checks
	const int nErrs = 29; // error 0 is unused
	int SeriousErrorCount = 0;
	int TotalErrorCount = 0;
	vector<double> EntryTime;
	vector<double> EntryNum;
	vector<bool> BadScalers;
	int EventNumPrev_good = 0;

	// Specify which error types to print during the loop over events
	vector<int> SeriousErrors = {1, 13, 14, 18, 19, 20, 21, 22, 23, 24};

	int EventError[nErrs] = {0};	// run-level error bools (true if any errors are found)
	int ErrorCount[nErrs] = {0};

	// LED variables
	int LEDMultipThreshold = 5;  // "multipThreshold" = "highestMultip" - "LEDMultipThreshold"
	int LEDSimpleThreshold = 15;  // used when LED frequency measurement is bad.

	int highestMultip=0, multipThreshold=0;
	double LEDWindow = 0.1;
	double LEDfreq=0, LEDperiod=0, LEDrms=0;
	bool badLEDFreq=false;
	TH1F *LEDDeltaT = new TH1F("LEDDeltaT","LEDDeltaT",100000,0,100); // 0.001 sec/bin
	bool IsLED = false, IsLEDPrev = false;
	bool LEDTurnedOff = false;

	// muon cut variables
	bool TimeCut = true;
	bool EnergyCut = false;

	// time variables
	double SBCOffset=0;
	bool ApproxTime = false;
	bool foundScalerJump = false;	// this is related to ApproxTime but should be kept separate.
	double deltaTadj = 0;
	double timeSBC=0, prevtimeSBC=0;
	int scalerJumps=0, skippedEvents=0, corruptScaler=0;

	double timeOffset=0, prevGoodTime=0, firstGoodScaler=0;
	int prevGoodEntry=0, simpleLEDCount=0;
	bool foundFirst=false, foundFirstScaler=false;
	MJVetoEvent first, prev, last;

	// custom SW Threshold (obtained from vetoThreshFinder)
	int swThresh[32] = {0};
	for (int i = 0; i < (int)thresholds.size(); i+=2)
		swThresh[thresholds[i]] = thresholds[i+1];

	// root output variables
	int errorCode=0, run=0, PlaneHitCount=0;
	long start=0, stop=0;
	double duration=0, livetime=0, xTime=0, x_deltaT=0, x_LEDDeltaT=0;
	bool useSimpleThreshold=false;
	vector<int> CoinType(32), CutType(32), PlaneHits(32), PlaneTrue(32);
	MJVetoEvent out;

	// used in loop 2
	MJVetoEvent prevLED;
	bool firstLED=false;
	int almostMissedLED=0;
	double xTimePrev=0, x_deltaTPrev=0, xTimePrevLED=0, xTimePrevLEDSimple=0, TSdifference=0;

	// initialize output file
	int loRun = ds->GetRunNumber(0), hiRun = 0;
	if (ds->GetNRuns() > 1) hiRun = ds->GetRunNumber(ds->GetNRuns()-1);

	char outputFile[200], runRange[200];
	sprintf(runRange,"%i.root",loRun);
	if (hiRun != 0) sprintf(runRange,"%i-%i.root",loRun,hiRun);
	sprintf(outputFile,"./output/veto_run%s",runRange);
	TFile *RootFile = new TFile(outputFile, "RECREATE");
	TTree *vetoTree = new TTree("vetoTree","MJD Veto Events");
	vetoTree->Branch("events","MJVetoEvent",&out,32000,1);	// TODO: organize these better
	vetoTree->Branch("errorCode",&errorCode,"errorCode/I");
	vetoTree->Branch("timeSBC",&timeSBC);
	vetoTree->Branch("LEDfreq",&LEDfreq);
	vetoTree->Branch("LEDrms",&LEDrms);
	vetoTree->Branch("multipThreshold",&multipThreshold);
	vetoTree->Branch("highestMultip",&highestMultip);
	vetoTree->Branch("LEDWindow",&LEDWindow);
	vetoTree->Branch("LEDMultipThreshold",&LEDMultipThreshold);
	vetoTree->Branch("LEDSimpleThreshold",&LEDSimpleThreshold); // TODO: add a bool signifiying this is in use
	vetoTree->Branch("useSimpleThreshold",&useSimpleThreshold);
	vetoTree->Branch("start",&start,"start/L");
	vetoTree->Branch("stop",&stop,"stop/L");
	vetoTree->Branch("duration",&duration);
	vetoTree->Branch("xTime",&xTime);
	vetoTree->Branch("x_deltaT",&x_deltaT);
	vetoTree->Branch("x_LEDDeltaT",&x_LEDDeltaT);
	vetoTree->Branch("CoinType",&CoinType);
	vetoTree->Branch("CutType",&CutType);
	vetoTree->Branch("PlaneHits",&PlaneHits);
	vetoTree->Branch("PlaneTrue",&PlaneTrue);
	vetoTree->Branch("PlaneHitCount",&PlaneHitCount);
	vetoTree->Branch("livetime",&livetime);

	// initialize input data
	TChain *v = ds->GetVetoChain();
	long vEntries = v->GetEntries();
	TTreeReader reader(v);
	TTreeReaderValue<unsigned int> vMult(reader, "mVeto");
	TTreeReaderValue<uint32_t> vBits(reader, "vetoBits");
	TTreeReaderValue<MGTBasicEvent> vEvt(reader,"vetoEvent");
	TTreeReaderValue<MJTRun> vRun(reader,"run");

	reader.Next();
	start = vRun->GetStartTime();
	stop = vRun->GetStopTime();
	duration = (double)(stop - start);
	reader.SetTree(v);  // resets the reader

	printf("\n======= Scanning run %i, %li entries, %.0f sec. =======\n",loRun,vEntries,duration);

	// ========= 1st loop over veto entries - Measure LED frequency. =========
	//
	// Goal is to measure the LED frequency, to be used in the second loop as a
	// time cut. This is done by finding the maximum bin of a delta-t histogram.
	// This section of the code uses a weak multiplicity threshold of 20 -- it
	// doesn't need to be exact, and should also work for runs where there were
	// only 24 panels installed.

	while(reader.Next())
	{
		long i = reader.GetCurrentEntry();
		int run = vRun->GetRunNumber();

		MJVetoEvent veto;
		veto.SetSWThresh(swThresh);
		errorCode = veto.WriteEvent(i,&*vRun,&*vEvt,*vBits,run,true);

		if (!veto.GetBadScaler()) {
			BadScalers.push_back(0);
			xTime = veto.GetTimeSec();
		}
		else {
			BadScalers.push_back(1);
			xTime = ((double)i / vEntries) * duration; // this breaks if we have corrupted duration
		}
		EntryNum.push_back(i);
		EntryTime.push_back(xTime);

		if (foundFirst && veto.GetError(1)) {
			foundFirst = false;
		}
		if (!foundFirstScaler && !veto.GetError(4)){
			foundFirstScaler = true;
			firstGoodScaler = veto.GetTimeSec();
		}
		if (CheckEventErrors(veto,prev,first,prevGoodEntry,errorCode,EventError)){
			skippedEvents++;
			continue;
		}
		if (!foundFirst && veto.GetTimeSBC() > 0 && veto.GetTimeSec() > 0 && !veto.GetError(4)) {
			first = veto;
			foundFirst = true;
		}
		if (veto.GetMultip() > highestMultip) {
			highestMultip = veto.GetMultip();
		}
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
	if (duration == 0) {
		cout << "Corrupted duration.  Last good timestamp: " << prevGoodTime-firstGoodScaler << endl;
		duration = prevGoodTime-firstGoodScaler;
	}
	livetime = duration - (first.GetTimeSec() - firstGoodScaler);

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
		cout << "Warning! No multiplicity > 15 events.  LED may be off.\n";
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
		}
		else {
			LEDperiod = 9999;
			badLEDFreq = true;
		}
	}
	if (LEDperiod > 9 || LEDperiod < 5 || badLEDFreq) {
		ErrorCount[25]++;
		EventError[25] = true;
	}

	// ========= 2nd loop over entries - Error checks =========

	reader.SetTree(v); // reset the reader
	while(reader.Next())
	{
		long i = reader.GetCurrentEntry();
		int run = vRun->GetRunNumber();

		// this time we don't skip anything.
		MJVetoEvent veto;
		veto.SetSWThresh(swThresh);
		errorCode = veto.WriteEvent(i,&*vRun,&*vEvt,*vBits,run,true);
		CheckEventErrors(veto,prev,first,prevGoodEntry,errorCode,EventError);
		for (int j=0; j<nErrs; j++) if (EventError[j]==1) ErrorCount[j]++;

		// find event time
		if (!veto.GetBadScaler())
		{
			xTime = veto.GetTimeSec();
			if(run > 8557 && veto.GetTimeSBC() < 2000000000)
				timeSBC = veto.GetTimeSBC() - SBCOffset;
		}
		else if (run > 8557 && veto.GetTimeSBC() < 2000000000)
			xTime = veto.GetTimeSBC() - SBCOffset;
		else
			xTime = InterpTime(i,EntryTime,EntryNum,BadScalers);

		EntryTime[i] = xTime;	// replace entry with the more accurate one

		// Print errors to screen
		bool PrintError = false;
		for (auto i : SeriousErrors){
			if (EventError[i]) {
				PrintError = true;
				break;
			}
		}
		// 1, 13, 14, 18, 19, 20, 21, 22, 23, 24 - must match vector "SeriousErrors"
		if (PrintError)
		{
			cout << "\nSerious errors found in entry " << i << ":\n";
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
					 << "  Prev TSdifference " << TSdifference
					 << "  Scaler DeltaT " << veto.GetTimeSec()-prev.GetTimeSec()
					 << "\n    Scaler Index " << veto.GetScalerIndex()
					 << "  Previous Scaler Index " << prev.GetScalerIndex()
					 << "  Scaler Time " << veto.GetTimeSec()
					 << "  SBC Time " << timeSBC << "\n";
				TSdifference = veto.GetTimeSec() - timeSBC;
			}
			if (EventError[19]) {
				cout << "EventError[19] Scaler Event Count Reset. "
					 << "  Scaler Index " << veto.GetScalerIndex()
					 << "  SEC " << veto.GetSEC()
					 << "  Previous SEC " << prev.GetSEC() << "\n";
			}
			if (EventError[20]) {
				cout << "EventError[20] Scaler Event Count Jump."
					 << "    xTime " << xTime
					 << "  Scaler Index " << veto.GetScalerIndex()
					 << "  SEC " << veto.GetSEC()
					 << "  Previous SEC " << prev.GetSEC() << "\n";
			}
			if (EventError[21]) {
				cout << "EventError[21] QDC1 Event Count Reset."
					 << "  Scaler Index " << veto.GetScalerIndex()
					 << "  QEC1 " << veto.GetQEC()
					 << "  Previous QEC1 " << prev.GetQEC() << "\n";
			}
			if(EventError[22]) {
				cout << "EventError[22] QDC 1 Event Count Jump."
					 << "  xTime " << xTime
					 << "  QDC 1 Index " << veto.GetQDC1Index()
					 << "  QEC 1 " << veto.GetQEC()
					 << "  Previous QEC 1 " << prev.GetQEC() << "\n";
			}
			if (EventError[23]) {
				cout << "EventError[23] QDC2 Event Count Reset."
					 << "  Scaler Index " << veto.GetScalerIndex()
					 << "  QEC2 " << veto.GetQEC2()
					 << "  Previous QEC2 " << prev.GetQEC2() << "\n";
			}
			if(EventError[24]) {
				cout << "EventError[24] QDC 2 Event Count Jump."
					 << "  xTime " << xTime
					 << "  QDC 2 Index " << veto.GetQDC2Index()
					 << "  QEC 2 " << veto.GetQEC2()
					 << "  Previous QEC 2 " << prev.GetQEC2() << "\n";
			}
		}

		TSdifference = veto.GetTimeSec() - timeSBC;
		prevtimeSBC = timeSBC;
		timeSBC = 0;
		prev = veto;
		last = prev;
		prevGoodEntry = i;

		// Reset error bools each entry
		for (int j=0; j<nErrs; j++) EventError[j]=false;
	}
	cout << "=================== End error scan. ======================\n";

	// ===== Summary: Calculate total errors and total serious errors ======

	for (int i = 1; i < nErrs; i++)
	{
			// Ignore 10 and 11, these will always be present
			// as long as the veto counters are not reset at the beginning of runs.
			if (i != 10 && i != 11) TotalErrorCount += ErrorCount[i];

			// make sure LED being off is counted as serious
			if (i == 25 && ErrorCount[i] > 0) SeriousErrorCount += ErrorCount[i];

			// count up the serious errors in the vector
			for (auto j : SeriousErrors)
				if (i == j)
					SeriousErrorCount += ErrorCount[i];
	}
	cout << "Serious errors found :: " << SeriousErrorCount << endl;
	if (SeriousErrorCount > 0)
	{
		cout << "Total Errors : " << TotalErrorCount << endl;
		if (duration != livetime)
			cout << "Run duration (" << duration << " sec) doesn't match live time: " << livetime << endl;

		for (int i = 1; i < nErrs; i++)
		{
			if (ErrorCount[i] > 0)
			{
				if (i != 25)
					cout << "  Error[" << i <<"]: " << ErrorCount[i] << " events ("
						 << 100*(double)ErrorCount[i]/vEntries << " %)\n";
		 		else if (i == 25)
				{
					cout << "  EventError[25]: Bad LED rate: " << LEDfreq << "  Period: " << LEDperiod << endl;
					if (LEDperiod > 0.1 && (abs(duration/LEDperiod) - simpleLEDCount) > 5)
					{
						cout << "   Simple LED count: " << simpleLEDCount
							 << "  Expected: " << (int)duration/(int)LEDperiod << endl;
					}
				}
			}
		}
		cout << "\n  For reference, \"serious\" error types are: ";
		for (auto i : SeriousErrors) cout << i << " ";
		cout << "\n  Please report these to the veto group.\n";
		cout << "================= End veto error report. =================\n";
	}
	if (errorCheckOnly) return;

	// ========= 3rd loop over veto entries - Find muons! Write ROOT output! =========

	cout << "================= Scanning for muons ... =================\n";
	reader.SetTree(v); // reset the reader
	prev.Clear();
	skippedEvents = 0;
	cout << "multiplicity threshold: " << multipThreshold << endl;

	while(reader.Next())
	{
		long i = reader.GetCurrentEntry();
		int run = vRun->GetRunNumber();

		MJVetoEvent veto;
		veto.SetSWThresh(swThresh);
		errorCode = veto.WriteEvent(i,&*vRun,&*vEvt,*vBits,run,true);
		CheckEventErrors(veto,prev,first,prevGoodEntry,errorCode,EventError);
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
			xTime = InterpTime(i,EntryTime,EntryNum,BadScalers);
			ApproxTime = true;
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
		// These will NOT be present in the final ROOT output, obviously.
		// TODO: implement some kind of "garbage collector",
		// or a variable signifying that this event is missing in the ROOT output.
		if (CheckEventErrors(veto,prev,first,prevGoodEntry,errorCode,EventError))
		{
			skippedEvents++;
			// do the end-of-run reset
			TSdifference = veto.GetTimeSec() - timeSBC;
			prevtimeSBC = timeSBC;
			timeSBC = 0;
			prev = veto;
			last = prev;
			prevGoodEntry = i;
			continue;
		}

		// LED / Time Cut

		// IsLED = true;
		TimeCut = false;
		LEDTurnedOff = ErrorCount[25];
		x_deltaT = xTime - xTimePrevLED;

		// TODO: right now this enforces a hard multiplicity cut.
		// Need to evaluate if the frequency-based multiplicity cut will
		// recover high-multiplicity muon events without mis-identifying LED's as muons.
		if (veto.GetMultip() < multipThreshold && !LEDTurnedOff) {
			// IsLED = false;
			TimeCut = true;
		}
		else if (LEDTurnedOff)
			TimeCut = true;

		/*
		// if (!LEDTurnedOff && !badLEDFreq && fabs(LEDperiod - x_deltaT) < LEDWindow && veto.GetMultip() > multipThreshold)
		// {
		// 	TimeCut = false;
		// 	IsLED = true;
		// }
		// // almost missed a high-multiplicity event somehow ...
		// // often due to skipping previous events.
		// else if (!LEDTurnedOff && !badLEDFreq && fabs(LEDperiod - x_deltaT) >= (LEDperiod - LEDWindow)
		// 				&& veto.GetMultip() > multipThreshold && i > 1)
		// {
		// 	TimeCut = false;
		// 	IsLED = true;
		// 	almostMissedLED++;
		// 	cout << "Almost missed LED:\n";
		// 	printf("Current: %-3li  m %-3i LED? %i t %-6.2f LEDP %-5.2f  XDT %-6.2f LEDP-XDT %-6.2f\n"
		// 		,i,veto.GetMultip(),IsLED,xTime,LEDperiod,x_deltaT,LEDperiod-x_deltaT);
		// }
		// else TimeCut = true;
		// // Grab first LED
		// if (!LEDTurnedOff && !firstLED && veto.GetMultip() > multipThreshold) {
		// 	printf("Found first LED.  i %-2li m %-2i t %-5.2f thresh:%i  LEDoff:%i\n\n",i,veto.GetMultip(),xTime,multipThreshold,LEDTurnedOff);
		// 	IsLED=true;
		// 	firstLED=true;
		// 	TimeCut=false;
		// 	x_deltaT = -1;
		// }
		// // If frequency measurement is bad, revert to standard multiplicity cut
		// if (badLEDFreq && veto.GetMultip() >= LEDSimpleThreshold){
		// 	IsLED = true;
		// 	TimeCut = false;
		// }
		// // Simple x_LEDDeltaT uses the multiplicity-only threshold, veto.GetMultip() > multipThreshold.
		// x_LEDDeltaT = xTime - xTimePrevLEDSimple;
		//
		// // If LED is off, all events pass time cut.
		// if (LEDTurnedOff) {
		// 	IsLED = false;
		// 	TimeCut = true;
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
		// printf("Entry %li  Time %-6.2f  QDC %-5i  Mult %i  hits>500 %i  Loff? %i  T/LCut %i  ECut %i  \n",i,xTime,veto.GetTotE(),veto.GetMultip(),over500Count,LEDTurnedOff,TimeCut,EnergyCut);


		// Hit Pattern "Cut": Map hits above SW threshold to planes and count the hits.
		PlaneHitCount = 0;
		for (int k = 0; k < 12; k++) {
			PlaneTrue[k] = 0;
			PlaneHits[k]=0;
		}
		for (int k = 0; k < 32; k++) {
			if (veto.GetQDC(k) > veto.GetSWThresh(k)) {
				if (PanelMap(k)==0) { PlaneTrue[0]=1; PlaneHits[0]++; }			// 0: Lower Bottom
				else if (PanelMap(k)==1) { PlaneTrue[1]=1; PlaneHits[1]++; }		// 1: Upper Bottom
				else if (PanelMap(k)==2) { PlaneTrue[2]=1; PlaneHits[2]++; }		// 3: Inner Top
				else if (PanelMap(k)==3) { PlaneTrue[3]=1; PlaneHits[3]++; }		// 4: Outer Top
				else if (PanelMap(k)==4) { PlaneTrue[4]=1; PlaneHits[4]++; }		// 5: Inner North
				else if (PanelMap(k)==5) { PlaneTrue[5]=1; PlaneHits[5]++; }		// 6: Outer North
				else if (PanelMap(k)==6) { PlaneTrue[6]=1; PlaneHits[6]++; }		// 7: Inner South
				else if (PanelMap(k)==7) { PlaneTrue[7]=1; PlaneHits[7]++; }		// 8: Outer South
				else if (PanelMap(k)==8) { PlaneTrue[8]=1; PlaneHits[8]++; }		// 9: Inner West
				else if (PanelMap(k)==9) { PlaneTrue[9]=1; PlaneHits[9]++; }		// 10: Outer West
				else if (PanelMap(k)==10) { PlaneTrue[10]=1; PlaneHits[10]++; }	// 11: Inner East
				else if (PanelMap(k)==11) { PlaneTrue[11]=1; PlaneHits[11]++; }	// 12: Outer East
			}
		}
		for (int k = 0; k < 12; k++) if (PlaneTrue[k]) PlaneHitCount++;

		// Muon Identification
		// Use EnergyCut, TimeCut, and the Hit Pattern to identify them sumbitches.
		for (int r = 0; r < 32; r++) {CoinType[r]=0; CutType[r]=0;}
		if (TimeCut && EnergyCut)
		{
			int type = 0;
			bool a=0,b=0,c=0;
			CoinType[0] = true;

			if (PlaneTrue[0] && PlaneTrue[1] && PlaneTrue[2] && PlaneTrue[3]) {
				CoinType[1]==true;
				a=true;
				type=1;
			}
			if ((PlaneTrue[0] && PlaneTrue[1]) && ((PlaneTrue[2] && PlaneTrue[3]) || (PlaneTrue[4] && PlaneTrue[5])
					|| (PlaneTrue[6] && PlaneTrue[7]) || (PlaneTrue[8] && PlaneTrue[9]) || (PlaneTrue[10] && PlaneTrue[11]))){
				CoinType[2] = true;
				b=true;
				type = 2;
			}
			if ((PlaneTrue[2] && PlaneTrue[3]) && ((PlaneTrue[4] && PlaneTrue[5]) || (PlaneTrue[6] && PlaneTrue[7])
					|| (PlaneTrue[8] && PlaneTrue[9]) || (PlaneTrue[10] && PlaneTrue[11]))) {
				CoinType[3] = true;
				c=true;
				type = 3;
			}
			if ((a && b)||(a && c)||(b && c)) type = 4;

			char hitType[200];
			if (type==0) sprintf(hitType,"2+ panels");
			if (type==1) sprintf(hitType,"vertical");
			if (type==2) sprintf(hitType,"side+bottom");
			if (type==3) sprintf(hitType,"top+sides");
			if (type==4) sprintf(hitType,"compound");
			printf("Hit: %-12s Entry %-4li  Time %-6.2f  QDC %-5i  Mult %i  T/LCut %i  ECut %i  ApxT %i\n", hitType,i,xTime,veto.GetTotE(),veto.GetMultip(),TimeCut,EnergyCut,ApproxTime);
		}

		// Output

		CutType[0] = LEDTurnedOff;
		CutType[1] = EnergyCut;
		CutType[2] = ApproxTime;
		CutType[3] = TimeCut;
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
	}

	printf("\n===================== End of Scan. =====================\n");
	if (skippedEvents > 0) cout << "processVetoData skipped " << skippedEvents << " events.\n";

	vetoTree->Write("",TObject::kOverwrite);
	RootFile->Close();
	cout << "Wrote ROOT file.\n";
}
// ====================================================================================
// ====================================================================================

bool CheckEventErrors(MJVetoEvent veto, MJVetoEvent prev, MJVetoEvent first, long prevGoodEntry, int errorCode, int *ErrorArray)
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
	*/

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

	int EventErrors[29] = {0};

	// Errors 1-18 are checked automatically when we call MJVetoEvent::WriteEvent
	veto.UnpackErrorCode(errorCode,EventErrors);
	for (int q=0; q<18; q++) {
		if (EventErrors[q]==1 && (q==1||q==2||q==3||q==5||q==9||q==13||q==14))
			skip = true;
	}

	if (veto.GetBadScaler() && (veto.GetRun() < 8557 || veto.GetTimeSBC() > 2000000000))
		EventErrors[25] = true;

	int entry = veto.GetEntry();
	int prevEntry = prev.GetEntry();
	int firstGoodEntry = first.GetEntry();

	if (ErrorArray!=NULL) memcpy(ErrorArray,EventErrors,sizeof(EventErrors));

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

	// Check errors 18-25
	for (int q=18; q<26; q++) {
		if (EventErrors[q]==1 && (q==18||q==19||q==20||q==21||q==22||q==23||q==24))
			skip = true;
	}

	if (ErrorArray!=NULL) memcpy(ErrorArray,EventErrors,sizeof(EventErrors));
	return skip;
}

int findPanelThreshold(TH1F *qdcHist)
{
	// Place threshold 35 qdc above pedestal location.
	int firstNonzeroBin = qdcHist->FindFirstBinAbove(1,1);
	qdcHist->GetXaxis()->SetRange(firstNonzeroBin-10,firstNonzeroBin+50);
	//qdcHist->GetXaxis()->SetRangeUser(0,500); //alternate method of finding pedestal
	int bin = qdcHist->GetMaximumBin();
	if (firstNonzeroBin == -1) return -1;
	double xval = qdcHist->GetXaxis()->GetBinCenter(bin);
	return xval+35;
}

double InterpTime(int entry, vector<double> times, vector<double> entries, vector<bool> badScaler)
{
	if ((times.size() != entries.size()) || times.size() != badScaler.size())
	{
		cout << "Vectors are different sizes!\n";
		if (entry >= (int)times.size())
			cout << "Entry is larger than number of entries in vector!\n";
		return -1;
	}

	double iTime = 0;
	double lower = 0;
	double upper = 0;
	if (!badScaler[entry]) iTime = times[entry];
	else
	{
		for (int i = entry; i < (int)entries.size(); i++)
			if (badScaler[i] == 0) { upper = times[i]; break; }
		for (int j = entry; j > 0; j--)
			if (badScaler[j] == 0) { lower = times[j]; break; }
		iTime = (upper + lower)/2.0;
	}
	return iTime;
}

int PanelMap(int i)
{
	// For tagging plane-based coincidences.
	// This uses a zero-indexed map.

	// 0: Lower Bottom
	// 1: Upper Bottom
	// 2: Inner Top
	// 3: Outer Top
	// 4: Inner North
	// 5: Outer North
	// 6: Inner South
	// 7: Outer South
	// 8: Inner West
	// 9: Outer West
	// 10: Inner East
	// 11: Outer East

	if 		(i == 0) return 0;  // L-bot 1
	else if (i == 1) return 0;
	else if (i == 2) return 0;
	else if (i == 3) return 0;
	else if (i == 4) return 0;
	else if (i == 5) return 0;  // L-bot 6

	else if (i == 6) return 1;  // U-bot 1
	else if (i == 7) return 1;
	else if (i == 8) return 1;
	else if (i == 9) return 1;
	else if (i == 10) return 1;
	else if (i == 11) return 1; // U-bot 6

	else if (i == 17) return 3; // Top outer
	else if (i == 18) return 3; // Top outer
	else if (i == 20) return 2; // Top inner
	else if (i == 21) return 2; // Top inner

	else if (i == 15) return 5; // North outer
	else if (i == 16) return 5; // North outer
	else if (i == 19) return 4; // North inner
	else if (i == 23) return 4; // North inner

	else if (i == 24) return 6; // South inner
	else if (i == 25) return 7; // South outer
	else if (i == 26) return 6; // South inner
	else if (i == 27) return 7; // South outer

	else if (i == 12) return 8; // West inner
	else if (i == 13) return 8; // West inner
	else if (i == 14) return 9; // West outer
	else if (i == 22) return 9; // West outer

	else if (i == 28) return 10; // East inner
	else if (i == 29) return 11; // East outer
	else if (i == 30) return 10; // East inner
	else if (i == 31) return 11; // East outer

	else return -1;
}