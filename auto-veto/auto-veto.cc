// auto-veto.cc
// Runs during production, creates ROOT files of MJD veto data.
// C. Wiseman, 10/24/16.

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
void processVetoData(GATDataSet *ds, vector<int> thresholds);

int findPanelThreshold(TH1F *qdcHist);
bool vEventErrorCheck(MJVetoEvent veto, int entry, int isGood, bool verbose);
bool CheckForAllErrors(MJVetoEvent veto, int entry, int isGood, MJVetoEvent prev, int EventNumPrev_good, MJVetoEvent first, int errornum, int firstgoodentry);
int PanelMap(int i);

int main(int argc, char** argv)
{
	// NOTE: running over multiple runs is not recommended at this point due to lack of testing.
	if (argc < 1) {
		cout << "Usage: ./auto-veto [run number] [optional: upper run number]\n";
		return 0;
	}
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

	// Verify that the veto data exists in the given run.
	bool goodFile = vetoFileCheck(&ds);
	if (!goodFile) {
		cout << "Found corrupted files!  Exiting ...\n";
		return 0;
	}
	// Find the QDC Pedestal location in each channel.
	// Give a threshold that is 20 QDC above this location, and optionally output a plot
	// that confirms this choice.
	gStyle->SetOptStat(0);
	vector<int> thresholds = vetoThreshFinder(&ds,true);

	// Tag muon and LED events in veto data,
	// and output a ROOT file for further analysis.
	processVetoData(&ds,thresholds);
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

	// standard veto initialization block
	TChain *v = ds->GetVetoChain();
	long vEntries = v->GetEntries();
	TTreeReader reader(v);
	TTreeReaderValue<unsigned int> vMult(reader, "mVeto");
	TTreeReaderValue<uint32_t> vBits(reader, "vetoBits");
	TTreeReaderValue<MGTBasicEvent> vEvt(reader,"vetoEvent");
	TTreeReaderValue<MJTRun> vRun(reader,"run");

	// Initially set all thresholds to 1, causing all entries to have a multiplicity of 32
	int def[32] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
	int isGood = 0;
	long skippedEvents = 0;
	while (reader.Next())
	{
		long i = reader.GetCurrentEntry();
		int run = vRun->GetRunNumber();
		MJVetoEvent veto;
		veto.SetSWThresh(def);
    	isGood = veto.WriteEvent(i,&*vRun,&*vEvt,*vBits,run,true);
    	if (vEventErrorCheck(veto,i,isGood,false)) {
    		skippedEvents++;
    		continue;
    	}
    	for (int q = 0; q < 32; q++) {
    		hLowQDC[q]->Fill(veto.GetQDC(q));
    		hFullQDC[q]->Fill(veto.GetQDC(q));
		}
	}
	if (skippedEvents > 0) printf("Skipped %li of %li entries.\n",skippedEvents,vEntries);

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

void processVetoData(GATDataSet *ds, vector<int> thresholds)
{
	// LED Cut Parameters (C-f "Display Cut Parameters" below.)
	double LEDWindow = 0.1;
	int LEDMultipThreshold = 10;  // "multipThreshold" = "highestMultip" - "LEDMultipThreshold"
	int LEDSimpleThreshold = 20;  // used when LED frequency measurement is bad.

	// Custom SW Threshold (obtained from vetoThreshFinder)
	int swThresh[32] = {0};
	for (int i = 0; i < (int)thresholds.size(); i+=2)
		swThresh[thresholds[i]] = thresholds[i+1];

	// Set up output file
	int loRun = ds->GetRunNumber(0), hiRun = 0;
	if (ds->GetNRuns() > 1) hiRun = ds->GetRunNumber(ds->GetNRuns()-1);

	char outputFile[200], runRange[200];
	sprintf(runRange,"%i.root",loRun);
	if (hiRun != 0) sprintf(runRange,"%i-%i.root",loRun,hiRun);
	sprintf(outputFile,"./output/veto-run%s",runRange);

	TFile *RootFile = new TFile(outputFile, "RECREATE");
  	TH1::AddDirectory(kFALSE);
	TTree *vetoTree = new TTree("vetoTree","MJD Veto Events");
	int isGood=0, run=0, PlaneHitCount=0, highestMultip=0, multipThreshold=0;
	long rEntry=0, start=0, stop=0;
	double duration=0, LEDfreq=0, LEDrms=0, xTime=0, x_deltaT=0, x_LEDDeltaT=0, timeSBC=0;
	vector<int> CoinType(32), CutType(32), PlaneHits(32), PlaneTrue(32);
	MJVetoEvent out;
	vetoTree->Branch("events","MJVetoEvent",&out,32000,1);
	vetoTree->Branch("rEntry",&rEntry,"rEntry/L");
	vetoTree->Branch("timeSBC",&timeSBC);
	vetoTree->Branch("LEDfreq",&LEDfreq);
	vetoTree->Branch("LEDrms",&LEDrms);
	vetoTree->Branch("multipThreshold",&multipThreshold);
	vetoTree->Branch("highestMultip",&highestMultip);
	vetoTree->Branch("LEDWindow",&LEDWindow);
	vetoTree->Branch("LEDMultipThreshold",&LEDMultipThreshold);
	vetoTree->Branch("LEDSimpleThreshold",&LEDSimpleThreshold); // TODO: add a bool signifiying this is in use
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

	// standard veto initialization block
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
	reader.SetTree(v);  // reset the reader

	printf("\n======= Scanning run %i, %li entries, %.0f sec. =======\n",loRun,vEntries,duration);
	cout << "start: " << start << "  stop: " << stop << endl;

	// ========= 1st loop over veto entries - Measure LED frequency. =========
	//
	// Goal is to measure the LED frequency, to be used in the second loop as a
	// time cut. This is done by finding the maximum bin of a delta-t histogram.
	// This section of the code uses a weak multiplicity threshold of 20 -- it
	// doesn't need to be exact, and should also work for runs where there were
	// only 24 panels installed.

	highestMultip=0;
	double timeOffset=0;
	int scalerJumps=0, firstGoodEntry=0;
	long skippedEvents=0, corruptScaler=0;
	bool badLEDFreq=false, foundFirst=false;
	MJVetoEvent first, prev;

	char hname[200];
	sprintf(hname,"LEDDeltaT_run%i",run);
	TH1F *LEDDeltaT = new TH1F(hname,hname,100000,0,100); // 0.001 sec/bin

	while(reader.Next())
	{
		long i = reader.GetCurrentEntry();
		int run = vRun->GetRunNumber();

		MJVetoEvent veto;
		veto.SetSWThresh(swThresh);
    	isGood = veto.WriteEvent(i,&*vRun,&*vEvt,*vBits,run,true);
    	if (vEventErrorCheck(veto,i,isGood,false)) {
    		skippedEvents++;
    		continue;
    	}
    	if (veto.GetBadScaler())
			corruptScaler++;
    	if (veto.GetMultip() > highestMultip && veto.GetMultip() < 33)
			highestMultip = veto.GetMultip();

		if (isGood && !foundFirst && veto.GetTimeSBC()>0.01 && veto.GetTimeSec()>0.01 && !veto.GetBadScaler())
		{
			first = veto;		// save first good entry info
			foundFirst = true;
			firstGoodEntry = i;
			timeOffset = veto.GetTimeSec();
		}
		if (veto.GetMultip() >= 20)
			LEDDeltaT->Fill(veto.GetTimeSec()-prev.GetTimeSec()); // simple LED tag

		prev = veto;
	}
	if (skippedEvents > 0) printf("Skipped %li of %li entries.\n",skippedEvents,vEntries);
	if (corruptScaler > 0) printf("Corrupt scaler: %li of %li entries (%.2f%%) .\n"
		,corruptScaler,vEntries,100*(double)corruptScaler/vEntries);

	// Find the SBC offset
	double SBCOffset = first.GetTimeSBC() - first.GetTimeSec();
	printf("First good entry: %i  Scaler %.2f  SBC %.2f  SBCOffset %.2f\n"
		,firstGoodEntry,first.GetTimeSec(),first.GetTimeSBC(),SBCOffset);

	// Find the LED frequency
	bool LEDTurnedOff = false;
	if (highestMultip < 20) {
		printf("Warning!  LED's may be off!\n");
		LEDTurnedOff = true;
	}
	LEDrms = 0;
	LEDfreq = 0;
	int dtEntries = LEDDeltaT->GetEntries();
	if (dtEntries > 0) {
		int maxbin = LEDDeltaT->GetMaximumBin();
		LEDDeltaT->GetXaxis()->SetRange(maxbin-100,maxbin+100); // looks at +/- 0.1 seconds of max bin.
		LEDrms = LEDDeltaT->GetRMS();
		if (LEDrms==0) LEDrms = 0.1;
		LEDfreq = 1/LEDDeltaT->GetMean();
	}
	else {
		printf("Warning! No multiplicity > 20 events!!\n");
		LEDrms = 9999;
		LEDfreq = 9999;
		LEDTurnedOff = true;
	}
	double LEDperiod = 1/LEDfreq;

	// Display LED Cut parameters
	multipThreshold = highestMultip - LEDMultipThreshold;
	printf("HM: %i LED_f: %.8f LED_t: %.8f RMS: %8f\n",highestMultip,LEDfreq,1/LEDfreq,LEDrms);
	printf("LED window: %.2f  Multip Threshold: %i\n",LEDWindow,multipThreshold);
	if (LEDperiod > 9 || vEntries < 100) {
		badLEDFreq = true;
		printf("Warning: LED period is %.2f, total entries: %li.  Can't use it in the time cut!\n",LEDperiod,vEntries);
	}
	delete LEDDeltaT;

	// ========= 2nd loop over veto entries - Find muons! =========

	reader.SetTree(v); // reset the reader
	prev.Clear();
	MJVetoEvent prevLED;
	bool firstLED=false;
	int almostMissedLED=0;
	double xTimePrev=0, x_deltaTPrev=0, xTimePrevLED=0, xTimePrevLEDSimple=0, TSdifference=0;

	while(reader.Next())
	{
		long i = reader.GetCurrentEntry();
		int run = vRun->GetRunNumber();

		MJVetoEvent veto;
		veto.SetSWThresh(swThresh);
		isGood = veto.WriteEvent(i,&*vRun,&*vEvt,*vBits,run,true);
		rEntry = i;
		timeSBC = veto.GetTimeSBC()-SBCOffset;

		// 0: Time of event and skipping if necessary.
		// Employ alternate methods if the scaler is corrupted.
		// Should implement an estimate of the error when alternate methods are used.

		bool ApproxTime = false;
		xTime = -1;
		if (!veto.GetBadScaler())
		{
			xTime = veto.GetTimeSec() - timeOffset;

			// Find scaler jumps and adjust xTime by "TSdifference"
			// TSdifference starts at 0 at the beginning of the run.
			if (veto.GetTimeSec() != 0 && veto.GetTimeSBC() !=0 && SBCOffset != 0 && i>=firstGoodEntry)
			{
				double sbc = veto.GetTimeSBC() - SBCOffset;
				double diff = veto.GetTimeSec() - sbc;

				// if (fabs(fabs(diff) - TSdifference) > 1)	// andrew's original method (11472 - fails)
				if (fabs(diff-TSdifference) > 1)	// clint's method (11472 bkwds - OK)
				{
					scalerJumps++;
					TSdifference = diff;
					printf("i %li  Scaler Jump! Adjusting all following timestamps by: %.2f\n",i,diff);
					printf("   diff (scaler-sbc) %.2f  TSdiff %.2f  diff-TSdiff %.2f\n",diff,TSdifference,diff-TSdifference);
				}
			}

			// modify xTime by the running difference in timestamps
			xTime -= TSdifference;

    		// printf("i %i  scaler %.2f  sbc %.2f  xTime %.2f\n"
    			// ,i,veto.GetTimeSec(),veto.GetTimeSBC()-SBCOffset,xTime);
		}
		else if (run > 8557 && veto.GetTimeSBC() < 2000000000) {
			xTime = veto.GetTimeSBC() - SBCOffset;
			ApproxTime = true;
		}
    	else {
    		xTime = ((double)i / vEntries) * duration;
    		ApproxTime = true;
    	}

    	// Skip events after the event time is calculated.
    	if (vEventErrorCheck(veto,i,isGood,false))
    	{
    		printf("Skipping Entry %li.  Errors: ",i);

    		for (int j=0; j<18; j++) if (veto.GetError(j)==1)
    		{
    			cout << j << " ";
    		}
    		cout << endl;
    		// cout << "\n \t Full event summary: " << endl;
    		// veto.Print();

    		// do the end-of-run reset
    		// if (veto.GetMultip() > multipThreshold) {
				// xTimePrevLEDSimple = xTime;
			// }
			// IsLEDPrev = IsLED;
			// prev = veto;
			xTimePrev = xTime;
			x_deltaTPrev = x_deltaT;
    		continue;
    	}

		// 1. LED Cut
		// TRUE if an event PASSES (i.e. is physics.)  FALSE if an event is an LED.
		// If LED's are turned off or the frequency measurement is bad, we revert
		// to a simple multiplicity threshold.

		bool TimeCut = true;
		bool IsLED = false;

		x_deltaT = xTime - xTimePrevLED;
		if (!LEDTurnedOff && !badLEDFreq && fabs(LEDperiod - x_deltaT) < LEDWindow && veto.GetMultip() > multipThreshold)
		{
			TimeCut = false;
			IsLED = true;
		}

		// almost missed a high-multiplicity event somehow ...
		// often due to skipping previous events.
		else if (!LEDTurnedOff && !badLEDFreq && fabs(LEDperiod - x_deltaT) >= (LEDperiod - LEDWindow) && veto.GetMultip() > multipThreshold)
		{
			TimeCut = false;
			IsLED = true;
			almostMissedLED++;
			cout << "Almost missed LED:\n";

			// check this entry
			printf("Current: %-3li  m %-3i LED? %i t %-6.2f LEDP %-5.2f  XDT %-6.2f LEDP-XDT %-6.2f\n"
				,i,veto.GetMultip(),IsLED,xTime,LEDperiod,x_deltaT,LEDperiod-x_deltaT);

			// check previous entry
			// printf("Previous: %-3li  m %-3i LED? %i t %-6.2f LEDP %-5.2f  XDT %-6.2f LEDP-XDT %-6.2f LEDW %-6.2f\n"
				// ,i-1,prev.GetMultip(),IsLEDPrev,xTimePrev,LEDperiod,x_deltaTPrev,LEDperiod-x_deltaTPrev,LEDWindow);

			printf("Bools: IsLED %i  TimeCut %i  LEDTurnedOff %i  badLEDFreq %i\n"
				,IsLED,TimeCut,LEDTurnedOff,badLEDFreq);
		}
		else TimeCut = true;

		// Grab first LED
		if (!LEDTurnedOff && !firstLED && veto.GetMultip() > multipThreshold) {
			printf("Found first LED.  i %-2li m %-2i t %-5.2f\n\n",i,veto.GetMultip(),xTime);
			IsLED=true;
			firstLED=true;
			TimeCut=false;
			x_deltaT = -1;
		}

		// If frequency measurement is bad, revert to standard multiplicity cut
		if (badLEDFreq && veto.GetMultip() >= LEDSimpleThreshold){
			IsLED = true;
			TimeCut = false;
		}
		// Simple x_LEDDeltaT uses the multiplicity-only threshold, veto.GetMultip() > multipThreshold.
		x_LEDDeltaT = xTime - xTimePrevLEDSimple;

		// If LED is off, all events pass time cut.
		if (LEDTurnedOff) {
			IsLED = false;
			TimeCut = true;
		}
		// // Check output
		// printf("%-3li  m %-3i LED? %i t %-6.2f LEDP %-5.2f  XDT %-6.2f LEDP-XDT %-6.2f\n"
		// 	,i,veto.GetMultip(),IsLED,xTime,LEDperiod,x_deltaT,LEDperiod-x_deltaT);

    	// 2: Energy (Gamma) Cut
    	// The measured muon energy threshold is QDC = 500.
    	// Set TRUE if at least TWO panels are over 500.

    	bool EnergyCut = false;

    	int over500Count = 0;
    	for (int q = 0; q < 32; q++) {
    		if (veto.GetQDC(q) > 500)
    			over500Count++;
    	}
    	if (over500Count >= 2) EnergyCut = true;

		// 3: Hit Pattern
		// Map hits above SW threshold to planes and count the hits.

		// reset
		PlaneHitCount = 0;
		for (int k = 0; k < 12; k++) {
			PlaneTrue[k] = 0;
			PlaneHits[k]=0;
		}
		for (int k = 0; k < 32; k++)
		{
			if (veto.GetQDC(k) > veto.GetSWThresh(k))
			{
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
		for (int k = 0; k < 12; k++) {
			if (PlaneTrue[k]) PlaneHitCount++;
		}

		// 4: Muon Identification
		// Use EnergyCut, TimeCut, and the Hit Pattern to identify them sumbitches.

		// reset
		for (int r = 0; r < 32; r++) {CoinType[r]=0; CutType[r]=0;}

		// Check output
		// printf("%-3li  m %-3i  t %-6.2f  XDT %-6.2f  LED? %i  TC %i  EC %i  QTot %i\n"
			// ,i,veto.GetMultip(),xTime,x_deltaT,IsLED,TimeCut,EnergyCut,veto.GetTotE());

		if (TimeCut && EnergyCut)
		{
			// 0. Everything that passes TimeCut and EnergyCut.
			// This is what goes into the DEMONSTRATOR veto cut.
			CoinType[0] = true;
			printf("Entry: %li  2+Panel Muon.  QDC: %i  Mult: %i  LED? %i  T: %-6.2f  XDT %-6.2f  LEDP-XDT %-6.2f\n",
				i,veto.GetTotE(),veto.GetMultip(),IsLED,xTime,x_deltaT,LEDperiod-x_deltaT);

			// 1. Definite Vertical Muons
			if (PlaneTrue[0] && PlaneTrue[1] && PlaneTrue[2] && PlaneTrue[3]) {
				CoinType[1] = true;
				printf("Entry: %li  Vertical Muon.  QDC: %i  Mult: %i  LED? %i  T: %-6.2f  XDT %-6.2f  LEDP-XDT %-6.2f\n",
					i,veto.GetTotE(),veto.GetMultip(),IsLED,xTime,x_deltaT,LEDperiod-x_deltaT);
			}

			// 2. Both top or side layers + both bottom layers.
			if ((PlaneTrue[0] && PlaneTrue[1]) && ((PlaneTrue[2] && PlaneTrue[3]) || (PlaneTrue[4] && PlaneTrue[5])
				|| (PlaneTrue[6] && PlaneTrue[7]) || (PlaneTrue[8] && PlaneTrue[9]) || (PlaneTrue[10] && PlaneTrue[11]))) {
				CoinType[2] = true;

				// show output if we haven't seen it from CT1 already
				if (!CoinType[1]) {
					printf("Entry: %li  Side+Bottom Muon.  QDC: %i  Mult: %i  LED? %i  T: %-6.2f  XDT %-6.2f  LEDP-XDT %-6.2f\n",
						i,veto.GetTotE(),veto.GetMultip(),IsLED,xTime,x_deltaT,LEDperiod-x_deltaT);
				}
			}

			// 3. Both Top + Both Sides
			if ((PlaneTrue[2] && PlaneTrue[3]) && ((PlaneTrue[4] && PlaneTrue[5]) || (PlaneTrue[6] && PlaneTrue[7])
				|| (PlaneTrue[8] && PlaneTrue[9]) || (PlaneTrue[10] && PlaneTrue[11]))) {
				CoinType[3] = true;

				// show output if we haven't seen it from CT1 or CT2 already
				if (!CoinType[1] && !CoinType[2]) {
					printf("Entry: %li  Top+Sides Muon.  QDC: %i  Mult: %i  LED? %i  T: %-6.2f  XDT %-6.2f  LEDP-XDT %-6.2f\n",
						i,veto.GetTotE(),veto.GetMultip(),IsLED,xTime,x_deltaT,LEDperiod-x_deltaT);
				}
			}
		}
		// 5: Output
		// Write the ROOT file containing all the real data.

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
		if (IsLED) {
			prevLED = veto;
			xTimePrevLED = xTime;
		}
		if (veto.GetMultip() > multipThreshold) {
			xTimePrevLEDSimple = xTime;
		}
		// IsLEDPrev = IsLED;
		prev = veto;
		xTimePrev = xTime;
		x_deltaTPrev = x_deltaT;
	}

	printf("\n===================== End of Scan. =====================\n");
	if (almostMissedLED > 0) cout << "\nWarning, almost missed " << almostMissedLED << " LED events.\n";
	if (scalerJumps > 0) cout << "\nWarning, found " << scalerJumps << " scaler jumps.\n";

	vetoTree->Write("",TObject::kOverwrite);
	RootFile->Close();
}

// ====================================================================================
// ====================================================================================

// Check the 18 built-in error types in a MJVetoEvent object to see if veto entry is worth analyzing further
bool vEventErrorCheck(MJVetoEvent veto, int entry, int isGood, bool verbose)
{
	bool badError = false;

	if (isGood != 1) {
		int error[29] = {0};
		veto.UnpackErrorCode(isGood,error);

		// search for particular errors and set a flag.
		for (int q=0; q<28; q++)
		{
			// 4: don't skip bad-scaler events
			// 10: P3K93: don't skip "event count doesn't match ROOT entry" events
			// 7 & 11: don't skip hw count mismatches or scaler != qdc event count.  Due to problems with continuous running mode.

			// if (q!=4 && q!=10 && error[q] == 1) badError = true;
			if (q != 4 && q != 7 && q != 10 && q != 11 && q != 12 && error[q] == 1) badError = true;
		}


		if (badError) {
			if (verbose == true) {
				cout << "Skipped Entry: " << entry << endl;
				veto.Print();
				cout << endl;
			}
			return badError;
		}
	}
	else return badError;	// no bad errors found

	return false;
}

bool CheckForAllErrors(MJVetoEvent veto, int entry, int isGood, MJVetoEvent prev, int EventNumPrev_good, MJVetoEvent first, int errornum, int firstgoodentry) //for all entries
{
	//now finds all errors except 25, 26, 27
	//first good entry can't find errors 18, 20, 22, 24
	if (errornum > 28 || errornum == 0) cout << "Error " << errornum << " Not Defined. Only Errors 1-28 are defined.\n";

	int error[29] = {0};
	double SBCOffset = first.GetTimeSBC() - first.GetTimeSec();
	double SBCTime = veto.GetTimeSBC() - SBCOffset;
	double SBCTimePrev = prev.GetTimeSBC() - SBCOffset;
	veto.UnpackErrorCode(isGood,error);

	// search for particular errors and set a flag.
	if (veto.GetBadScaler() && (veto.GetRun() < 8557 || veto.GetTimeSBC() > 2000000000)) error[28] = true;
	if (firstgoodentry > -1){	//firstgood entry found, if false: can't find errors 18-24
		if (veto.GetTimeSec() > 0 && SBCTime > 0 && SBCOffset !=0 && !veto.GetError(1) && entry > firstgoodentry && fabs((veto.GetTimeSec() - prev.GetTimeSec())-(SBCTime - SBCTimePrev)) > 2) error[18] = true;
		if (veto.GetSEC() == 0 && entry != 0 && entry > firstgoodentry) error[19] = true;
		if (abs(veto.GetSEC() - prev.GetSEC()) > entry-EventNumPrev_good && entry > firstgoodentry && veto.GetSEC()!=0) error[20] = true;
		if (veto.GetQEC() == 0 && entry != 0 && entry > firstgoodentry && !veto.GetError(1)) error[21] = true;
		if (abs(veto.GetQEC() - prev.GetQEC()) > entry-EventNumPrev_good && entry > firstgoodentry && veto.GetQEC() != 0) error[22] = true;
		if (veto.GetQEC2() == 0 && entry != 0 && entry > firstgoodentry && !veto.GetError(1)) error[23] = true;
		if (abs(veto.GetQEC2() - prev.GetQEC2()) > entry-EventNumPrev_good && entry > firstgoodentry && veto.GetQEC2() != 0) error[24] = true;
	}

	return error[errornum];
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