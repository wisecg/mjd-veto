// Finds muons.  Optionally writes some output.
// Clint Wiseman, USC/Majorana
// 3/9/2016

#include "vetoScan.hh"

using namespace std;

void muFinder(string Input, int *thresh, bool root, bool list)
{
	// LED Cut Parameters (C-f "Display Cut Parameters" below.)
	double LEDWindow = 0.1;
	int LEDMultipThreshold = 10;  // "multipThreshold" = "highestMultip" - "LEDMultipThreshold"
	int LEDSimpleThreshold = 20;  // used when LED frequency measurement is bad.

	// Custom SW Threshold (obtained from vetoThreshFinder)
	int swThresh[32] = {0};
	if (thresh != NULL) {
		cout << "muFinder is using these SW thresholds: " << endl;
		memcpy(swThresh,thresh,sizeof(swThresh));
		for (int j=0;j<32;j++) cout << j << ":" << swThresh[j] << " ";
		cout << endl;
	}
	else {
		cout << "thresh is NULL!  Using default." << endl;
		for (int j=0;j<32;j++) swThresh[j] = 500;
	}

	// Input a list of run numbers
	ifstream InputList(Input.c_str());
	if(!InputList.good()) {
		cout << "Couldn't open " << Input << endl;
		return;
	}

	// Set up output files
	string Name = Input;
	Name.erase(Name.find_last_of("."),string::npos);
	Name.erase(0,Name.find_last_of("\\/")+1);

	// Output 1: Text file muon list (used in skim files)
	string outName = "./output/MuonList_"+Name+".txt";
	ofstream MuonList(outName.c_str());
	if (!list) MuonList.close();

	// Output 2: ROOT output
	int isGood;
	int run;
	long rEntry;
	long start;
	long stop;
	long prevStopTime = 0;
	double duration;
	int CoinType[32];
	int CutType[32];
	int PlaneHits[12];
	int PlaneTrue[12];
	int PlaneHitCount = 0;
	int highestMultip = 0;
	int multipThreshold = 0;
	double LEDfreq = 0;
	double LEDrms = 0;
	double xTime = 0;
	double x_deltaT = 0;
	double x_LEDDeltaT = 0;
	double timeSBC = 0;
	Char_t OutputFile[200];
	sprintf(OutputFile,"./output/%s.root",Name.c_str());
	TFile *RootFile = new TFile(OutputFile, "RECREATE");
  	TH1::AddDirectory(kFALSE); // Global flag: "When a (root) file is closed, all histograms in memory associated with this file are automatically deleted."
	TTree *vetoEvent = new TTree("vetoEvent","MJD Veto Events");
	MJVetoEvent out;
	if (root) {
		vetoEvent->Branch("events","MJVetoEvent",&out,32000,1);
		vetoEvent->Branch("rEntry",&rEntry,"rEntry/L");
		vetoEvent->Branch("timeSBC",&timeSBC);
		vetoEvent->Branch("LEDfreq",&LEDfreq);
		vetoEvent->Branch("LEDrms",&LEDrms);
		vetoEvent->Branch("multipThreshold",&multipThreshold);
		vetoEvent->Branch("highestMultip",&highestMultip);
		vetoEvent->Branch("LEDWindow",&LEDWindow);
		vetoEvent->Branch("LEDMultipThreshold",&LEDMultipThreshold);
		vetoEvent->Branch("LEDSimpleThreshold",&LEDSimpleThreshold); // add a bool signifiying it's in use
		vetoEvent->Branch("start",&start,"start/L");
		vetoEvent->Branch("stop",&stop,"stop/L");
		vetoEvent->Branch("duration",&duration);
		vetoEvent->Branch("xTime",&xTime);
		vetoEvent->Branch("x_deltaT",&x_deltaT);
		vetoEvent->Branch("x_LEDDeltaT",&x_LEDDeltaT);
		vetoEvent->Branch("CoinType[32]",CoinType,"CoinType[32]/I");
		vetoEvent->Branch("CutType[32]",CutType,"CutType[32]/I");
		vetoEvent->Branch("PlaneHits[12]",PlaneHits,"PlaneHits[12]/I");
		vetoEvent->Branch("PlaneTrue[12]",PlaneTrue,"PlaneTrue[12]/I");
		vetoEvent->Branch("PlaneHitCount",&PlaneHitCount);
	}

	// Loop over files.
	int JumpCount = 0;	// scaler jump counter
	while(!InputList.eof()){

		// initialize
		InputList >> run;
		GATDataSet *ds = new GATDataSet(run);
		TChain *v = ds->GetVetoChain();
		long vEntries = v->GetEntries();
		MJTRun *vRun = new MJTRun();
		MGTBasicEvent *vEvent = new MGTBasicEvent();
		unsigned int mVeto = 0;
		uint32_t vBits = 0;
		v->SetBranchAddress("run",&vRun);
		v->SetBranchAddress("mVeto",&mVeto);
		v->SetBranchAddress("vetoEvent",&vEvent);
		v->SetBranchAddress("vetoBits",&vBits);
		v->GetEntry(0);
		start = (long)vRun->GetStartTime();
		stop = (long)vRun->GetStopTime();
		duration = ds->GetRunTime()/CLHEP::second;

		printf("\n======= Scanning run %i, %li entries, %.0f sec. =======\n",run,vEntries,duration);
		cout << "start: " << start << "  stop: " << stop << endl;

		// ========= 1st loop over veto entries - Measure LED frequency. =========
		//
		// Goal is to measure the LED frequency, to be used in the second loop as a
		// time cut. This is done by finding the maximum bin of a delta-t histogram.
		// This section of the code uses a weak multiplicity threshold of 20 -- it
		// doesn't need to be exact, and should also work for runs where there were
		// only 24 panels installed.
		//
		bool badLEDFreq = false;
		MJVetoEvent prev;
		char hname[200];
		sprintf(hname,"LEDDeltaT_run%i",run);
		TH1F *LEDDeltaT = new TH1F(hname,hname,100000,0,100); // 0.001 sec/bin
		highestMultip = 0;	// try to predict how many panels there are for this run.
		long skippedEvents = 0;
		long corruptScaler = 0;
		bool foundFirst = false;
		int firstGoodEntry = 0;
		MJVetoEvent first;
		highestMultip=0;
		for (long i = 0; i < vEntries; i++)
		{
			v->GetEntry(i);
			MJVetoEvent veto;
			veto.SetSWThresh(swThresh);
	    	isGood = veto.WriteEvent(i,vRun,vEvent,vBits,run,true);
    		if (CheckForBadErrors(veto,i,isGood,false)) {
    			skippedEvents++;
    			continue;
    		}

	    	if (veto.GetBadScaler()) corruptScaler++;

	    	if (veto.GetMultip() > highestMultip && veto.GetMultip() < 33) {
	    		highestMultip = veto.GetMultip();
	    		cout << "Finding highest multiplicity: " << highestMultip << "  entry: " << i << endl;
	    	}

	    	// Save the first good entry number for the SBC offset time
			if (isGood && !foundFirst && veto.GetTimeSBC()>0.01 && veto.GetTimeSec()>0.01 && !veto.GetBadScaler()) {
				first = veto;
				foundFirst = true;
				firstGoodEntry = i;
			}

			// if (i > 200 && i < 350)
			// printf("scaler %.2f  sbc %.2f\n",veto.GetTimeSec(),veto.GetTimeSBC());


	    	// Very simple LED tag.
			if (veto.GetMultip() >= 20) {
				LEDDeltaT->Fill(veto.GetTimeSec()-prev.GetTimeSec());
			}
			prev = veto;
		}
		// Find the SBC offset
		double SBCOffset = first.GetTimeSBC() - first.GetTimeSec();
		printf("First good entry: %i  Scaler %.2f  SBC %.2f  SBCOffset %.2f\n"
			,firstGoodEntry,first.GetTimeSec(),first.GetTimeSBC(),SBCOffset);

		// Find the LED frequency
		if (skippedEvents > 0) printf("Skipped %li of %li entries.\n",skippedEvents,vEntries);
		// if (corruptScaler > 0) printf("Corrupt scaler: %li of %li entries (%.2f%%) .\n"
			// ,corruptScaler,vEntries,100*(double)corruptScaler/vEntries);

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
		//
		prev.Clear();
		MJVetoEvent prevLED;
		double xTimePrev = 0;
		double x_deltaTPrev = 0;
		double xTimePrevLED = 0;
		double xTimePrevLEDSimple = 0;
		bool firstLED = false;
		// bool IsLEDPrev = false;
		int almostMissedLED = 0;
		double TSdifference = 0;
		for (long i = 0; i < vEntries; i++)
		// for (long i = 250; i < 300; i++)
		{
			v->GetEntry(i);
			rEntry = i;	// save ROOT entry in output
			MJVetoEvent veto;
			veto.SetSWThresh(swThresh);
    	isGood = veto.WriteEvent(i,vRun,vEvent,vBits,run,true);
			timeSBC = veto.GetTimeSBC()-SBCOffset;

    	//----------------------------------------------------------
			// 0: Time of event and skipping if necessary.
			// Employ alternate methods if the scaler is corrupted.
			// Should implement an estimate of the error when alternate methods are used.
			//
			bool ApproxTime = false;

			xTime = -1;

			if (!veto.GetBadScaler())
			{
				xTime = veto.GetTimeSec();

				// Find scaler jumps and adjust xTime by "TSdifference"
				// TSdifference starts at 0 at the beginning of the run.
				if (veto.GetTimeSec() != 0 && veto.GetTimeSBC() !=0 && SBCOffset != 0 && i>=firstGoodEntry)
				{
					double sbc = veto.GetTimeSBC() - SBCOffset;
					double diff = veto.GetTimeSec() - sbc;

					// if (fabs(fabs(diff) - TSdifference) > 1)	// andrew's original method (11472 - fails)
					if (fabs(diff-TSdifference) > 1)	// clint's method (11472 bkwds - OK)
					{
						JumpCount++;
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
	    	if (CheckForBadErrors(veto,i,isGood,false))
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

			//----------------------------------------------------------
			// 1. LED Cut
			//
			// TRUE if an event PASSES (i.e. is physics.)  FALSE if an event is an LED.
			//
			// If LED's are turned off or the frequency measurement is bad, we revert
			// to a simple multiplicity threshold.
			//
			bool TimeCut = true;
			bool IsLED = false;

			// Set Cut
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

			//----------------------------------------------------------
	    	// 2: Energy (Gamma) Cut
	    	// The measured muon energy threshold is QDC = 500.
	    	// Set TRUE if at least TWO panels are over 500.
	    	//
	    	bool EnergyCut = false;

	    	int over500Count = 0;
	    	for (int q = 0; q < 32; q++) {
	    		if (veto.GetQDC(q) > 500)
	    			over500Count++;
	    	}
	    	if (over500Count >= 2) EnergyCut = true;

			//----------------------------------------------------------
			// 3: Hit Pattern
			// Map hits above SW threshold to planes and count the hits.
			//

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

			//----------------------------------------------------------
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

				// Other coincidence types can be found by parsing the ROOT output.
			}

			//----------------------------------------------------------
			// 5: Output
			// The skim file used to take a text file of muon candidate events.
			// Additionally, write the ROOT file containing all the real data.
			//

			// Write a text file
			if (list) {
				char buffer[200];
				if (CoinType[1] || CoinType[0]) {
					int type;
					if (CoinType[0]) type = 1;
					if (CoinType[1]) type = 2;
					sprintf(buffer,"%i %li %.8f %i %i\n",run,start,xTime,type,veto.GetBadScaler());
					MuonList << buffer;
				}
				// This is Jason's TYPE 3: flag runs with gaps since the last stop time.
				if ((start - prevStopTime) > 10 && i == 0) {
					sprintf(buffer,"%i %li 0.0 3 0\n",run,start);
					MuonList << buffer;
				}
			}

			// Assign all bools calculated to the int array CutType[32];
			CutType[0] = LEDTurnedOff;
			CutType[1] = EnergyCut;
			CutType[2] = ApproxTime;
			CutType[3] = TimeCut;
			CutType[4] = IsLED;
			CutType[5] = firstLED;
			CutType[6] = badLEDFreq;

			// Write ROOT output
			if (root) {
				out = veto;
				vetoEvent->Fill();
			}

			// Reset for next entry
			//----------------------------------------------------------
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

	    // End of run summaries.
		if (almostMissedLED > 0) cout << "\nWarning, almost missed " << almostMissedLED << " LED events.\n";

	    // done with this run.
		delete ds;
		prevStopTime = stop;
	}

	printf("\n===================== End of Scan. =====================\n");

	if (JumpCount > 0) cout << "\nWarning, found " << JumpCount << " scaler jumps.\n";

	if (list) MuonList.close();
	if (root) vetoEvent->Write();
	RootFile->Close();
}
