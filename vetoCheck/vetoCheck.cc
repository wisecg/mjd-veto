// vetoCheck for automatic error checking.
// A. Lopez, C. Wiseman
// 9/12/2016
//
// Input: List of run numbers.
// Output: text file & terminal.

#include <iostream>
#include <fstream>

#include "TH1.h"

#include "MJVetoEvent.hh"
#include "GATDataSet.hh"

using namespace std;

bool CheckForBadErrors(MJVetoEvent veto, int entry, int isGood, bool verbose);
int FindQDCThreshold(TH1F *qdcHist, int panel, bool deactivate);
double InterpTime(int entry, vector<double> times, vector<double> entries, vector<bool> badScaler);
void vetoCheck(string Input, int *thresh);

int main(int argc, char* argv[])
{
	if (argc < 2) {
		cout << "Usage:\n ./vetoCheck [run list] ( Optional: [threshold name] )\n\n";
		return 1;
	}
	string file = argv[1];

	// Set first QDC software threshold (it is revised for the second run)
	int thresh[32];
	fill(thresh, thresh + 32, 400);

	// Run performance checks
	vetoCheck(file, thresh);
}

// Search for particular errors and set a flag.
// Ignore errors 4, 7, 10, 11, and 12.
bool CheckForBadErrors(MJVetoEvent veto, int entry, int isGood, bool verbose)
{
	bool badError = false;
	if (isGood != 1)
	{
		int error[18] = {0};
		veto.UnpackErrorCode(isGood,error);

		for (int q=0; q<18; q++)
			if (q != 4 && q != 7  && q != 10 && q != 11 && q != 12 && error[q] == 1)
				badError=true;

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

// Place threshold 35 qdc above pedestal location.
// Also check if panel is active (any entries over QDC = 300)
int FindQDCThreshold(TH1F *qdcHist, int panel, bool deactivate)
{
	int lastNonzeroBin = qdcHist->FindLastBinAbove(1,1);
	if (lastNonzeroBin < 300 && lastNonzeroBin > 0 && deactivate)
	{
		cout << "Panel: " << panel << "  Found last nonzero bin: " << lastNonzeroBin
			 << "  No QDC entries over 300! Deactivating panel ...\n";
		return 9999;
	}
	else if (lastNonzeroBin == -1 && deactivate)
	{
		cout << "Error!  Last nonzero bin == -1\n";
		return 9999;
	}
	int firstNonzeroBin = qdcHist->FindFirstBinAbove(1,1);
	qdcHist->GetXaxis()->SetRange(firstNonzeroBin-10,firstNonzeroBin+50);
	//qdcHist->GetXaxis()->SetRangeUser(0,500); //alternate method of finding pedestal
	int bin = qdcHist->GetMaximumBin();
	if (firstNonzeroBin == -1)
	{
		cout << "ERROR: Panel " << panel << " -- First Nonzero Bin " << firstNonzeroBin
			 << "  Max Bin: " << bin << "  Range: " << firstNonzeroBin-10 << " to " << firstNonzeroBin+50 << endl;
	}
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
		{
			if (badScaler[i] == 0) { upper = times[i]; break; }
		}
		for (int j = entry; j > 0; j--)
		{
			if (badScaler[j] == 0) { lower = times[j]; break; }
		}

		iTime = (upper + lower)/2.0;
	}

	return iTime;
}

void vetoCheck(string Input, int *thresh)
{
	// input a list of run numbers
	ifstream InputList(Input.c_str());
	if(!InputList.good()) {
    	cout << "Couldn't open " << Input << endl;
    	return;
    }
    int filesScanned = 0;	// 1-indexed.

    // output a TXT file
	string Name = Input;
	Name.erase(Name.find_last_of("."),string::npos);
	Name.erase(0,Name.find_last_of("\\/")+1);

	// Char_t OutputFile[200];
	// sprintf(OutputFile,"./vPerf_%s.txt",Name.c_str());
	// ofstream errordat;
	// errordat.open(OutputFile);

	const int nErrs = 27;
	int globalErrorCount[nErrs] = {0};
	int globalRunsWithErrors[nErrs] = {0};
	int RunsWithBadLED = 0;
	int runThresh[32] = {0};	// run-by-run threshold
	int prevThresh[32] = {0};
	long totEntries = 0;
	double totDuration = 0;
	double totLivetime = 0;
	int totErrorCount = 0;
	int globalRunsWithSeriousErrors = 0;
	int globalSeriousErrorCount = 0;
	char hname[50];

	//define lastprevrun vetoevent holder
	MJVetoEvent lastprevrun;	//DO NOT CLEAR

	// ==========================loop over input files==========================
	//
	while(!InputList.eof())
	{
		int run = 0;
		InputList >> run;
		filesScanned++;

		GATDataSet *ds = new GATDataSet(run);
		TChain *v = ds->GetVetoChain();
		long vEntries = v->GetEntries();

		cout << "\n======= Scanning run " << run << ", " << vEntries << " entries. =======\n";

		MJTRun *vRun = new MJTRun();
		MGTBasicEvent *vEvent = new MGTBasicEvent();
		unsigned int mVeto = 0;
		uint32_t vBits = 0;
		v->SetBranchAddress("run",&vRun);
		v->SetBranchAddress("mVeto",&mVeto);
		v->SetBranchAddress("vetoEvent",&vEvent);
		v->SetBranchAddress("vetoBits",&vBits);
		totEntries += vEntries;
		v->GetEntry(0);
		time_t start = vRun->GetStartTime();
		time_t stop = vRun->GetStopTime();
		double duration = (double)(stop - start);
		double livetime = 0;
		totDuration += duration;

		int ErrorCountPerRun[nErrs] = {0};
		vector<double> LocalEntryTime;
		vector<double> LocalEntryNum;
		vector<bool> LocalBadScalers;

		sprintf(hname,"%d_LEDDeltaT",run);
		TH1D *LEDDeltaT = new TH1D(hname,hname,100000,0,100); // 0.001 sec/bin

		TH1F *hRunQDC[32];
		for (int i = 0; i < 32; i++) {
			sprintf(hname,"hRunQDC%d",i);
			hRunQDC[i] = new TH1F(hname,hname,500,0,500);
		}

		MJVetoEvent prev;
		MJVetoEvent first;
		MJVetoEvent last;
		bool foundFirst = false;
		int firstGoodEntry = 0;
		int pureLEDcount = 0;
		bool errorRunBools[nErrs] = {0};
		double xTime = 0;
		double lastGoodTime = 0;
		double SBCOffset = 0;

		// ====================== First loop over entries =========================
		for (int i = 0; i < vEntries; i++)
		{
			v->GetEntry(i);
			MJVetoEvent veto;

			// use the result from FindQDCThreshold
			if (filesScanned > 1) veto.SetSWThresh(runThresh);
			else veto.SetSWThresh(thresh);

	    	int isGood = veto.WriteEvent(i,vRun,vEvent,vBits,run,true); // true: force-write event with errors.

	    	// count up error types
	    	if (isGood != 1)
	    	{
	    		for (int j=0; j<18; j++) if (veto.GetError(j)==1)
	    		{
					ErrorCountPerRun[j]++;
	    			errorRunBools[j]=true;
	    		}
	    	}

	    	// find event time and fill vectors
			if (!veto.GetBadScaler()) {
				LocalBadScalers.push_back(0);
				xTime = veto.GetTimeSec();
			}
			else {
				LocalBadScalers.push_back(1);
				xTime = ((double)i / vEntries) * duration; // this breaks if we have corrupted duration
			}

	    	// fill vectors (the time vectors are revised in the second loop)
			LocalEntryNum.push_back(i);
			LocalEntryTime.push_back(xTime);

			// skip bad entries (true = print contents of skipped event)
	    	if (CheckForBadErrors(veto,i,isGood,false)) continue;

    		// save the first good entry number for the SBC offset
			if (!foundFirst && veto.GetTimeSBC() > 0 && veto.GetTimeSec() > 0 && !veto.GetError(4)) {
				first = veto;
				foundFirst = true;
				firstGoodEntry = i;
			}

			// check if first good entry isn't acutally first good entry
			if (foundFirst && veto.GetError(1))
				foundFirst = false;

	    	// very simple LED tag (fMultip is number of channels above QDC threshold)
			if (veto.GetMultip() > 15) {
				LEDDeltaT->Fill(veto.GetTimeSec()-prev.GetTimeSec());
				pureLEDcount++;
			}

			// end of loop
			prev = veto;
			lastGoodTime = xTime;
			veto.Clear();
		}

		SBCOffset = first.GetTimeSBC() - first.GetTimeSec();

		cout << "First good entry: " << firstGoodEntry
			 << "  Scaler index " << first.GetScalerIndex()
			 << "  Scaler time: " << first.GetTimeSec()
			 << "  SBC time: " << first.GetTimeSBC()-SBCOffset << endl;

		if (duration == 0)
		{
			cout << "Corrupted duration.  Last good timestamp: " << lastGoodTime-first.GetTimeSec() << endl;
			duration = lastGoodTime-first.GetTimeSec();
			totDuration += duration;
		}
		if (stop == 0)
			livetime = lastGoodTime - first.GetTimeSec();
		else
			livetime = duration;
		totLivetime += livetime;
		cout << "Livetime : " << livetime << "  Total livetime: " << totLivetime << endl;

		// find the LED frequency
		double LEDrms = 0;
		double LEDfreq = 0;
		int dtEntries = LEDDeltaT->GetEntries();
		if (dtEntries > 0) {
			int maxbin = LEDDeltaT->GetMaximumBin();
			LEDDeltaT->GetXaxis()->SetRange(maxbin-100,maxbin+100); // looks at +/- 0.1 seconds of max bin.
			LEDrms = LEDDeltaT->GetRMS();
			LEDfreq = 1/LEDDeltaT->GetMean();
		}
		else {
			cout << "Warning! No multiplicity > 15 events!!\n";
			LEDrms = 9999;
			LEDfreq = 9999;
		}
		double LEDperiod = 1/LEDfreq;
		cout << "LED rate: " << LEDfreq << "  Period: " << LEDperiod << " (histo method)\n";
		delete LEDDeltaT;

		// set a flag for "bad LED" (usually a short run causes it)
		// and replace the period with the "simple" one if possible
		bool badLEDFreq = false;
		if (LEDperiod > 9 || vEntries < 100)
		{
			cout << "Warning: Short run.\n";
			if (pureLEDcount > 3) {
				cout << "   From histo method, LED freq is " << LEDfreq
					 << "  Using approximate rate: " << pureLEDcount/duration << endl;
				LEDperiod = duration/pureLEDcount;
			}
			else {
				cout << "   Warning: LED info is corrupted!  Will not use LED period information for this run.\n";
				LEDperiod = 9999;
				badLEDFreq = true;
			}
		}
		if (badLEDFreq) RunsWithBadLED++;
		if (LEDperiod > 9 || LEDperiod < 5) {
			ErrorCountPerRun[25]++;
			errorRunBools[25] = true;
		}

		// ====================== Second loop over entries =========================
		double STime = 0;
		double STimePrev = 0;
		int SIndex = 0;
		int SIndexPrev = 0;
		int EventNumPrev_good = 0;
		int EventNum = 0;
		double SBCTime = 0;
		double SBCTimePrev = 0;
		double TSdifference = 0; // a running total of the time difference between the scaler and SBC timestamps
		prev.Clear();
		pureLEDcount = 0;

		for (int i = 0; i < vEntries; i++)
		{
			EventNum = i;

			// this time we don't skip anything until all errors are checked.
			// we also skip setting QDC thresholds b/c we don't need multiplicity in loop 2.
			v->GetEntry(i);
			MJVetoEvent veto;
	    	int isGood = veto.WriteEvent(i,vRun,vEvent,vBits,run,true);	// true: force-write event with errors.

	    	// find event time
			if (!veto.GetBadScaler())
			{
				xTime = veto.GetTimeSec();
				STime = veto.GetTimeSec();
				SIndex = veto.GetScalerIndex();
				if(run > 8557 && veto.GetTimeSBC() < 2000000000)
					SBCTime = veto.GetTimeSBC() - SBCOffset;
			}
			else if (run > 8557 && veto.GetTimeSBC() < 2000000000)
				xTime = veto.GetTimeSBC() - SBCOffset;
			else
			{
				xTime = InterpTime(i,LocalEntryTime,LocalEntryNum,LocalBadScalers);
				cout << "Entry " << i << ": Bad scaler and SBC times, interpolating : " << xTime << endl;
			}
			LocalEntryTime[i] = xTime;	// replace entry with the more accurate one

			if (veto.GetError(1))
			{
				cout << "Error[1] Missing Packet. Run " << run << "  Entry " << i << endl
					 << "    Scaler index " << veto.GetScalerIndex()
					 << "  Scaler Time " << veto.GetTimeSec()
					 << "  SBC Time " << veto.GetTimeSBC() << "\n\n";
			}

			if (veto.GetSEC() == 0 && i != 0 && i > firstGoodEntry)
			{
				cout << "Error[19] SEC Reset. Run " << run << "  Entry " << i << endl
					 << "    Scaler Index " << veto.GetScalerIndex()
					 << "  SEC " << veto.GetSEC()
					 << "  Previous SEC " << prev.GetSEC() << "\n\n";
				ErrorCountPerRun[19]++;
				errorRunBools[19] = true;
			}

			if (veto.GetQEC() == 0 && i != 0 && i > firstGoodEntry)
			{
				cout << "Error[21] QEC1 Reset.  Run " << run << "  Entry " << i << endl
					 << "    Scaler Index " << veto.GetScalerIndex()
					 << "  QEC1 " << veto.GetQEC()
					 << "  Previous QEC1 " << prev.GetQEC() << "\n\n";
				ErrorCountPerRun[21]++;
				errorRunBools[21] = true;
			}

			if (veto.GetQEC2() == 0 && i != 0 && i > firstGoodEntry)
			{
				cout << "Error[23] QEC1 Reset.  Run " << run << "  Entry " << i << endl
					 << "    Scaler Index " << veto.GetScalerIndex()
					 << "  QEC2 " << veto.GetQEC2()
					 << "  Previous QEC2 " << prev.GetQEC2() << "\n\n";
				ErrorCountPerRun[23]++;
				errorRunBools[23] = true;
			}

			if (abs(veto.GetSEC() - prev.GetSEC()) > EventNum-EventNumPrev_good && i > firstGoodEntry && veto.GetSEC()!=0)
			{
				cout << "Error[20] SEC Change.  Run " << run << "  Entry " << i << endl
					 << "    xTime " << xTime
					 << "  Scaler Index " << veto.GetScalerIndex()
					 << "  SEC " << veto.GetSEC()
					 << "  Previous SEC " << prev.GetSEC() << "\n\n";
				ErrorCountPerRun[20]++;
				errorRunBools[20] = true;
			}

			if(abs(veto.GetQEC() - prev.GetQEC()) > EventNum-EventNumPrev_good && i > firstGoodEntry && veto.GetQEC() != 0)
			{
				cout << "Error[22] QEC1 Jump.  Run " << run << "  Entry " << i << endl
					 << "    xTime " << xTime
					 << "  QDC 1 Index " << veto.GetQDC1Index()
					 << "  QEC 1 " << veto.GetQEC()
					 << "  Previous QEC 1 " << prev.GetQEC() << "\n\n";
				ErrorCountPerRun[22]++;
				errorRunBools[22] = true;
			}

			if(abs(veto.GetQEC2() - prev.GetQEC2()) > EventNum-EventNumPrev_good && i > firstGoodEntry && veto.GetQEC2() != 0)
			{
				cout << "Error[24] QEC2 Jump.  Run " << run << "  Entry " << i << endl
					 << "    xTime " << xTime
					 << "  QDC 2 Index " << veto.GetQDC2Index()
					 << "  QEC 2 " << veto.GetQEC2()
					 << "  Previous QEC 2 " << prev.GetQEC2() << "\n\n";
				ErrorCountPerRun[24]++;
				errorRunBools[24] = true;
			}

			if (STime > 0 && SBCTime >0 && SBCOffset != 0 && !veto.GetError(1)  && i > firstGoodEntry)
			{
				if( fabs((STime - STimePrev) - (SBCTime - SBCTimePrev)) > 2)
				{
					cout << "Error[18] Scaler/SBC Desynch.  Run " << run << "  Entry " << i
						 << "\n    DeltaT (adjusted) " << fabs(fabs(STime - SBCTime) - TSdifference)
						 << "  DeltaT " << fabs(STime - SBCTime)
						 << "  Prev TSdifference " << TSdifference
						 << "  Scaler DeltaT " << STime-STimePrev
						 << "\n    Scaler Index " << SIndex
						 << "  Previous Scaler Index " << SIndexPrev
						 << "  Scaler Time " << STime
						 << "  SBC Time " << SBCTime << "\n\n";
					TSdifference = STime - SBCTime;
					ErrorCountPerRun[18]++;
					errorRunBools[18] = true;
				}
			}

			// save a few things
			TSdifference = STime - SBCTime;
			STimePrev = STime;
			SBCTimePrev = SBCTime;
			SIndexPrev = SIndex;
			STime = 0;
			SBCTime = 0;
			SIndex = 0;

			// skip bad entries before filling QDC
	    	if (CheckForBadErrors(veto,i,isGood,false)) continue; // (true = print contents of skipped event)
	    	for (int j = 0; j < 32; j++) {
				hRunQDC[j]->Fill(veto.GetQDC(j));
			}

			// end of loop
			prev = veto;
			EventNumPrev_good = EventNum; //save last good event number to search for unexpected SEC/QEC changes
			EventNum = 0;
			if (i == vEntries-1){
				last = veto;
				lastprevrun = last;
			}
			veto.Clear();
		}

		// Calculate the run-by-run threshold location.
		// Throw a warning if a pedestal shifts by more than 5%.
		for (int c = 0; c < 32; c++)
		{
			runThresh[c] = FindQDCThreshold(hRunQDC[c],c,false);
			double ratio = (double)runThresh[c]/prevThresh[c];
			if (filesScanned !=1 && (ratio > 1.1 || ratio < 0.9))
			{
				cout << "Warning! Found pedestal shift! Panel: " << c
					 << "  Previous " << prevThresh[c] << "  This run " << runThresh[c] << endl;
				errorRunBools[26] = true;
				ErrorCountPerRun[26]++;
			}
			// save threshold for next scan
			prevThresh[c] = runThresh[c];
		}

		// ADD TO # of RUNS WITH ERROR COUNT
		for (int q = 0; q < nErrs; q++) {
			if (errorRunBools[q]) {
				globalRunsWithErrors[q]++;
			}
		}

		// clear memory
		for (int c=0;c<32;c++) delete hRunQDC[c];

		// calculate # of errors & # of serious errors this run
		int errorsthisrun = 0;
		int seriouserrorsthisrun = 0;
		for (int i = 1; i < nErrs; i++){
			if (i != 10 && i != 11){
				errorsthisrun += ErrorCountPerRun[i];
			}
			if (i == 1 || i > 17) {
				seriouserrorsthisrun += ErrorCountPerRun[i];
				globalSeriousErrorCount += ErrorCountPerRun[i];
			}
		}

		// if this run has errors add to the global counts.
		// globalRunsWithErrors[0] counts how many runs have > 1 error  (print to feresa)
		if (errorsthisrun != 0) globalRunsWithErrors[0]++;
		if (seriouserrorsthisrun != 0) globalRunsWithSeriousErrors++;
		for (int i = 0; i < nErrs; i++) {
			if (ErrorCountPerRun[i] > 0) {
				globalErrorCount[i] += ErrorCountPerRun[i];
				totErrorCount += ErrorCountPerRun[i];
			}
		}

		// ====================== End of run scan ======================
		bool displayRunErrors = false;
		for (int &n : ErrorCountPerRun) if (n > 0) displayRunErrors = true;
		if (displayRunErrors)
		{
			cout << "  End of run scan.  Duration : " << duration << "  Live time: " << livetime << endl;
			for (int i = 0; i < 27; i++){
				if (ErrorCountPerRun[i] > 0)
					cout << "  Error " << i << ": " << ErrorCountPerRun[i] << endl;
			}
			cout << "First Event:"
			     << "\n  SEC " << first.GetSEC() << "  QEC " << first.GetQEC() << "  QEC2 " << first.GetQEC2()
				 << "\n  scaler " << first.GetTimeSec() << "  SBC " << first.GetTimeSBC()-SBCOffset
				 << "\n  scaler index " << first.GetScalerIndex() << endl;
			 cout << "Last Event:"
	 		     << "\n  SEC " << last.GetSEC() << "  QEC " << last.GetQEC() << "  QEC2 " << last.GetQEC2()
	 			 << "\n  scaler " << last.GetTimeSec() << "  SBC " << last.GetTimeSBC()-SBCOffset
	 			 << "\n  scaler index " << last.GetScalerIndex()
				 << "\n  Scaler / SBC time difference at end of run: " << TSdifference << endl;
		}

		// end of run cleanup
		LocalBadScalers.clear();
		LocalEntryNum.clear();
		LocalEntryTime.clear();
	}

	cout << "\n\n================= End of scan. =====================\n"
		 << filesScanned << " runs, " << totEntries
		 << " Total entries,  Total duration (livetime):  "
		 << totDuration << " (" << totLivetime << ") seconds.\n";

	// check for errors and print summary if we find them.
	bool foundErrors = false;
	for (int i = 0; i < nErrs; i++) {
		if (globalErrorCount[i] > 0) {
			foundErrors = true;
			break;
		}
	}

	if (foundErrors)
	{
		cout << "\nError summary:\n"
			 << "Total errors: " << totErrorCount << endl
			 << "Runs with >1 error: " << globalRunsWithErrors[0] << endl
			 << "Number of runs with bad LED: " << RunsWithBadLED << endl
			 << "Number of serious errors: " << globalSeriousErrorCount << endl
			 << "Number of runs w/ serious errors: " << globalRunsWithSeriousErrors << endl;

		for (int i = 0; i < nErrs; i++)
			if (globalErrorCount[i] > 0)
			{
				foundErrors = true;
				cout << i << ": " << globalRunsWithErrors[i] << " runs ("
					 << 100*(double)globalRunsWithErrors[i]/filesScanned << " %)  ";

				if (i != 25)
					cout << "\t" << globalErrorCount[i] << " events ("
					     << 100*(double)globalErrorCount[i]/totEntries << " %)\n";
				else cout << endl;
			}

		// cout << "\nFor reference, error types are: \n"
		//  	 << "1. Missing channels (< 32 veto datas in event) \n"
		//  	 << "2. Extra Channels (> 32 veto datas in event) \n"
		//  	 << "3. Scaler only (no QDC data) \n"
		//  	 << "4. Bad Timestamp: FFFF FFFF FFFF FFFF \n"
		//  	 << "5. QDCIndex - ScalerIndex != 1 or 2 \n"
		//  	 << "6. Duplicate channels (channel shows up multiple times) \n"
		//  	 << "7. HW Count Mismatch (SEC - QEC != 1 or 2) \n"
		//  	 << "8. MJTRun run number doesn't match input file \n"
		//  	 << "9. MJTVetoData cast failed (missing QDC data) \n"
		//  	 << "10. Scaler EventCount doesn't match ROOT entry \n"
		//  	 << "11. Scaler EventCount doesn't match QDC1 EventCount \n"
		//  	 << "12. QDC1 EventCount doesn't match QDC2 EventCount \n"
		//  	 << "13. Indexes of QDC1 and Scaler differ by more than 2 \n"
		//  	 << "14. Indexes of QDC2 and Scaler differ by more than 2 \n"
		//  	 << "15. Indexes of either QDC1 or QDC2 PRECEDE the scaler index \n"
		//  	 << "16. Indexes of either QDC1 or QDC2 EQUAL the scaler index \n"
		//  	 << "17. Unknown Card is present. \n"
		//  	 << "18. Scaler & SBC Timestamp Desynch. \n"
		//  	 << "19. Scaler Event Count reset. \n"
		//  	 << "20. Scaler Event Count increment by > +1. \n"
		//  	 << "21. QDC1 Event Count reset. \n"
		//  	 << "22. QDC1 Event Count increment by > +1. \n"
		//  	 << "23. QDC2 Event Count reset. \n"
		//  	 << "24. QDC2 Event Count increment > +1. \n"
		//  	 << "25. LED frequency very low/high. \n"
		//  	 << "26. Threshold shift by more than 5%. \n";
		//
		// cout << "Serious errors are 1, 18, 19, 20, 21, 22, 23, 24, 25, 26." << endl;
	}

	// errordat.close();
}