// vetoCheck for automatic error checking.
// A. Lopez, C. Wiseman
// 9/12/2016
//
// To be used in auto-processing, run by run.
// Takes one input argument, a run number.
//
// Known Error types:
// 1. Missing channels (< 32 veto datas in event)
// 2. Extra Channels (> 32 veto datas in event)
// 3. Scaler only (no QDC data)
// 4. Bad Timestamp: FFFF FFFF FFFF FFFF
// 5. QDCIndex - ScalerIndex != 1 or 2
// 6. Duplicate channels (channel shows up multiple times)
// 7. HW Count Mismatch (SEC - QEC != 1 or 2)
// 8. MJTRun run number doesn't match input file
// 9. MJTVetoData cast failed (missing QDC data)
// 10. Scaler EventCount doesn't match ROOT entry
// 11. Scaler EventCount doesn't match QDC1 EventCount
// 12. QDC1 EventCount doesn't match QDC2 EventCount
// 13. Indexes of QDC1 and Scaler differ by more than 2
// 14. Indexes of QDC2 and Scaler differ by more than 2
// 15. Indexes of either QDC1 or QDC2 PRECEDE the scaler index
// 16. Indexes of either QDC1 or QDC2 EQUAL the scaler index
// 17. Unknown Card is present.
// 18. Scaler & SBC Timestamp Desynch.
// 19. Scaler Event Count reset.
// 20. Scaler Event Count increment by > +1.
// 21. QDC1 Event Count reset.
// 22. QDC1 Event Count increment by > +1.
// 23. QDC2 Event Count reset.
// 24. QDC2 Event Count increment > +1.
// 25. LED frequency very low/high or corrupted.
// 26. QDC threshold not found
// 27. No events above QDC threshold
// 28. Used interpolated time
//
// Serious errors are: 1, 13, 14, 18, 19, 20, 21, 22, 23, 24, 25

#include <iostream>
#include <fstream>
#include "TH1.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TLine.h"
#include "MJVetoEvent.hh"
#include "GATDataSet.hh"

using namespace std;

bool CheckForBadErrors(MJVetoEvent veto, int entry, int isGood, bool verbose);
double InterpTime(int entry, vector<double> times, vector<double> entries, vector<bool> badScaler);
int FindQDCThreshold(TH1F *qdcHist);
void vetoCheck(int run, bool draw);

int main(int argc, char* argv[])
{
	if (argc < 2) {
		cout << "Usage:\n ./vetoCheck [run number] ([-d] draws qdc plot)\n\n";
		return 1;
	}
	int run = atoi(argv[1]);

	bool draw = false;
	string opt = "";
	if (argc > 2) opt = argv[2];
	if (argc > 2 && opt == "-d")
		draw = true;

	vetoCheck(run,draw);
}

void vetoCheck(int run, bool draw)
{
	const int nErrs = 29; // error 0 is unused
	int SeriousErrorCount = 0;
	int TotalErrorCount = 0;
	bool badLEDFreq = false;

	// Specify which error types to print during the loop over events
	vector<int> SeriousErrors = {1, 13, 14, 18, 19, 20, 21, 22, 23, 24};

	GATDataSet *ds = new GATDataSet(run);
	TChain *v = ds->GetVetoChain();
	long vEntries = v->GetEntries();

	cout << "===== Scanning veto data, run " << run << ", " << vEntries << " entries. ======\n";

	// Suppress the "Error in <TClass::LoadClassInfo>" messages
	gROOT->ProcessLine( "gErrorIgnoreLevel = 3001;");

	MJTRun *vRun = new MJTRun();
	MGTBasicEvent *vEvent = new MGTBasicEvent();
	unsigned int mVeto = 0;
	uint32_t vBits = 0;
	v->SetBranchAddress("run",&vRun);
	v->SetBranchAddress("mVeto",&mVeto);
	v->SetBranchAddress("vetoEvent",&vEvent);
	v->SetBranchAddress("vetoBits",&vBits);
	v->GetEntry(0);
	time_t start = vRun->GetStartTime();
	time_t stop = vRun->GetStopTime();
	double duration = (double)(stop - start);
	double livetime = 0;

	int ErrorCount[nErrs] = {0};
	vector<double> EntryTime;
	vector<double> EntryNum;
	vector<bool> BadScalers;

	char hname[50];
	sprintf(hname,"%d_LEDDeltaT",run);
	TH1D *LEDDeltaT = new TH1D(hname,hname,100000,0,100); // 0.001 sec/bin
	TH1F *hRunQDC[32];
	for (int i = 0; i < 32; i++) {
		sprintf(hname,"hRunQDC%d",i);
		hRunQDC[i] = new TH1F(hname,hname,4200,0,4200);
	}

	MJVetoEvent prev;
	MJVetoEvent first;
	MJVetoEvent last;
	bool foundFirst = false;
	bool foundFirstSTS = false;
	int firstGoodEntry = 0;
	int pureLEDcount = 0;
	bool Error[nErrs] = {0};
	double xTime = 0;
	double lastGoodTime = 0;
	double SBCOffset = 0;
	double firstGoodSTS = 0;	//to accurately calculate duration and livetime if start/stop = 0

	// ====================== First loop over entries =========================
	for (int i = 0; i < vEntries; i++)
	{
		v->GetEntry(i);
		MJVetoEvent veto;

		// Set QDC software threshold (used for multiplicity calculation)
		int thresh[32];
		fill(thresh, thresh + 32, 400);
		veto.SetSWThresh(thresh);

		// true: force-write an event with errors.
    	int isGood = veto.WriteEvent(i,vRun,vEvent,vBits,run,true);

    	// find event time and fill vectors
		if (!veto.GetBadScaler()) {
			BadScalers.push_back(0);
			xTime = veto.GetTimeSec();
		}
		else {
			BadScalers.push_back(1);
			xTime = ((double)i / vEntries) * duration; // this breaks if we have corrupted duration
		}

    	// fill vectors (the time vectors are revised in the second loop)
		EntryNum.push_back(i);
		EntryTime.push_back(xTime);

		// check if first good entry isn't acutally first good entry
		if (foundFirst && veto.GetError(1))
			foundFirst = false;

		//check for first good scaler time stamp
		if (!foundFirstSTS && !veto.GetError(4)){
			foundFirstSTS = true;
			firstGoodSTS = veto.GetTimeSec();
		}

		// skip bad entries (true = print contents of skipped event)
    	if (CheckForBadErrors(veto,i,isGood,false)) continue;

		// save the first good entry number for the SBC offset
		if (!foundFirst && veto.GetTimeSBC() > 0 && veto.GetTimeSec() > 0 && !veto.GetError(4)) {
			first = veto;
			foundFirst = true;
			firstGoodEntry = i;
		}

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

	if (duration == 0)
	{
		cout << "Corrupted duration.  Last good timestamp: " << lastGoodTime-firstGoodSTS << endl;
		duration = lastGoodTime-firstGoodSTS;
	}
	livetime = duration - (first.GetTimeSec() - firstGoodSTS);

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
		cout << "Warning! No multiplicity > 15 events.  LED may be off.\n";
		LEDrms = 9999;
		LEDfreq = 9999;
		badLEDFreq = true;
	}
	double LEDperiod = 1/LEDfreq;
	delete LEDDeltaT;
	if (LEDperiod > 9 || vEntries < 100)
	{
		cout << "Warning: Short run.\n";
		if (pureLEDcount > 3) {
			cout << "   From histo method, LED freq is " << LEDfreq
				 << "  Using approximate rate: " << pureLEDcount/duration << endl;
			LEDperiod = duration/pureLEDcount;
		}
		else {
			LEDperiod = 9999;
			badLEDFreq = true;
		}
	}
	if (LEDperiod > 9 || LEDperiod < 5 || badLEDFreq) {
		ErrorCount[25]++;
		Error[25] = true;
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
    	veto.WriteEvent(i,vRun,vEvent,vBits,run,true);	// true: force-write event with errors.

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
			xTime = InterpTime(i,EntryTime,EntryNum,BadScalers);
		 	Error[28] = true;
 			ErrorCount[28]++;
		}
		EntryTime[i] = xTime;	// replace entry with the more accurate one

		// Check for errors
		for (int j=0; j<18; j++) if (veto.GetError(j)==1)
		{
			ErrorCount[j]++;
			Error[j]=true;
		}
		Error[18] = (bool)((STime > 0 && SBCTime > 0 && SBCOffset != 0 && !veto.GetError(1) && i > firstGoodEntry) && fabs((STime - STimePrev) - (SBCTime - SBCTimePrev)) > 2);
		Error[19] = (bool)(veto.GetSEC() == 0 && i != 0 && i > firstGoodEntry);
		Error[20] = (bool)(abs(veto.GetSEC() - prev.GetSEC()) > EventNum-EventNumPrev_good && i > firstGoodEntry && veto.GetSEC()!=0);
		Error[21] = (bool)(veto.GetQEC() == 0 && i != 0 && i > firstGoodEntry && !veto.GetError(1));
		Error[22] = (bool)(abs(veto.GetQEC() - prev.GetQEC()) > EventNum-EventNumPrev_good && i > firstGoodEntry && veto.GetQEC() != 0);
		Error[23] = (bool)(veto.GetQEC2() == 0 && i != 0 && i > firstGoodEntry && !veto.GetError(1));
		Error[24] = (bool)(abs(veto.GetQEC2() - prev.GetQEC2()) > EventNum-EventNumPrev_good && i > firstGoodEntry && veto.GetQEC2() != 0);

		// Print errors to screen
		bool PrintError = false;
		for (auto i : SeriousErrors){
			if (Error[i]) {
				PrintError = true;
				break;
			}
		}
		// 1, 13, 14, 18, 19, 20, 21, 22, 23, 24 - must match vector "SeriousErrors"
		if (PrintError)
		{
			cout << "\nSerious errors found in entry " << i << ":\n";
			if (Error[1]) {
				cout << "Error[1] Missing Packet."
					 << "  Scaler index " << veto.GetScalerIndex()
					 << "  Scaler Time " << veto.GetTimeSec()
					 << "  SBC Time " << veto.GetTimeSBC() << "\n";
			}
			if (Error[13]) {
				cout << "Error[13] ORCA packet indexes of QDC1 and Scaler differ by more than 2."
					 << "\n    Scaler Index " << veto.GetScalerIndex()
					 << "  QDC1 Index " << veto.GetQDC1Index()
					 << "\n    Previous scaler Index " << prev.GetScalerIndex()
					 << "  Previous QDC1 Index " << prev.GetQDC1Index() << endl;
				ErrorCount[13]++;
			}
			if (Error[14]) {
				cout << "Error[14] ORCA packet indexes of QDC2 and Scaler differ by more than 2."
					 << "\n    Scaler Index " << veto.GetScalerIndex()
					 << "  QDC2 Index " << veto.GetQDC2Index()
					 << "\n    Previous scaler Index " << prev.GetScalerIndex()
					 << "  Previous QDC2 Index " << prev.GetQDC2Index() << endl;
				ErrorCount[14]++;
			}
			if (Error[18]) {
				cout << "Error[18] Scaler/SBC Desynch."
					 << "\n    DeltaT (adjusted) " << STime - SBCTime - TSdifference
					 << "  DeltaT " << STime - SBCTime
					 << "  Prev TSdifference " << TSdifference
					 << "  Scaler DeltaT " << STime-STimePrev
					 << "\n    Scaler Index " << SIndex
					 << "  Previous Scaler Index " << SIndexPrev
					 << "  Scaler Time " << STime
					 << "  SBC Time " << SBCTime << "\n";
				TSdifference = STime - SBCTime;
				ErrorCount[18]++;
			}
			if (Error[19]) {
				cout << "Error[19] Scaler Event Count Reset. "
					 << "  Scaler Index " << veto.GetScalerIndex()
					 << "  SEC " << veto.GetSEC()
					 << "  Previous SEC " << prev.GetSEC() << "\n";
				ErrorCount[19]++;
			}
			if (Error[20]) {
				cout << "Error[20] Scaler Event Count Jump."
					 << "    xTime " << xTime
					 << "  Scaler Index " << veto.GetScalerIndex()
					 << "  SEC " << veto.GetSEC()
					 << "  Previous SEC " << prev.GetSEC() << "\n";
				ErrorCount[20]++;
			}
			if (Error[21]) {
				cout << "Error[21] QDC1 Event Count Reset."
					 << "  Scaler Index " << veto.GetScalerIndex()
					 << "  QEC1 " << veto.GetQEC()
					 << "  Previous QEC1 " << prev.GetQEC() << "\n";
				ErrorCount[21]++;
			}
			if(Error[22]) {
				cout << "Error[22] QDC 1 Event Count Jump."
					 << "  xTime " << xTime
					 << "  QDC 1 Index " << veto.GetQDC1Index()
					 << "  QEC 1 " << veto.GetQEC()
					 << "  Previous QEC 1 " << prev.GetQEC() << "\n";
				ErrorCount[22]++;
			}
			if (Error[23]) {
				cout << "Error[23] QDC2 Event Count Reset."
					 << "  Scaler Index " << veto.GetScalerIndex()
					 << "  QEC2 " << veto.GetQEC2()
					 << "  Previous QEC2 " << prev.GetQEC2() << "\n";
				ErrorCount[23]++;
			}
			if(Error[24]) {
				cout << "Error[24] QDC 2 Event Count Jump."
					 << "  xTime " << xTime
					 << "  QDC 2 Index " << veto.GetQDC2Index()
					 << "  QEC 2 " << veto.GetQEC2()
					 << "  Previous QEC 2 " << prev.GetQEC2() << "\n";
				ErrorCount[24]++;
			}
			cout << endl;
		}

		TSdifference = STime - SBCTime;
		STimePrev = STime;
		SBCTimePrev = SBCTime;
		SIndexPrev = SIndex;
		STime = 0;
		SBCTime = 0;
		SIndex = 0;
		prev = veto;
		last = prev;
		EventNumPrev_good = EventNum; //save last good event number to search for unexpected SEC/QEC changes
		EventNum = 0;

		// Reset error bools each entry
		for (int j=0; j<nErrs; j++){
			Error[j]=false;
		}	
		
		
		// Skip bad entries before filling QDC.
		if (PrintError) {
			veto.Clear();
			continue;
		}
		
		for (int j = 0; j < 32; j++) hRunQDC[j]->Fill(veto.GetQDC(j));
	}

	// Find QDC threshold and make sure we have counts above pedestal.
	// If -d option is used, print a graph.
	TCanvas *can = new TCanvas("can","veto QDC thresholds, panels 1-32",0,0,800,600);
	can->Divide(8,4,0,0);
	for (int i = 0; i < 32; i++)
	{
		int thresh = FindQDCThreshold(hRunQDC[i]);
		if (thresh == -1) {
			Error[26] = true;
			ErrorCount[26]++;
			cout << "Error[26]: QDC threshold for panel " << i << " not found." << endl;
		}
		int CountsAboveThresh = hRunQDC[i]->Integral(thresh,4200);
		if (CountsAboveThresh == 0) {
			Error[27] = true;
			ErrorCount[27]++;
			cout << "Error[27]: No counts above QDC threshold " << thresh << " for panel " << i << endl;
		}

		if (draw) {
			can->cd(i+1);
			TVirtualPad *vpad0 = can->cd(i+1); vpad0->SetLogy();

			hRunQDC[i]->GetXaxis()->SetRange(0,500);
			hRunQDC[i]->Draw();

			double ymax = hRunQDC[i]->GetMaximum();
			TLine *line = new TLine(thresh,0,thresh,ymax+10);
			line->SetLineColor(kRed);
			line->SetLineWidth(2.0);
			line->Draw();
		}
	}

	sprintf(hname,"QDC_run%i.pdf",run);
	if (draw) can->Print(hname);

	// Calculate total errors and total serious errors
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
	cout << "=================== End veto scan. ======================\n";
	if (SeriousErrorCount > 0)
	{
		cout << "Total Errors : " << TotalErrorCount << endl;
		cout << "Serious Errors : " << SeriousErrorCount << endl;

		if (duration != livetime)
			cout << "Run duration (" << duration << " sec) doesn't match live time: " << livetime << endl;

		for (int i = 0; i < 27; i++)
		{
			if (ErrorCount[i] > 0)
			{
				if (i != 25)
					cout << "  Error[" << i <<"]: " << ErrorCount[i] << " events ("
						 << 100*(double)ErrorCount[i]/vEntries << " %)\n";
		 		else if (i == 25)
				{
					cout << "  Error[25]: Bad LED rate: " << LEDfreq << "  Period: " << LEDperiod << endl;
					if (LEDperiod > 0.1 && (abs(duration/LEDperiod) - pureLEDcount) > 5)
					{
						cout << "   LED count: " << pureLEDcount
							 << "  Expected :" << (int)duration/(int)LEDperiod << endl;
					}
				}
			}
		}

		cout << "\n  For reference, \"serious\" error types are: ";
		for (auto i : SeriousErrors) cout << i << " ";
		cout << "\n  Please report these to the veto group.\n";

		cout << "================= End veto error report. =================\n";
	}
}

// ================================================================================
// ================================================================================

// Check the 18 built-in error types in a MJVetoEvent object
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
int FindQDCThreshold(TH1F *qdcHist)
{
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