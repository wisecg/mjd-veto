// auto-veto.cc
// Runs during production, creates ROOT files of MJD veto data.
// C. Wiseman, A. Lopez, 10/24/16.
//
// NOTE: The scans are split across a few different loops over the events in the run.
// This is done to increase the flexibility of the code, since it checks many different
// quantities.  The performance hit is minimal, since the size of the veto trees
// is relatively small.

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"
#include "TFile.h"
#include "MJVetoEvent.hh"
#include "GATDataSet.hh"
#include "MGTEvent.hh"
#include "MGVDigitizerData.hh"

using namespace std;

vector<int> MeasurePanelThresholds(TChain *vetoChain, string outputDir, bool makePlots=false);
void ProcessVetoData(TChain *vetoChain, vector<int> thresholds, string outputDir, bool errorCheckOnly=false);

void SetCardNumbers(int runNum, int &card1, int &card2);
int FindThreshold(TH1D *qdcHist, int threshVal, int panel, int runNum);
int PlaneMap(int qdcChan, int runNum=0);
bool CheckErrors(MJVetoEvent veto, MJVetoEvent prev, vector<int> &ErrorVec);
bool CheckErrors(MJVetoEvent veto, MJVetoEvent prev);

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
	GATDataSet ds;
	string runPath = ds.GetPathToRun(run,GATDataSet::kBuilt);
	// string runPath = "./stage/OR_run"+std::to_string(run)+".root"; // manually set path

  TChain *vetoChain = new TChain("VetoTree");
	if (!vetoChain->Add(runPath.c_str())){
		cout << "File doesn't exist.  Exiting ...\n";
		return 1;
	}

	printf("\n========= Processing run %i ... %lli entries. =========\n",run,vetoChain->GetEntries());
	cout << "Path: " << runPath << endl;

	// Find the QDC pedestal location in each channel.
	// Set a software threshold value above this location,
	// and optionally output plots that confirm this choice.
	vector<int> thresholds = MeasurePanelThresholds(vetoChain, outputDir, draw);

	// Check for data quality errors,
	// tag muon and LED events in veto data,
	// and output a ROOT file for further analysis.
	ProcessVetoData(vetoChain,thresholds,outputDir,errorCheckOnly);

	printf("=================== Done processing. ====================\n\n");
	return 0;
}

vector<int> MeasurePanelThresholds(TChain *vetoChain, string outputDir, bool makePlots)
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
	long skippedEvents = 0;

	// MJVetoEvent variables, with run-based card numbers
	int card1=0, card2=0;
	SetCardNumbers(runNum,card1,card2);
	MJVetoEvent veto(card1,card2);
	MJVetoEvent prev, first;
  // printf("VME Slots: QDC1 %i  QDC2 %i\n",card1,card2);

	while (reader.Next())
	{
		long i = reader.GetCurrentEntry();
		int run = vRun->GetRunNumber();
		veto.SetSWThresh(def);
    veto.WriteEvent(i,&*vRun,&*vEvt,*vBits,run,true);
		if (CheckErrors(veto,prev))
		{
    	skippedEvents++;
			// do the end of event resets before continuing
			prev = veto;
			veto.Clear();
    		continue;
    	}
    	for (int q = 0; q < 32; q++) {
    		hLowQDC[q]->Fill(veto.GetQDC(q));
    		hFullQDC[q]->Fill(veto.GetQDC(q));
		}
		// save previous entries for the event error check
		prev = veto;
		veto.Clear();
	}
	if (skippedEvents > 0) printf("MeasurePanelThresholds skipped %li of %li entries.\n",skippedEvents,vEntries);

	int thresh[32] = {9999};
	for (int i = 0; i < 32; i++)
	{
		thresh[i] = FindThreshold(hLowQDC[i],threshVal,i,runNum);
		thresholds.push_back(i);
		thresholds.push_back(thresh[i]);
	}
	// cout << "Found thresholds: " << endl;
	// for (int i = 0; i < 32; i++)
		// cout << i << " " << thresh[i] << endl;

  if (makePlots)
	{
  	// re-scan with the found thresholds to make a multiplicity plot
  	reader.SetTree(vetoChain);  // resets the reader
  	while (reader.Next())
  	{
  		long i = reader.GetCurrentEntry();
  		int run = vRun->GetRunNumber();
  		veto.SetSWThresh(thresh);
      veto.WriteEvent(i,&*vRun,&*vEvt,*vBits,run,true);
  		if (CheckErrors(veto,prev))
  		{
      	skippedEvents++;
  			// do the end of event resets before continuing
  			prev = veto;
  			veto.Clear();
      		continue;
      	}
  		hMultip->Fill(veto.GetMultip());
  		// save previous entries for the event error check
  		prev = veto;
  		veto.Clear();
  	}
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

		c1->Print(TString::Format("%s/veto-%i-qdc.pdf",outputDir.c_str(),runNum));
		c2->Print(TString::Format("%s/veto-%i-qdcThresh.pdf",outputDir.c_str(),runNum));
		c3->Print(TString::Format("%s/veto-%i-multip.pdf",outputDir.c_str(),runNum));
	}
	return thresholds;
}

void ProcessVetoData(TChain *vetoChain, vector<int> thresholds, string outputDir, bool errorCheckOnly)
{
	// QDC software threshold (obtained from MeasurePanelThresholds)
	int swThresh[32] = {0};
	for (int i = 0; i < (int)thresholds.size(); i+=2)
		swThresh[thresholds[i]] = thresholds[i+1];

	// LED variables
	int LEDMultipThreshold=5;  // "multipThreshold" = "highestMultip" - "LEDMultipThreshold"
	int LEDSimpleThreshold=10;  // used when LED frequency measurement is bad.
	int highestMultip=0, multipThreshold=0;
	double LEDfreq=0, LEDperiod=0;
	bool badLEDFreq=false;
	bool LEDTurnedOff = false;
	int simpleLEDCount=0;
	bool useSimpleThreshold=false;

	// Error check variables
	const int nErrs=29; // error 0 is unused
  vector<int> SeriousErrors = {1, 13, 14, 18, 19, 20, 21, 22, 23, 24, 25, 26};
	int SeriousErrorCount = 0;
	int TotalErrorCount = 0;
	vector<int> Error(nErrs); // write this to ROOT tree
	vector<int> ErrorCount(nErrs); // don't write this, but keep it for the error summary.
	long skippedEvents=0;

	// time variables
  double xTime=0;
	double deltaScaler=0, deltaSBC=0;
	double firstGoodScaler=0, lastGoodScaler=0;
  double timePrevLED=0;
	long start=0, stop=0;
	double unixDuration=0, scalerDuration=0;
  int firstEntryAfterFlush=0;
  double jumpCorrection=0;
  double scalerOffset=0,scalerUnc=0;

	// muon ID variables
	bool LEDCut = true;
	bool EnergyCut = false;
	vector<int> CoinType(32), Plane(32);

	// initialize input data
	long vEntries = vetoChain->GetEntries();
	TTreeReader reader(vetoChain);
	TTreeReaderValue<unsigned int> vMult(reader, "mVeto");
	TTreeReaderValue<uint32_t> vBits(reader, "vetoBits");
	TTreeReaderValue<MGTBasicEvent> vEvt(reader,"vetoEvent");
	TTreeReaderValue<MJTRun> vRun(reader,"run");
  TTreeReaderValue<long> bTimeStart(reader,"fStartTime");
  TTreeReaderValue<long> bTimeStop(reader,"fStopTime");
	reader.SetEntry(0);
	int runNum = vRun->GetRunNumber();
	// start = vRun->GetStartTime(); // method 1
	// stop = vRun->GetStopTime();
  start = (*bTimeStart);  // method 2
  stop = (*bTimeStop);
	unixDuration = (double)(stop - start);
	reader.SetTree(vetoChain);  // resets the reader

	// MJVetoEvent variables, with run-based card numbers
	int card1=0, card2=0;
	SetCardNumbers(runNum,card1,card2);
	MJVetoEvent veto(card1,card2);
	MJVetoEvent sync, prev, prevLED, out;

	// initialize output file
	char outputFile[200];
	sprintf(outputFile,"%s/veto_run%i.root",outputDir.c_str(),runNum);
	TFile *RootFile = new TFile(outputFile, "RECREATE");
	TTree *vetoTree = new TTree("vetoTree","MJD Veto Events");
	// event info
	vetoTree->Branch("run",&runNum);
	vetoTree->Branch("vetoEvent","MJVetoEvent",&out,32000,1);
	// time variables
  vetoTree->Branch("xTime",&xTime);
  vetoTree->Branch("jumpCorrection",&jumpCorrection);
	vetoTree->Branch("deltaScaler",&deltaScaler);
	vetoTree->Branch("deltaSBC",&deltaSBC);
	vetoTree->Branch("timePrevLED",&timePrevLED);
	vetoTree->Branch("start",&start,"start/L");
	vetoTree->Branch("stop",&stop,"stop/L");
	vetoTree->Branch("unixDuration",&unixDuration);
  vetoTree->Branch("scalerDuration",&scalerDuration);
  vetoTree->Branch("scalerOffset",&scalerOffset);
  vetoTree->Branch("scalerUnc",&scalerUnc);
  vetoTree->Branch("firstEntryAfterFlush",&firstEntryAfterFlush);
	// LED variables
	vetoTree->Branch("LEDfreq",&LEDfreq);
	vetoTree->Branch("multipThreshold",&multipThreshold);
	vetoTree->Branch("highestMultip",&highestMultip);
	vetoTree->Branch("LEDMultipThreshold",&LEDMultipThreshold);
	vetoTree->Branch("LEDSimpleThreshold",&LEDSimpleThreshold);
	vetoTree->Branch("useSimpleThreshold",&useSimpleThreshold);
	// muon ID variables
	vetoTree->Branch("CoinType",&CoinType);
	vetoTree->Branch("Plane",&Plane);
	// error variables
	vetoTree->Branch("Errors",&Error);

	// Error "garbage event" tree
	TTree *skipTree = new TTree("skipTree","skipped veto events");
	skipTree->Branch("run",&runNum);
	skipTree->Branch("vetoEvent","MJVetoEvent",&out,32000,1);
	skipTree->Branch("Errors",&Error);
	skipTree->Branch("start",&start,"start/L");
	skipTree->Branch("stop",&stop,"stop/L");

	// ==================== 1st loop over veto entries  =================
  // Measure the LED frequency, find the highest-multiplicity entry,
  // identify buffer flush bursts (so we can ignore any scaler jumps
  // during the flush because deltaSBC is not trustworthy),
  // and compare the unix duration with the scaler duration.

  bool foundSyncEvent = false;
  bool foundBufferFlush = false;
	TH1D *LEDDeltaT = new TH1D("LEDDeltaT","LEDDeltaT",100000,0,100); // 0.001 sec/bin
	while(reader.Next())
	{
		long i = reader.GetCurrentEntry();
		veto.Clear();
		veto.SetSWThresh(swThresh);
		veto.WriteEvent(i,&*vRun,&*vEvt,*vBits,runNum,true);

    if (firstGoodScaler==0 && !veto.GetBadScaler())
      firstGoodScaler = veto.GetTimeSec();
    if (!veto.GetBadScaler()) lastGoodScaler = veto.GetTimeSec();

		if (CheckErrors(veto,prev,Error)){
			skippedEvents++;
      if (Error[25]) {
        foundBufferFlush = true;
        firstEntryAfterFlush = i;
        // cout << i << " Found buffer flush.  Index: " << veto.GetScalerIndex() << endl;
      }
			// do end of loop reset
			prev = veto;
			continue;
		}
    // Don't let the sync event be the first one, we need a Ge entry before and after it.
    else if (!foundSyncEvent && i > 1) {
      foundSyncEvent = true;
      sync = veto;
    }
		if (veto.GetMultip() > highestMultip)
      highestMultip = veto.GetMultip();

		if (veto.GetMultip() > LEDSimpleThreshold) {
			LEDDeltaT->Fill(veto.GetTimeSec()-prev.GetTimeSec());
			simpleLEDCount++;
		}
		// end of loop reset
		prev = veto;
	}
  if (foundBufferFlush) {
    firstEntryAfterFlush += 1;
    cout << "Warning: found buffer flush.  First entry after : " << firstEntryAfterFlush << endl;
    reader.SetEntry(firstEntryAfterFlush);
    sync.WriteEvent(firstEntryAfterFlush,&*vRun,&*vEvt,*vBits,runNum,true);
  }

  // ============== Loop 1-a: Scan built data for rough sync ==============

  // Scan the built data for this run to determine if there is a scaler offset.
  // Find times of Ge events whose packets are immediately before and after the "sync" event.
  bool applyOffset = false;
  if (foundSyncEvent)
  {
    GATDataSet *ds = new GATDataSet(runNum);
    TChain *builtChain = ds->GetBuiltChain(false);
    MGTEvent *evt=0;
    MGVDigitizerData *dig=0;
    builtChain->SetBranchAddress("event",&evt);
    double bTimeFirst=0, bTimeBefore=0, bTimeAfter=0;
    uint64_t bItr=0,bIndex=0;
    while (true)
    {
      builtChain->GetEntry(bItr);
      if (evt->GetNDigitizerData() > 0)
      {
        dig = evt->GetDigitizerData(0);
        bIndex = dig->GetIndex();
        if (bTimeFirst == 0) bTimeFirst = ((double)dig->GetTimeStamp())*1.e-8;
        if ((int)bIndex < sync.GetScalerIndex())
          bTimeBefore = ((double)dig->GetTimeStamp())*1.e-8;
        if ((int)bIndex > sync.GetScalerIndex())
          bTimeAfter = ((double)dig->GetTimeStamp())*1.e-8;
        // printf("%li  ind %lu  first %.3f  before %.3f  after %.3f  diff: %i\n", bItr,bIndex,bTimeFirst,bTimeBefore,bTimeAfter,(int)bIndex - sync.GetScalerIndex());
      }
      if ((int)bIndex > sync.GetScalerIndex()) break;
      if ((int)bItr > builtChain->GetEntries()-1) break;
      bItr++;
    }
    delete ds;
    double bVetoTime = (bTimeAfter + bTimeBefore)/2.;
    scalerOffset = bVetoTime-sync.GetTimeSec();
    scalerUnc = (bTimeAfter - bTimeBefore)/2.;
    printf("Syncing entry %i (packet %li) with Ge timestamps.\n",sync.GetEntry(),sync.GetScalerIndex());
    if (fabs(sync.GetTimeSec() - bVetoTime) > scalerUnc)
    {
      printf("Scaler (%.2f) out of sync with trigger card (%.2f) by %.3f +/- %.3f sec.\n", sync.GetTimeSec(),bVetoTime,scalerOffset,scalerUnc);
      applyOffset = true;
    }
    if (scalerUnc == bVetoTime) {
      cout << "Warning: Sync failed.\n";
      applyOffset=false;
    }
  }
  else cout << "Unable to sync veto and Ge clocks.\n";
  if (!applyOffset){
    scalerOffset=0;
    scalerUnc=0;
  }

  // =======================================================================
  cout << "===================== Veto Error Report =====================\n";

	// Determine run duration from start and stop packets
  scalerDuration = lastGoodScaler - firstGoodScaler;
  if (start == 0 || stop == 0)
  {
    cout << "Warning: missing start or stop packet.  Start: " << start << "  Stop: " << stop
         << "\n  Replacing unix duration with scaler duration.  First Scaler: " << firstGoodScaler << "  Last Scaler: " << lastGoodScaler
         << "\n  Setting unix duration to " << lastGoodScaler - firstGoodScaler
         << "\n  NOTE: scalerDuration - 3600 = " << scalerDuration - 3600 << "  (large excess indicates buffer flush problems)\n";
    unixDuration = scalerDuration;
  }

	// Error 27: QDC threshold not found
	// Error 28: No events above QDC threshold
	for (int i=0; i < 32; i++) if (swThresh[i] == 9999) {
		Error[27]=true;
		ErrorCount[27]++;
		cout << "Warning: Couldn't find QDC threshold for panel " << i << ". Set to 9999\n";
	}

	// Set LED multiplicity threshold, find LED frequency, and use alternate method if we have a short run.
	multipThreshold = highestMultip - LEDMultipThreshold;
	if (multipThreshold < 0) multipThreshold = 0;
	int dtEntries = LEDDeltaT->GetEntries();
	if (dtEntries > 0) {
		int maxbin = LEDDeltaT->GetMaximumBin();
		LEDDeltaT->GetXaxis()->SetRange(maxbin-100,maxbin+100); // looks at +/- 0.1 seconds of max bin.
		LEDfreq = 1/LEDDeltaT->GetMean();
	}
	else {
		cout << "Warning! No multiplicity > " << LEDSimpleThreshold << " events.  LED may be off.\n";
		LEDfreq = 9999;
		badLEDFreq = true;
	}
	LEDperiod = 1/LEDfreq;
	delete LEDDeltaT;
	if (LEDperiod > 9 || vEntries < 100) {
		cout << "Warning: Short run.\n";
		if (simpleLEDCount > 3) {
			cout << "  From delta-T histogram, LED frequency is " << LEDfreq << " Hz."
				  << "\n  Reverting to 'simple' rate: " <<  simpleLEDCount/unixDuration << " Hz.\n";
			LEDperiod = unixDuration/simpleLEDCount;
			useSimpleThreshold=true;
		}
		else {
			LEDperiod = 9999;
			badLEDFreq = true;
		}
	}
	// Error 26: LED frequency very low/high, corrupted, or LED's off.
	if (LEDperiod > 20 || LEDperiod < 0 || badLEDFreq) {
		ErrorCount[26]++;
		Error[26] = true;
	}

  // ========================================================================
	// ================ 2nd loop over entries - Error checks ==================
	// We don't skip any events, and we count the number of each type of error.

	reader.SetTree(vetoChain); // reset the reader
	std::fill(Error.begin(), Error.end(), 0); // reset error bools
	while(reader.Next())
	{
		long i = reader.GetCurrentEntry();
		veto.Clear();
		veto.SetSWThresh(swThresh);
		veto.WriteEvent(i,&*vRun,&*vEvt,*vBits,runNum,true);
		CheckErrors(veto,prev,Error);
		for (int j=0; j<nErrs; j++) if (Error[j]==1) ErrorCount[j]++;

		// Print errors to screen
		bool PrintError = false;
		for (auto i : SeriousErrors){
			if (Error[i]) {
				PrintError = true;
				break;
			}
		}
		if (PrintError)
		{
      if (Error[1] && Error[25])  cout << i << ":[1] Missing Packet & [25] Buffer Flush.";
      if (Error[1] && !Error[25]) cout << i << ":[1] Missing Packet.";
      if (!Error[1] && Error[25]) cout << i << ":[25] Buffer Flush.";
      if (Error[1] || Error[25])
        printf("  Index %li  Scaler %-5.2f  d(sca) %-5.3f  d(sbc) %-5.3f\n", veto.GetScalerIndex(),veto.GetTimeSec(),veto.GetTimeSec()-prev.GetTimeSec(),veto.GetTimeSBC()-prev.GetTimeSBC());

      if (Error[13])
				cout << i << ":[13] Indexes of QDC1 and Scaler differ by more than 2."
					 << "\n    Scaler Index " << veto.GetScalerIndex()
					 << "  QDC1 Index " << veto.GetQDC1Index()
					 << "\n    Previous scaler Index " << prev.GetScalerIndex()
					 << "  Previous QDC1 Index " << prev.GetQDC1Index() << endl;

      if (Error[14])
				cout << i << ":[14] Indexes of QDC2 and Scaler differ by more than 2."
					 << "\n    Scaler Index " << veto.GetScalerIndex()
					 << "  QDC2 Index " << veto.GetQDC2Index()
					 << "\n    Previous scaler Index " << prev.GetScalerIndex()
					 << "  Previous QDC2 Index " << prev.GetQDC2Index() << endl;

      if (Error[18])
				cout << i << ":[18] Scaler/SBC Desynch."
					  << "\n    Scaler " << veto.GetTimeSec() << "  SBC " << veto.GetTimeSBC()
					  << "\n    Delta(scaler) " << veto.GetTimeSec() - prev.GetTimeSec()
					  << "\n    Delta(sbc) " << veto.GetTimeSBC() - prev.GetTimeSBC()
            << "\n    Scaler jump correction: " << (veto.GetTimeSBC()-prev.GetTimeSBC()) - (veto.GetTimeSec()-prev.GetTimeSec())
					  << "\n    Adjusted time: "
					  << veto.GetTimeSec() + (veto.GetTimeSBC()-prev.GetTimeSBC()) - (veto.GetTimeSec()-prev.GetTimeSec()) << endl;

      if (Error[19])
				cout << i << ":[19] Scaler Event Count Reset. "
					 << "\n    Scaler Index " << veto.GetScalerIndex()
					 << "  SEC " << veto.GetSEC()
					 << "  Previous SEC " << prev.GetSEC() << "\n";

      if (Error[20])
				cout << i << ":[20] Scaler Event Count Jump."
					 << "\n    Scaler Time " << veto.GetTimeSec()
					 << "  Scaler Index " << veto.GetScalerIndex()
					 << "  Prev scaler time " << prev.GetTimeSec()
					 << "\n    SEC " << veto.GetSEC()
					 << "  Previous SEC " << prev.GetSEC() << "\n";

      if (Error[21])
				cout << i << ":[21] QDC1 Event Count Reset."
					 << "\n    Scaler Index " << veto.GetScalerIndex()
					 << "  QEC1 " << veto.GetQEC()
					 << "  Previous QEC1 " << prev.GetQEC() << "\n";
			if(Error[22])
				cout << i << ":[22] QDC 1 Event Count Jump."
					 << "\n    Scaler time " << veto.GetTimeSec()
					 << "  QDC 1 Index " << veto.GetQDC1Index()
					 << "  QEC 1 " << veto.GetQEC()
					 << "  Previous QEC 1 " << prev.GetQEC() << "\n";

      if (Error[23])
				cout << i << ":[23] QDC2 Event Count Reset."
					 << "\n    Scaler Index " << veto.GetScalerIndex()
					 << "  QEC2 " << veto.GetQEC2()
					 << "  Previous QEC2 " << prev.GetQEC2() << "\n";

      if(Error[24])
				cout << i << ":[24] QDC 2 Event Count Jump."
					 << "\n    Scaler time " << veto.GetTimeSec()
					 << "  QDC 2 Index " << veto.GetQDC2Index()
					 << "  QEC 2 " << veto.GetQEC2()
					 << "  Previous QEC 2 " << prev.GetQEC2() << "\n";
		}
		// end of event resets
		prev = veto;
		std::fill(Error.begin(), Error.end(), 0);
	}
	// Calculate total errors and total serious errors
	// Ignore Error 10 & 11 - the veto counters are not reset at the beginning of runs.
	for (int i = 1; i < nErrs; i++) {
			if (i != 10 && i != 11) TotalErrorCount += ErrorCount[i];
			for (auto j : SeriousErrors) if (i == j) SeriousErrorCount += ErrorCount[i];
	}
	cout << "Serious errors found :: " << SeriousErrorCount << endl;
	if (SeriousErrorCount > 0)
	{
		// cout << "Total Errors : " << TotalErrorCount << endl;
		for (int i = 1; i < nErrs; i++)
		{
			if (ErrorCount[i] > 0 && (i!=7 && i!=10 && i!=11))
			{
				if (i != 26) cout << "  Error[" << i <<"]: " << ErrorCount[i] << " events ("<< 100*(double)ErrorCount[i]/vEntries << " %)\n";
		 		else if (i == 26) {
					cout << "  Error[26]: Bad LED rate: " << LEDfreq << "  Period: " << LEDperiod << endl;
					if (LEDperiod > 0.1 && (abs(unixDuration/LEDperiod) - simpleLEDCount) > 5) {
						cout << "   Simple LED count: " << simpleLEDCount
							 << "  Expected: " << (int)(unixDuration/LEDperiod) << endl;
					}
				}
			}
		}
		// cout << "For reference, \"serious\" error types are: ";
		// for (auto i : SeriousErrors) cout << i << " ";
		// cout << "\nPlease report these to the veto group.\n";
	}
	if (errorCheckOnly) return;

	// ================ 3nd loop over entries - Find muons! =================
	// Determine event time, skip bad entries, and apply all cuts for muon ID.

  cout << "=================== Scanning for muons ... ==================\n";

	reader.SetTree(vetoChain); // reset the reader
	prev.Clear();
	skippedEvents = 0;
	printf("unixDuration %.0f sec  Highest mult. %i  LED threshold %i\n", unixDuration,highestMultip,multipThreshold);
	while(reader.Next())
	{
		long i = reader.GetCurrentEntry();
		veto.Clear();
		veto.SetSWThresh(swThresh);
		veto.WriteEvent(i,&*vRun,&*vEvt,*vBits,runNum,true);
		CheckErrors(veto,prev,Error);
		for (int j=0; j<nErrs; j++) if (Error[j]==1) ErrorCount[j]++;

    deltaScaler = veto.GetTimeSec()-prev.GetTimeSec();
    deltaSBC = veto.GetTimeSBC()-prev.GetTimeSBC();

    // Apply the results from the veto-ge sync to xTime
    xTime = veto.GetTimeSec();
    if (applyOffset) xTime += scalerOffset;

    // Scaler jump handling:
    // Calculate the jumpCorrection and save it to the ROOT output.
    // Ignore any scaler jumps that happen during a buffer flush,
    // because deltaSBC is not trustworthy.
    if (i > firstEntryAfterFlush && Error[18]) {
      jumpCorrection += deltaSBC - deltaScaler;
      cout << "Scaler jump found.  Applying jump correction: " << jumpCorrection << endl;
      xTime += jumpCorrection;
    }
		// debug block (don't delete!)
    // if (i > 715 && i < 720)
    // printf("%li  ind %li  e1 %i  e18 %i  e19 %i  scaler %-5.2f  dScaler %-5.2f  dSBC %-5.2f  jumpCor %-5.2f\n" ,i,veto.GetScalerIndex(),Error[1],Error[18],Error[19],veto.GetTimeSec(),deltaScaler,deltaSBC,jumpCorrection);

    // Skip bad events and fill the skipTree.
		if (CheckErrors(veto,prev,Error))
		{
			skippedEvents++;
			// do the end-of-event reset
			prev = veto;
			skipTree->Fill();
			veto.Clear();
			continue;
		}

		// LED Cut
		LEDCut = false;
		LEDTurnedOff = ErrorCount[26];
		if (veto.GetMultip() < multipThreshold && !LEDTurnedOff) LEDCut = true;
		else if (LEDTurnedOff) LEDCut = true;

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

    // debug block (don't delete!)
		// if (veto.GetMultip() < 27 && veto.GetMultip() > 1)
		// printf("Entry %li  Time %-6.2f  QDC %-5i  Mult %i  Ov500 %i  Loff? %i  LEDCut %i  ECut %i  \n",i,veto.GetTimeSec(),veto.GetTotE(),veto.GetMultip(),over500Count,LEDTurnedOff,LEDCut,EnergyCut);

		// Muon Identification:
		// Use EnergyCut, LEDCut, and the Hit Pattern to identify them sumbitches.
		for (int k = 0; k < 12; k++) Plane[k] = 0;
		for (int k = 0; k < 32; k++) {
			if (veto.GetQDC(k) > veto.GetSWThresh(k)) {
				if (PlaneMap(k,runNum)==0)       Plane[0]=1;  // 0: Lower Bottom
				else if (PlaneMap(k,runNum)==1)  Plane[1]=1;  // 1: Upper Bottom
				else if (PlaneMap(k,runNum)==2)  Plane[2]=1;  // 3: Top Inner
				else if (PlaneMap(k,runNum)==3)  Plane[3]=1;  // 4: Top Outer
				else if (PlaneMap(k,runNum)==4)  Plane[4]=1;  // 5: North Inner
				else if (PlaneMap(k,runNum)==5)  Plane[5]=1;  // 6: North Outer
				else if (PlaneMap(k,runNum)==6)  Plane[6]=1;  // 7: South Inner
				else if (PlaneMap(k,runNum)==7)  Plane[7]=1;  // 8: South Outer
				else if (PlaneMap(k,runNum)==8)  Plane[8]=1;  // 9: West Inner
				else if (PlaneMap(k,runNum)==9)  Plane[9]=1;  // 10: West Outer
				else if (PlaneMap(k,runNum)==10) Plane[10]=1; // 11: East Inner
				else if (PlaneMap(k,runNum)==11) Plane[11]=1; // 12: East Outer
				else if (PlaneMap(k,runNum)==-1)
					cout << "Error: Panel " << k << " was not installed for this run and should not be giving counts above threshold.\n";
			}
		}
    std::fill(CoinType.begin(), CoinType.end(), 0);  // reset
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
			// for (int i=0; i<12; i++) cout << "p" << i << ":" << Plane[i] << "  ";
			// cout << endl;

			// Type 1: vertical muons
			if (Plane[0] && Plane[1] && Plane[2] && Plane[3]) {
				CoinType[1]=true;
				a=true;
				type=1;
			}
			// Type 2: (both planes of a side) + both bottom planes
			if ((Plane[0] && Plane[1]) &&
				  ((Plane[4] && Plane[5])
          || (Plane[6] && Plane[7])
				  || (Plane[8] && Plane[9])
				  || (Plane[10] && Plane[11])))
			{
				CoinType[2] = true;
				b=true;
				type = 2;
			}
			// Type 3: (both top planes) + (both planes of a side)
			if ((Plane[2] && Plane[3]) &&
          ((Plane[4] && Plane[5])
          || (Plane[6] && Plane[7])
					|| (Plane[8] && Plane[9])
          || (Plane[10] && Plane[11]))) {
				CoinType[3] = true;
				c=true;
				type = 3;
			}
			// Type 4: compound hit (combination of types 1-3)
			if ((a && b)||(a && c)||(b && c)) type = 4;

			char hitType[200];
			if (type==0) sprintf(hitType,"2+ panels");
			if (type==1) sprintf(hitType,"vertical");
			if (type==2) sprintf(hitType,"side+bottom");
			if (type==3) sprintf(hitType,"top+sides");
			if (type==4) sprintf(hitType,"compound");

			// print the details of the hit
			printf("Hit: %-12s Entry %-4li Time %-6.2f  QDC %-5i  Mult %i  Ov500 %i  LEDoff %i\n", hitType,i,xTime,veto.GetTotE(),veto.GetMultip(),over500Count,LEDTurnedOff);
		}

		out = veto;
		vetoTree->Fill();
		// end of event resets
		prev = veto;
		if (veto.GetMultip() > multipThreshold) {
      if (!veto.GetBadScaler()) timePrevLED = veto.GetTimeSec();
      else timePrevLED = -1;
    }
	}
	if (skippedEvents > 0) printf("ProcessVetoData skipped %li of %li entries.\n",skippedEvents,vEntries);

  vetoTree->Write("",TObject::kOverwrite);
	skipTree->Write("",TObject::kOverwrite);
  cout << "Wrote ROOT file: " << outputFile << endl;

	RootFile->Close();
}

// ====================================================================================
// =================================VETO TOOL KIT======================================
// ====================================================================================

void SetCardNumbers(int runNum, int &card1, int &card2)
{
	if (runNum > 45000000)
		{ card1 = 11;  card2 = 18; }
	else
		{ card1 = 13;  card2 = 18; }
}

int FindThreshold(TH1D *qdcHist, int threshVal, int panel, int runNum)
{
	// Returns 9999 if a panels is deactivated or the threshold is not found.
	// This (intentionally) causes that panel to not contribute to multiplicity or total QDC.
	// This is the run-level error 27/28, and is checked between loop 1 and loop 2.

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
bool CheckErrors(MJVetoEvent veto, MJVetoEvent prev)
{
	vector<int> ErrorVec(29);
	std::fill(ErrorVec.begin(), ErrorVec.end(), 0);
	return CheckErrors(veto,prev,ErrorVec);
}
bool CheckErrors(MJVetoEvent veto, MJVetoEvent prev, vector<int> &ErrorVec)
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
	// s 25. Buffer flush error.

	// Run-level error checks (not checked in this function)
	// 26. LED frequency very low/high, corrupted, or LED's off.
	// 27. QDC threshold not found
	// 28. No events above QDC threshold
	*/

	vector<int> Errors(29);
	std::fill(Errors.begin(), Errors.end(), 0);

	// Errors 1-18 are checked automatically when we call MJVetoEvent::WriteEvent
	for (int i=0; i < 18; i++) Errors[i] = veto.GetError(i);

	for (int q=0; q<18; q++) {
		if (Errors[q]==1 && (q==1||q==2||q==3||q==5||q==6||q==9||q==13||q==14))
			skip = true;
	}

  bool foundBothQDC = (!veto.GetError(1) && !prev.GetError(1));

	if (foundBothQDC && veto.GetEntry() > 1 && veto.GetTimeSec() > 0 && veto.GetTimeSBC() > 0
			&& fabs((veto.GetTimeSec() - prev.GetTimeSec())-(veto.GetTimeSBC() - prev.GetTimeSBC())) > 1) {
				Errors[18] = true;
			}

	if (!veto.GetError(1) && veto.GetSEC() == 0 && veto.GetEntry() > 1)
		Errors[19] = true;

	if (foundBothQDC && veto.GetEntry() > 1 && abs(veto.GetSEC() - prev.GetSEC()) > veto.GetEntry()-prev.GetEntry() && veto.GetSEC()!=0)
		Errors[20] = true;

	if (!veto.GetError(1) && veto.GetQEC() == 0 && veto.GetEntry() >1)
		Errors[21] = true;

	if (foundBothQDC && veto.GetEntry() > 1 && abs(veto.GetQEC() - prev.GetQEC()) > veto.GetEntry()-prev.GetEntry() && veto.GetQEC() != 0)
		Errors[22] = true;

	if (!veto.GetError(1) && veto.GetQEC2() == 0 && veto.GetEntry() > 1)
		Errors[23] = true;

	if (foundBothQDC && abs(veto.GetQEC2() - prev.GetQEC2()) > veto.GetEntry()-prev.GetEntry() && veto.GetEntry() > 1 && veto.GetQEC2() != 0)
		Errors[24] = true;

  if (abs(veto.GetScalerIndex() - prev.GetScalerIndex()) == 1)
    Errors[25] = true;

	// Check errors 18-25 (don't accidentally change the run-level errors 26-28)
	for (int q=18; q<26; q++) {
		if (Errors[q]==1 && (q==18||q==19||q==20||q==21||q==22||q==23||q==24||q==25))
			skip = true;
	}

  // assign Errors to the output vetor
	ErrorVec = Errors;

  // debug: print error vector
  // for (int i=0; i < 28; i++) cout << i << ":" << Errors[i] << "  ";
	// cout << endl;

	return skip;
}
