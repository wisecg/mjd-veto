// skim-veto.cc
// Takes output of auto-veto.cc for veto analysis.
// C. Wiseman, 10/24/16

#include <iostream>
#include <fstream>
#include <numeric>
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include "MJVetoEvent.hh"

#include "MGVDigitizerData.hh"
#include "MGTEvent.hh"
#include "GATDataSet.hh"

using namespace std;

void GenerateVetoList(TChain *vetoTree);
void ScanBuiltData();
void ListRunOffsets(TChain *vetoTree);
void GenerateDisplayList(TChain *vetoTree);
void CalculateDeadTime(string MuonList, int dsNumber);

int PanelMap(int i, int runNum);

int main(int argc, char** argv)
{
	if (argc < 2) {
		cout << "Usage: ./skim-veto -r [run number]\n"
         << "       ./skim-veto [run list file]\n";
		return 0;
	}
  TChain *vetoTree = new TChain("vetoTree");

  int run = 0;
  char file[200];
  string opt1 = argv[1];
  if (opt1 == "-r") {
    run = stoi(argv[2]);
    sprintf(file,"./avout/veto_run%i.root",run);
    if (!vetoTree->Add(file)){
      cout << "File doesn't exist.  Exiting ... \n";
      return 1;
    }
  }
  else {
    ifstream runFile(argv[1]);
    // make sure there are no duplicate runs,
    // and all the run numbers are in increasing order.
    set<int> uniqueRuns;
    while (runFile >> run) uniqueRuns.insert(run);
    vector<int> runList(uniqueRuns.begin(), uniqueRuns.end());
    sort(runList.begin(), runList.end());

    for (auto i : runList) {
      sprintf(file,"./avout/veto_run%i.root",i);
      if (!vetoTree->Add(file)){
        cout << "File doesn't exist.  Continuing ... \n";
        continue;
      }
    }
  }
	printf("Found %lli entries.\n",vetoTree->GetEntries());

  ListRunOffsets(vetoTree);
	// GenerateVetoList(vetoTree);
	// GenerateDisplayList(vetoTree);
	// CalculateDeadTime("./output/MuonList_test.txt",1);
}

void ListRunOffsets(TChain *vetoTree)
{
  TH1D *hUnc = new TH1D("hUnc","hUnc",100,0,0.2);

  TTreeReader reader(vetoTree);
  TTreeReaderValue<int> runIn(reader,"run");
  TTreeReaderValue<MJVetoEvent> vetoEventIn(reader,"vetoEvent");
  TTreeReaderValue<double> xTimeIn(reader,"xTime");
  TTreeReaderValue<double> scalerOffsetIn(reader,"scalerOffset");
  TTreeReaderValue<double> scalerUncIn(reader,"scalerUnc");
  int runSave = 0;
  while(reader.Next())
  {
    if(runSave != *runIn) // run boundary condition
    {
      runSave = *runIn;
      MJVetoEvent veto = *vetoEventIn;
      printf("Run %i  Entry %i  xTime %-5.3f  Offset %-5.3f  Unc %-5.3f\n", veto.GetRun(),veto.GetEntry(),*xTimeIn,*scalerOffsetIn,*scalerUncIn);

      hUnc->Fill(*scalerUncIn);
    }
  }
  TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",800,600);
  hUnc->Draw();
  // hUnc->Fit("gaus");  // not really gaussian, probably flat with a big spike at 0
  c1->Print("./output/scalerUnc.pdf");
}

void GenerateVetoList(TChain *vetoTree)
{
	ofstream MuonList("./output/MuonList_test.txt");

	TTreeReader reader(vetoTree);
	TTreeReaderValue<MJVetoEvent> vetoEventIn(reader,"vetoEvent");
	TTreeReaderValue<Long64_t> start(reader,"start");
	TTreeReaderValue<Long64_t> stop(reader,"stop");
  TTreeReaderValue<double> jumpCorrection(reader,"jumpCorrection");
	TTreeReaderArray<int> CoinType(reader,"CoinType");	//[32]

	bool newRun=false;
	int prevRun=0;
	Long64_t prevStop=0;

	while(reader.Next())
	{
		MJVetoEvent veto = *vetoEventIn;
		int run = veto.GetRun();
		if (run != prevRun) newRun=true;
		else newRun = false;

		int type = 0;
		if (CoinType[0]) type=1;
		if (CoinType[1]) type=2;	// overrides type 1 if both are true
		if ((*start-prevStop) > 10 && newRun) type = 3;

		char muonList[200];
		if (type > 0) {
			if (type==1 || type==2)
				sprintf(muonList,"%i %lli %.3f %.3f%i %i\n",run,*start,veto.GetTimeSec(),*jumpCorrection,type,veto.GetBadScaler());
			else if (type==3)
				sprintf(muonList,"%i %lli 0.0 3 0\n",run,*start);
			cout << muonList;
			MuonList << muonList;
		}

		// end of entry, save the run and stop time
		prevStop = *stop;
		prevRun = run;
	}
	MuonList.close();
}

void ScanBuiltData()
{
  MJVetoEvent first;
  int runNum = 16991;
  double scalerOffset=0, scalerUnc=0; // deprecated ...

  cout << "Finding scaler offset ...\n";
  GATDataSet *ds = new GATDataSet(runNum);
  TChain *builtChain = ds->GetBuiltChain(false);
  MGTEvent *evt=0;
  MGVDigitizerData *dig=0;
  builtChain->SetBranchAddress("event",&evt);
  double bTimeFirst=0, bTimeBefore=0, bTimeAfter=0;
  uint64_t bItr=0,bIndex=0;
  while (true) {
    builtChain->GetEntry(bItr);
    // CAUTION: Are you sure the built data doesn't have events from the previous run at the beginning?
    if (evt->GetNDigitizerData() > 0)
    {
      dig = evt->GetDigitizerData(0);
      bIndex = dig->GetIndex();
      if (bTimeFirst == 0)
        bTimeFirst = ((double)dig->GetTimeStamp())*1.e-8;
      if ((int)bIndex < first.GetScalerIndex())
        bTimeBefore = ((double)dig->GetTimeStamp())*1.e-8;
      bTimeAfter = ((double)dig->GetTimeStamp())*1.e-8;

      // printf("%li  ind %lu  first %.1f  before %.1f  after %.1f  time-first %.1f\n", bItr,bIndex,bTimeFirst,bTimeBefore,bTimeAfter,bTimeBefore-bTimeFirst,bTimeAfter-bTimeFirst,bTimeAfter-bTimeFirst);
    }
    if ((int)bIndex > first.GetScalerIndex()) break;
    if ((int)bIndex > builtChain->GetEntries()) break;
    bItr++;
  }
  scalerOffset = (bTimeAfter + bTimeBefore)/2 - bTimeFirst;
  scalerUnc = (bTimeAfter - bTimeBefore)/2;
  // printf("First veto event, entry %li, index %li.  Ge times: start %.2f,  before %.2f,  after, %.2f.\n" ,first.GetEntry(),first.GetScalerIndex(),bTimeFirst,bTimeBefore,bTimeAfter);
  delete ds;
  printf("Scaler offset: %.2f +/- %.2f seconds.\n",scalerOffset,scalerUnc);

}

// this function used to be in auto-veto but I scrubbed it.  Didn't want to waste it completely, so it's here if needed.
int FindEventTime(double &xTime, double &sbcTime, MJVetoEvent &veto, MJVetoEvent &first,
	vector<double> &RawScalerTimes, vector<bool> &BadScalers, double duration, long vEntries)
{
	// returns 0: success  1: sbc time  2: interp time  3: entry time
	int timeMethod = 0;
	bool goodSBC = (veto.GetRun() > 8557 && veto.GetRun() < 4500000 && veto.GetTimeSBC() < 2000000000);
	bool vecSizes = (RawScalerTimes.size() == BadScalers.size());

	if (!veto.GetBadScaler())  // 0. best case (good scaler, good sbc)
	{
		xTime = veto.GetTimeSec() - first.GetTimeSec();
		sbcTime = veto.GetTimeSBC() - first.GetTimeSBC();
		timeMethod = 0;
	}
	else if (veto.GetBadScaler() && goodSBC)  // 1. next-best case (bad scaler but good sbc)
	{
		sbcTime = veto.GetTimeSBC() - first.GetTimeSBC();
		xTime = sbcTime;
		timeMethod = 1;
	}
	else if (veto.GetBadScaler() && !goodSBC && vecSizes)  // 2. try interpolating time
	{
		// TODO: if the scalers are mostly corrupted in a run, it may be better to put in a flag
		// that reverts to the entry time method.
		sbcTime=0;
		double loTime=0, upTime=0;
		if (!BadScalers[veto.GetEntry()]) xTime = RawScalerTimes[veto.GetEntry()];
		else
		{
			for (int i = veto.GetEntry(); i < (int)RawScalerTimes.size(); i++)
				if (BadScalers[i] == 0) { upTime = RawScalerTimes[i]; break; }
			for (int j = veto.GetEntry(); j > 0; j--)
				if (BadScalers[j] == 0) { loTime = RawScalerTimes[j]; break; }
			xTime = (upTime + loTime)/2.0;
		}
		timeMethod = 2;
	}
	else	// 3. last resort: use "entry" time
	{
		sbcTime=0;
		xTime = ((double)veto.GetEntry() / vEntries) * duration; // this breaks if we have corrupted duration
		timeMethod = 3;
	}
	return timeMethod;
}


void GenerateDisplayList(TChain *vetoTree)
{
	// Format:  (run) (unix start time) (time within run) (entry) (type) qdc1 ... qdc32
	ofstream DisplayList("./output/MuonDisplay_test.txt");

	TTreeReader reader(vetoTree);
	TTreeReaderValue<MJVetoEvent> events(reader,"events");
	TTreeReaderValue<Long64_t> start(reader,"start");
	TTreeReaderValue<double> xTime(reader,"xTime");
	TTreeReaderArray<int> CoinType(reader,"CoinType");	//[32]

	while(reader.Next())
	{
		long i = reader.GetCurrentEntry();
		MJVetoEvent veto = *events;

		int type = 0;
		if (CoinType[0]) type=1;
		if (CoinType[1]) type=2;	// overrides type 1 if both are true

		char display[200];
		if (type>0) // all events passing TimeCut & EnergyCut
		{
			sprintf(display,"%i  %li  %lli  %.3f  ",events->GetRun(),i,*start,*xTime);
			DisplayList << display;
			for (int j=0; j<32; j++)
			{
				if (events->GetQDC(j) >= events->GetSWThresh(j))
					DisplayList << events->GetQDC(j) << " ";
				else
					DisplayList << 0 << " ";
			}
			DisplayList << endl;
		}
	}
}

void CalculateDeadTime(string MuonList, int dsNumber)
{
	ifstream InputList(MuonList.c_str());
	if(!InputList.good()) {
    	cout << "Couldn't open " << MuonList << endl;
    	return;
    }

	long utc;
	bool badScaler;
	int run, type, numBadScalers=0;
	double timeBefore=0, timeAfter=0, badScalerWindow=0, deadTime=0, deadBadScalerTime=0, hitTime=0;

	if (dsNumber != 4){
		timeBefore = .0002;		// 0.2 ms before
		timeAfter = 1;			// 1 sec after
		badScalerWindow = 8;  	// +/- 8 sec
	}
	else if (dsNumber == 4){
		// for now, use a larger window for M2 because we are less sure of clock sync
		timeBefore = 2;
		timeAfter = 2;
		badScalerWindow = 8;
	}

	while(true)
	{
		InputList >> run >> utc >> hitTime >> type >> badScaler;
		if (InputList.eof()) break;

		if (type==1 || type==2)
		{
			if (badScaler) {
				deadTime += 2*badScalerWindow;
				deadBadScalerTime += 2*badScalerWindow;
			}
			else deadTime += timeBefore + timeAfter;
		}
		if (type == 3) {
			if (badScaler) {
				deadTime += badScalerWindow;
				deadBadScalerTime += badScalerWindow;
			}
			else deadTime += timeAfter;
		}
		if (badScaler) numBadScalers++;
	}
	printf("Dead time due to veto: %.2f seconds.\n",deadTime);
	if (numBadScalers > 0) printf("Bad scalers: %i, %.2f of %.2f sec (%.2f%%)\n", numBadScalers,deadBadScalerTime,deadTime,((double)deadBadScalerTime/deadTime)*100);
}

// ==========================================

int PanelMap(int qdcChan, int runNum)
{
	// For Dave and Bradley.
	// Using the 32-panel configuration as "standard", map the QDC index to a physical panel location for all of the run ranges.
	// The figure "Veto Panels, View From The Top" in Veto System Change Log is the map I use for ALL configurations.
	// key: "qdc channel"  value: "panel location"

	map<int,int> panels;

	// 32-panel config (default) - began 7/10/15
	if (runNum < 45000000 && runNum >= 3057) {
		map<int,int> tempMap =
		{{0,1},  {1,2},  {2,3},  {3,4}, {4,5}, {5,6},
		 {6,7},  {7,8},  {8,9},  {9,10}, {10,11}, {11,12},
		 {12,13}, {13,14}, {14,15}, {15,16},
		 {16,17}, {17,18}, {18,19}, {19,20},
		 {20,21}, {21,22}, {22,23}, {23,24},
		 {24,25}, {25,26}, {26,27}, {27,28},
		 {28,29}, {29,30}, {30,31}, {31,32}};
		panels = tempMap;
	}
	// 1st prototype config (24 panels)
	else if (runNum >= 45000509 && runNum <= 45004116) {
		map<int,int> tempMap =
		{{0,1},  {1,2},  {2,3},  {3,4}, {4,5}, {5,6},
		 {6,7},  {7,8},  {8,9},  {9,10}, {10,11}, {11,12},
		 {12,21}, {13,22}, {14,15}, {15,16},
		 {16,17}, {17,18}, {18,19}, {19,20},
		 {20,13}, {21,14}, {22,23}, {23,24},
		 {24,-1}, {25,-1}, {26,-1}, {27,-1},
		 {28,-1}, {29,-1}, {30,-1}, {31,-1}};
		 panels = tempMap;
	}
	// 2nd prototype config (24 panels)
	else if (runNum >= 45004117 && runNum <= 45008659) {
		map<int,int> tempMap =
		{{0,1},  {1,2},  {2,3},  {3,4}, {4,5}, {5,6},
		 {6,7},  {7,8},  {8,9},  {9,10}, {10,11}, {11,12},
		 {12,13}, {13,14}, {14,15}, {15,16},
		 {16,17}, {17,18}, {18,19}, {19,20},
		 {20,21}, {21,22}, {22,23}, {23,24},
		 {24,-1}, {25,-1}, {26,-1}, {27,-1},
		 {28,-1}, {29,-1}, {30,-1}, {31,-1}};
		 panels = tempMap;
	}
	// 1st module 1 (P3JDY) config, 6/24/15 - 7/7/15
	else if (runNum > 0 && runNum <=3056){
		map<int,int> tempMap =
		{{0,1},  {1,2},  {2,3},  {3,4}, {4,5}, {5,6},
		 {6,7},  {7,8},  {8,9},  {9,10}, {10,11}, {11,12},
		 {12,13}, {13,14}, {14,15}, {15,16},
		 {16,17}, {17,18}, {18,19}, {19,20},
		 {20,21}, {21,22}, {22,23}, {23,24},
		 {24,-1}, {25,-1}, {26,-1}, {27,-1},
		 {28,-1}, {29,-1}, {30,-1}, {31,-1}};
		 panels = tempMap;
	}
	else {
		cout << "Panel map not known for this run number!\n";
		return -1;
	}
	int panel = -1;
	auto search = panels.find(qdcChan);
	if(search != panels.end()) panel=search->second;
	return panel;
}
