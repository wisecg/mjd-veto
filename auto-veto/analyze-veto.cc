// analyze-veto.cc
// Takes output of auto-veto.cc for veto analysis.
// C. Wiseman, 10/24/16

#include <iostream>
#include <fstream>
#include <numeric>
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TChain.h"
#include "TH1.h"
#include "TFile.h"
#include "MJVetoEvent.hh"

using namespace std;

void GenerateVetoList(TChain *vetoTree);
void GenerateDisplayList(TChain *vetoTree);
void CalculateDeadTime(string MuonList, int dsNumber);

int PanelMap(int i, int runNum);

int main(int argc, char** argv)
{
	if (argc < 1) {
		cout << "Usage: ./analyze-veto [run number] [optional: upper run number]\n";
		return 0;
	}
	int run = stoi(argv[1]);
	int hiRun = run;
	if (argc > 2) {
		hiRun = stoi(argv[2]);
		printf("Analyzing runs %i through %i ...\n",run,hiRun);
	}
	else printf("Analyzing run %i ... \n",run);

	TChain *vetoTree = new TChain("vetoTree");
	vector<int> runs(hiRun - run + 1);
	iota(runs.begin(), runs.end(), run);
	for (auto i : runs) {
		char file[200];
		sprintf(file,"./output/veto-run%i.root",i);
		vetoTree->Add(file);
	}
	printf("Found %lli entries.\n",vetoTree->GetEntries());
	// alternately, could just do something like
	// vetoTree->Add("./output/*.root")
	// or maybe need to modify the options to take in a run list

	GenerateVetoList(vetoTree);
	GenerateDisplayList(vetoTree);
	CalculateDeadTime("./output/MuonList_test.txt",1);
}

void GenerateVetoList(TChain *vetoTree)
{
	ofstream MuonList("./output/MuonList_test.txt");

	TTreeReader reader(vetoTree);
	TTreeReaderValue<MJVetoEvent> events(reader,"events");
	TTreeReaderValue<Long64_t> start(reader,"start");
	TTreeReaderValue<Long64_t> stop(reader,"stop");
	TTreeReaderValue<double> xTime(reader,"xTime");
	TTreeReaderArray<int> CoinType(reader,"CoinType");	//[32]

	bool newRun=false;
	int prevRun=0;
	Long64_t prevStop=0;

	while(reader.Next())
	{
		MJVetoEvent veto = *events;
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
				sprintf(muonList,"%i %lli %.8f %i %i\n",run,*start,*xTime,type,veto.GetBadScaler());
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
