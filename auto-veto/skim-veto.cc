// skim-veto.cc
// Takes output of auto-veto.cc for veto analysis.
// C. Wiseman, 10/24/16

#include <iostream>
#include <iomanip>
#include <fstream>
#include <numeric>
#include <cmath>
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
void GenerateDisplayList(TChain *vetoTree);
void CalculateDeadTime(string MuonList, int dsNumber);
void GetRunInfo();
void GenerateDS4MuonList();
void ListRunOffsets(TChain *vetoTree);
int PanelMap(int i, int runNum);

int main(int argc, char** argv)
{
	if (argc < 2) {
		cout << "Usage: ./skim-veto -r [run number]\n"
         << "       ./skim-veto [run list file]\n"
         << "       ./skim-veto -4 (ds4 list gen)";
		return 0;
	}
  TChain *vetoTree = new TChain("vetoTree");

  int run = 0;
  char file[200];
  string opt1 = argv[1];
  if (opt1 == "-r") {
    run = stoi(argv[2]);
    sprintf(file,"./avout/DS3/veto_run%i.root",run);
    if (!vetoTree->Add(file)){
      cout << "File doesn't exist.  Exiting ... \n";
      return 1;
    }
  }
  else if (opt1 == "-4"){
    // GetRunInfo();
    GenerateDS4MuonList();
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
      sprintf(file,"./avout/DS3/veto_run%i.root",i);
      if (!vetoTree->Add(file)){
        cout << "File doesn't exist.  Continuing ... \n";
        continue;
      }
    }
  }
	// printf("Found %lli entries.\n",vetoTree->GetEntries());
  // ListRunOffsets(vetoTree);
	// GenerateVetoList(vetoTree);
	// GenerateDisplayList(vetoTree);
	// CalculateDeadTime("./output/MuonList_test.txt",1);
}

void GenerateVetoList(TChain *vetoTree)
{
  // OBSOLETE: Need to use the auto-veto files.

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
  // OBSOLETE:  Need to use the auto-veto files.

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
  double lastRunTS = 0;
  while(reader.Next())
  {
    if(runSave != *runIn) // run boundary condition
    {
      runSave = *runIn;
      MJVetoEvent veto = *vetoEventIn;
      printf("Run %i  Entry %i  xTime %-5.3f  Offset %-5.3f  Unc %-5.3f\n", veto.GetRun(),veto.GetEntry(),*xTimeIn,*scalerOffsetIn,*scalerUncIn);
      hUnc->Fill(*scalerUncIn);
      if (*xTimeIn < lastRunTS) cout << "Clock reset, run " << veto.GetRun() << endl;
      lastRunTS = *xTimeIn;
    }
  }
  TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",800,600);
  hUnc->Draw();
  // hUnc->Fit("gaus");  // not really gaussian, probably flat with a big spike at 0
  c1->Print("./output/scalerUnc.pdf");
}

void GetRunInfo()
{
  // TODO: When GetStartTimeStamp becomes available,
  // need to regenerate these lists to use it.
  
  int run = 0;
  ifstream runFile("./runs/ds3-complete.txt");
  ofstream infoFile("./runs/ds3-runInfo.txt");
  // ifstream runFile("./runs/ds4-complete.txt");
  // ofstream infoFile("./runs/ds4-runInfo.txt");
  vector<int> runList;
  vector<double> digStarts;
  vector<double> unixStarts;
  vector<double> unixStops;
  while (runFile >> run) runList.push_back(run);
  for (auto i : runList)
  {
    GATDataSet *ds = new GATDataSet(i);
    TChain *gatChain = ds->GetGatifiedChain(false);
    TTreeReader gatReader(gatChain);
    TTreeReaderValue<double> startTimeIn(gatReader, "startTime");
    TTreeReaderValue<double> stopTimeIn(gatReader, "stopTime");
    TTreeReaderValue< vector<double> > timestampIn(gatReader, "timestamp");
    gatReader.SetEntry(0);
    double firstGretinaTS = (*timestampIn)[0]*1.e-8;
    digStarts.push_back(firstGretinaTS);
    unixStarts.push_back(*startTimeIn);
    unixStops.push_back(*stopTimeIn);
    printf("%i  %-8.3f  %li  %li\n",i,firstGretinaTS,(long)(*startTimeIn),(long)(*stopTimeIn));
    infoFile << i << "  " << setprecision(15) << firstGretinaTS << "  "
                << (long)(*startTimeIn) << "  " << (long)(*stopTimeIn) << endl;
    delete ds;
  }
  runFile.close();
  infoFile.close();
}

void GenerateDS4MuonList()
{
  int run;
  double digStart, unixStart, unixStop;

  // load run info
  ifstream ds3infoFile("./runs/ds3-runInfo.txt");
  vector<int> ds3runList;
  vector<double> ds3digStarts;
  vector<double> ds3unixStarts;
  vector<double> ds3unixStops;
  while (ds3infoFile >> run >> digStart >> unixStart >> unixStop){
    ds3runList.push_back(run);
    ds3digStarts.push_back(digStart);
    ds3unixStarts.push_back(unixStart);
    ds3unixStops.push_back(unixStop);
  }
  ds3infoFile.close();
  ifstream ds4infoFile("./runs/ds4-runInfo.txt");
  vector<int> ds4runList;
  vector<double> ds4digStarts;
  vector<double> ds4unixStarts;
  vector<double> ds4unixStops;
  while (ds4infoFile >> run >> digStart >> unixStart >> unixStop){
    ds4runList.push_back(run);
    ds4digStarts.push_back(digStart);
    ds4unixStarts.push_back(unixStart);
    ds4unixStops.push_back(unixStop);
  }
  ds4infoFile.close();

  // load veto data
  TChain *vetoTree = new TChain("vetoTree");
  for (auto i : ds3runList) {
    if (!vetoTree->Add(TString::Format("./avout/DS3/veto_run%i.root",i))){
      cout << "File doesn't exist.  Continuing ... \n";
      continue;
    }
  }
  vector<int> muRuns;
  vector<double> muRunTStarts_s;
  vector<double> muTimes;
  vector<int> muTypes;
  vector<double> muUncert;
  TTreeReader vetoReader(vetoTree);
  TTreeReaderValue<MJVetoEvent> vetoEventIn(vetoReader,"vetoEvent");
  TTreeReaderValue<int> vetoRunIn(vetoReader,"run");
	TTreeReaderValue<Long64_t> vetoStart(vetoReader,"start");
	TTreeReaderValue<Long64_t> vetoStop(vetoReader,"stop");
	TTreeReaderValue<double> xTime(vetoReader,"xTime");
  TTreeReaderValue<double> scalerUnc(vetoReader,"scalerUnc");
	TTreeReaderArray<int> CoinType(vetoReader,"CoinType");	//[32]
  bool newRun=false;
	int prevRun=0;
	Long64_t prevStop=0;
	while(vetoReader.Next())
	{
    MJVetoEvent veto = *vetoEventIn;
    int run = *vetoRunIn;
		if (run != prevRun) newRun=true;
		else newRun = false;
		int type = 0;
		if (CoinType[0]) type=1;
		if (CoinType[1]) type=2;	// overrides type 1 if both are true
		if ((*vetoStart-prevStop) > 10 && newRun) type = 3;
    if (type > 0){
      muRuns.push_back(run);
      muRunTStarts_s.push_back(*vetoStart);
      muTypes.push_back(type);
      if (type!=3) muTimes.push_back(*xTime);
      else muTimes.push_back(*xTime); // time of the first veto entry in the run
      if (!veto.GetBadScaler()) muUncert.push_back(*scalerUnc);
      else muUncert.push_back(8.0); // uncertainty for corrupted scalers
    }
		prevStop = *vetoStop;  // end of entry, save the run and stop time
		prevRun = run;
	}

  // Convert the ds-3 muon list into ds-4 muon list.
  vector<int> ds4muRuns;
  vector<double> ds4muRunTStarts;
  vector<double> ds4muTimes;
  vector<int> ds4muTypes;
  vector<double> ds4muUncert;
  for (int i = 0; i < (int)muRuns.size(); i++)
  {
    auto it = find(ds3runList.begin(), ds3runList.end(), muRuns[i]);
    if (it != ds3runList.end())
    {
      auto idx = distance(ds3runList.begin(), it);
      double t_global_mu1 = (muTimes[i] - ds3digStarts[idx]) + ds3unixStarts[idx];

      for (int idx2 = 0; idx2 < (int)ds4runList.size(); idx2++)
      {
        if (t_global_mu1 > ds4unixStarts[idx2] && t_global_mu1 < ds4unixStops[idx2])
        {
          double t_mu2 = (t_global_mu1 - ds4unixStarts[idx2]) + ds4digStarts[idx2];  // relative time + digitizer start time
          double t_mu2_uncert = sqrt(pow(1.e-8,2) + pow(1.e-8,2) + 2*pow(1.,2) + pow(muUncert[i],2));
          ds4muRuns.push_back(ds4runList[idx2]);
          ds4muRunTStarts.push_back(ds4unixStarts[idx2]);
          ds4muTimes.push_back(t_mu2);
          ds4muUncert.push_back(t_mu2_uncert);
          ds4muTypes.push_back(muTypes[i]);
          printf("%lu  glob %li  %i (%-8.2fs)  %i (%-8.2fs)  ds4loc %-6.2f +/- %-4.2f\n", ds4muRuns.size(), (long)t_global_mu1, ds4runList[idx2], t_global_mu1-ds4unixStarts[idx2], ds3runList[idx], muTimes[i]-ds3digStarts[idx], t_mu2, t_mu2_uncert);
          break;
        }
      }
    }
    else cout << "Couldn't find this run in the muon list.\n";
  }
  // check muon list
  cout << ds4muRuns.size() << " of " << muRuns.size() << " DS-3 muon candidates persisted in DS-4.\n";
  // for (int i = 0; i < (int)ds4muRuns.size(); i++)
    // printf("%i  %i  %.0f  %.3f  %.3f\n",ds4muRuns[i],ds4muTypes[i],ds4muRunTStarts[i],ds4muTimes[i],ds4muUncert[i]);


  // Create the dictionaries for use in $GATDIR/Apps/DataSetInfo.hh

  // cout << "vector<int> ds4muRuns = {";
  // for (int i = 0; i < (int)ds4muRuns.size()-1; i++) cout << ds4muRuns[i] << ", ";
  // cout << ds4muRuns[(int)ds4muRuns.size()-1] << "};\n\n";
  //
  // cout << "vector<int> ds4muTypes = {";
  // for (int i = 0; i < (int)ds4muTypes.size()-1; i++) cout << ds4muTypes[i] << ", ";
  // cout << ds4muTypes[(int)ds4muTypes.size()-1] << "};\n\n";
  //
  // cout << "vector<int> ds4muRunTStarts = {";
  // for (int i = 0; i < (int)ds4muRunTStarts.size()-1; i++) cout << (long)ds4muRunTStarts[i] << ", ";
  // cout << (long)ds4muRunTStarts[(int)ds4muRunTStarts.size()-1] << "};\n\n";
  //
  // cout << "vector<int> ds4muTimes = {";
  // for (int i = 0; i < (int)ds4muTimes.size()-1; i++) cout << setprecision(9) << ds4muTimes[i] << ", ";
  // cout << setprecision(9) << ds4muTimes[(int)ds4muTimes.size()-1] << "};\n\n";
  //
  // cout << "vector<int> ds4muUncert = {";
  // for (int i = 0; i < (int)ds4muUncert.size()-1; i++) cout << setprecision(9) << ds4muUncert[i] << ", ";
  // cout << setprecision(9) << ds4muUncert[(int)ds4muUncert.size()-1] << "};\n\n";
}
