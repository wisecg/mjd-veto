// A very, very stripped down version of
// skim_mjd_data, for the purpose of checking
// muon-ge coincidences.
// C. Wiseman, 2016/11/19

#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TChain.h"
#include "TCut.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include "TEntryList.h"
#include "TROOT.h"

#include "GATDataSet.hh"
#include "MJTChannelMap.hh"
#include "TClonesArray.h"
#include "MGTEvent.hh"
#include "GATUtils.hh"
#include "MJVetoEvent.hh"

#include "DataSetInfo.hh"

using namespace std;
using namespace CLHEP;

void LoadDataSet(GATDataSet& ds, int dsNumber, size_t iRunSeq);
void LoadRun(GATDataSet& ds, size_t iRunSeq);

int main(int argc, const char** argv)
{
  if(argc < 3 || argc > 7) {
    cout << "Usage for single file: " << argv[0] << " -f [runNum] (output path)" << endl;
    cout << "Usage for data sets: " << argv[0] << " [dataset number] [runseq] (output path)" << endl;
    cout << "Usage for run lists: " << argv[0] << " -l [dataset number] [path to txt list] (output path)" << endl;
    return 1;
  }

  GATDataSet ds;
  TChain *vetoChain = new TChain("vetoTree");
  string outputPath = "";
  int noDS = -999999;
  int dsNumber = noDS;
  int runSeq=0;
  bool singleFile = false;
  bool runList = false;
  TString flags ="";
  int i = 1;
  while(!TString(argv[i]).IsDigit()){
    flags += TString(argv[i]);
    i++;
  }
  if(flags.Contains("f")){
    singleFile = true;
    runSeq = atoi(argv[i]);
    i++;
    if(runSeq >= 2580 && runSeq <= 6963) dsNumber = 0;
    else if(runSeq >= 9422 && runSeq <= 14502) dsNumber = 1;
    else if(runSeq >= 14503 && runSeq <= 15892) dsNumber = 2;
    else if(runSeq >= 16797 && runSeq <= 17980) dsNumber = 3;
    else if(runSeq >= 60000802 && runSeq <= 60001888) dsNumber = 4;
    else {
      cout << "Error: I don't know what dataset run " << runSeq << " is from." << endl;
      return 1;
    }
    cout << "Loading run " << runSeq << endl;
    ds.AddRunNumber(runSeq);
    if (!vetoChain->Add(TString::Format("./avout/veto_run%i.root",runSeq))){
      cout << "Veto files not found.  Exiting ...\n";
      return 1;
    }
  }
  else if (flags.Contains("l")){
    runList = true;
    dsNumber = atoi(argv[i]);
    i++;
    ifstream runFile(argv[i]);
    int rundummy;
    // make sure we only load unique, sequential runs
    set<int> uniqueRuns;
    while (runFile >> rundummy) uniqueRuns.insert(rundummy);
    vector<int> runList(uniqueRuns.begin(), uniqueRuns.end());
    sort(runList.begin(), runList.end());
    for (auto i : runList) {
      cout << "Loading run " << i << endl;
      ds.AddRunNumber(i);
      if (!vetoChain->Add(TString::Format("./avout/veto_run%i.root",i))){
        cout << "Veto files not found.  Exiting ...\n";
        return 1;
      }
    }
    i++;
  }
  else {
    // TODO: add a way to get the vetoChain for an entire dataset
    dsNumber = atoi(argv[i]);
    i++;
    runSeq = atoi(argv[i]);
    i++;
    cout << "loading dataset " << dsNumber << " run sequence " << runSeq << endl;
    LoadDataSet(ds, dsNumber, runSeq);
  }
  cout << "Found " << vetoChain->GetEntries() << " veto entries.\n";
  if(argc > i) outputPath += string(argv[i]);

  // set up dataset
  TChain* gatChain = ds.GetGatifiedChain(false);
  TTreeReader gatReader(gatChain);

  // set up input chain value readers

  // run level variables and indices
  TTreeReaderValue<unsigned int> gatrevIn(gatReader, "gatrev");
  TTreeReaderValue<double> runIn(gatReader, "run");
  double runSave = -1;

  // ID variables
  TTreeReaderValue< vector<double> > channelIn(gatReader, "channel");
  TTreeReaderValue< vector<int> > detIDIn(gatReader, "detID");
  TTreeReaderValue< vector<string> > detNameIn(gatReader, "detName");
  TTreeReaderValue< vector<int> > posIn(gatReader, "P");
  TTreeReaderValue< vector<int> > detIn(gatReader, "D");
  TTreeReaderValue< vector<int> > cryoIn(gatReader, "C");
  TTreeReaderValue< vector<bool> > isEnrIn(gatReader, "isEnr");
  TTreeReaderValue< vector<bool> > isNatIn(gatReader, "isNat");

  // time variables
  TTreeReaderValue<double> startTimeIn(gatReader, "startTime");
  TTreeReaderValue<double> stopTimeIn(gatReader, "stopTime");
  TTreeReaderValue< vector<double> > timestampIn(gatReader, "timestamp");
  TTreeReaderValue< vector<double> > timeMTIn(gatReader, "timeMT");
  TTreeReaderValue< vector<int> > dateMTIn(gatReader, "dateMT");

  // energy variables
  TTreeReaderValue< vector<double> > trapENFDBSGCalIn(gatReader, "trapENFDBSGCal");
  TTreeReaderValue< vector<double> > trapECalIn(gatReader, "trapECal");

  // data cleaning variables
  const int kPinghanPulserMask = 0x1 << 1; // pinghan pulsers
  TTreeReaderValue<unsigned int> eventDC1BitsIn(gatReader, "EventDC1Bits");

  // temporary: get first timestamp in chain
  gatReader.SetEntry(0);
  double firstGretinaTS = (*timestampIn)[0];
  gatReader.SetTree(gatChain); // reset the reader

  // load muon data
  cout << "Loading muon data..." << endl;
  vector<int> muRuns;
  vector<double> muRunTStarts_s;
  vector<double> muTimes_s;
  vector<int> muTypes;
  vector<bool> badScalers;
  TTreeReader vetoReader(vetoChain);
	TTreeReaderValue<MJVetoEvent> vetoEventIn(vetoReader,"vetoEvent");
	TTreeReaderValue<Long64_t> vetoStart(vetoReader,"start");
	TTreeReaderValue<Long64_t> vetoStop(vetoReader,"stop");
	TTreeReaderValue<double> xTime(vetoReader,"xTime");
	TTreeReaderArray<int> CoinType(vetoReader,"CoinType");	//[32]
  bool newRun=false;
	int prevRun=0;
	Long64_t prevStop=0;
	while(vetoReader.Next())
	{
		MJVetoEvent veto = *vetoEventIn;
		int run = veto.GetRun();
		if (run != prevRun) newRun=true;
		else newRun = false;
		int type = 0;
		if (CoinType[0]) type=1;
		if (CoinType[1]) type=2;	// overrides type 1 if both are true
		if ((*vetoStart-prevStop) > 10 && newRun) type = 3;
    if (type > 0){
      muRuns.push_back(run);
      muRunTStarts_s.push_back(*vetoStart);
      if (type!=3) muTimes_s.push_back(*xTime);
      else muTimes_s.push_back(0.0);
      muTypes.push_back(type);
      if (type!=3) badScalers.push_back(veto.GetBadScaler());
      else badScalers.push_back(0);
    }
		prevStop = *vetoStop;  // end of entry, save the run and stop time
		prevRun = run;
	}
  delete vetoChain;
  size_t iMu = 0;  // set up an iterator for the veto list
  size_t nMu = muTimes_s.size();
  if(nMu == 0) {
    cout << "couldn't load mu data" << endl;
    return 0;
  }
  // check muon list output
  for (int i = 0; i < (int)muRuns.size(); i++)
    cout << muRuns[i] << "  " << muRunTStarts_s[i]  << "  " << muTimes_s[i] << "  " << muTypes[i] << "  " << badScalers[i] << endl;

  // set up output file and tree
  cout << "creating output file ...\n";
  string filename = TString::Format("skimDS%d_", dsNumber).Data();
  char runStr[10];
  sprintf(runStr, "%d", runSeq);
  if(singleFile) filename += "run";
  if(runList) sprintf(runStr,"list");
  filename += runStr;
  filename += ".root";
  if(outputPath != "") filename = outputPath + "/" + filename;
  TFile *fOut = TFile::Open(filename.c_str(), "recreate");
  TTree* skimTree = new TTree("skimTree", "skimTree");

  // set up output chain branches

  // run level variables and indices
  int run = 0;
  skimTree->Branch("run", &run, "run/I");
  int iEvent = 0;
  skimTree->Branch("iEvent", &iEvent, "iEvent/I");
  vector<int> iHit;
  skimTree->Branch("iHit", &iHit);

  // ID variables
  vector<int> channel;
  skimTree->Branch("channel", &channel);
  vector<int> pos;
  skimTree->Branch("P", &pos);
  vector<int> det;
  skimTree->Branch("D", &det);
  vector<string> detName;
  skimTree->Branch("detName", &detName);
  vector<int> cryo;
  skimTree->Branch("C", &cryo);
  vector<int> gain;
  skimTree->Branch("gain", &gain);
  vector<int> detID;
  skimTree->Branch("detID", &detID);
  vector<bool> isEnr;
  skimTree->Branch("isEnr", &isEnr);
  vector<bool> isNat;
  skimTree->Branch("isNat", &isNat);
  vector<bool> isGood;
  skimTree->Branch("isGood", &isGood);

  // time variables
  double runTime_s = 0;
  double startTime = 0;
  double startTime0 = 0;
  double stopTime = 0;
  if(dsNumber == 0) {
    runTime_s = 4121110;
    startTime0 = 1435687000; // start time of run 2580
  }
  else if(dsNumber == 1) {
    runTime_s = 5282280; //was 4728790; before blinding
    startTime0 = 1452643100; // start time of run 9422
  }
  else if(dsNumber == 3) {
    runTime_s = 2584470;//1549710;
    startTime0 = 1472169600; // start time of run 16797
  }
  else if(dsNumber == 4) {
    runTime_s = 2060020;//1460390;
    startTime0 = 1472169600; // start time of run 60000802
  }
  skimTree->Branch("startTime", &startTime, "startTime/D");
  skimTree->Branch("startTime0", &startTime0, "startTime0/D");
  skimTree->Branch("runTime_s", &runTime_s, "runTime_s/D");
  skimTree->Branch("stopTime", &stopTime, "stopTime/D");
  vector<double> tloc_s;
  skimTree->Branch("tloc_s", &tloc_s);
  vector<double> time_s;
  skimTree->Branch("time_s", &time_s);
  vector<double> timeMT;
  skimTree->Branch("timeMT", &timeMT);
  vector<int> dateMT;
  skimTree->Branch("dateMT", &dateMT);

  // energy variables
  vector<double> trapENFDBSGCal;
  skimTree->Branch("trapENFDBSGCal", &trapENFDBSGCal);
  vector<double> trapECal;
  skimTree->Branch("trapECal", &trapECal);

  // analysis cut variables
  int mH = 0;
  skimTree->Branch("mH", &mH, "mH/I");
  int mL = 0;
  skimTree->Branch("mL", &mL, "mL/I");

  // data cleaning variables
  unsigned int eventDC1Bits = 0;
  skimTree->Branch("EventDC1Bits", &eventDC1Bits, "eventDC1Bits/i");

  // veto variables
  vector<double> dtmu_s;
  skimTree->Branch("dtmu_s", &dtmu_s);
  vector<int> muType;
  vector<bool> badScaler;
  skimTree->Branch("muType", &muType);
  skimTree->Branch("badScaler", &badScaler);
  vector<bool> muVeto;
  skimTree->Branch("muVeto", &muVeto);
  map<int,bool> detIDIsBad;
  if(dsNumber == 0) {
    detIDIsBad[28474] = true;
    detIDIsBad[1426622] = true;
    detIDIsBad[28480] = true;
    detIDIsBad[1426980] = true;
    detIDIsBad[1426620] = true;
    detIDIsBad[1425370] = true;
  }
  if(dsNumber == 1) {
    detIDIsBad[1426981] = true;
    detIDIsBad[1426622] = true;
    detIDIsBad[28455] = true;
    detIDIsBad[28470] = true;
    detIDIsBad[28463] = true;
    detIDIsBad[28465] = true;
    detIDIsBad[28469] = true;
    detIDIsBad[28477] = true;
    detIDIsBad[1425751] = true;
    detIDIsBad[1425731] = true;
    detIDIsBad[1426611] = true;
  }
  if(dsNumber == 3) {
    detIDIsBad[1426981] = true;
    detIDIsBad[1426622] = true;
    detIDIsBad[28477] = true;
    detIDIsBad[1425731] = true;
    detIDIsBad[1426611] = true;
  }
  if(dsNumber == 4) {
    detIDIsBad[28595] = true;
    detIDIsBad[28461] = true;
    detIDIsBad[1428530] = true;
    detIDIsBad[28621] = true;
    detIDIsBad[28473] = true;
    detIDIsBad[1426651] = true;
    detIDIsBad[1429092] = true;
    detIDIsBad[1426652] = true;
    detIDIsBad[28619] = true;
  }
  map<int,bool> detIDIsVetoOnly;
  if(dsNumber == 0) {
    detIDIsVetoOnly[1425381] = true;
    detIDIsVetoOnly[1425742] = true;
  }
  if(dsNumber == 1) {
    detIDIsVetoOnly[28480] = true;
    detIDIsVetoOnly[1426621] = true;
  }
  if(dsNumber == 3) {
    detIDIsVetoOnly[28480] = true;
    detIDIsVetoOnly[28470] = true;
    detIDIsVetoOnly[28463] = true;
  }
  if(dsNumber == 4) {
    detIDIsVetoOnly[28459] = true;
    detIDIsVetoOnly[1426641] = true;
    detIDIsVetoOnly[1427481] = true;
    detIDIsVetoOnly[28456] = true;
    detIDIsVetoOnly[1427120] = true;
    detIDIsVetoOnly[1427121] = true;
  }

  // start loop over all events
  while(gatReader.Next()) {

    // stuff to do on run boundaries
    if(runSave != *runIn) {
      runSave = *runIn;
      firstGretinaTS = (*timestampIn)[0];
      // cout << "Processing run " << *runIn << ", "
      //      << skimTree->GetEntries() << " entries saved so far"
      //      << endl;
      skimTree->Write("", TObject::kOverwrite);
    }

    // Skip this event if it is a pulser event as identified by Pinghan
    if(*eventDC1BitsIn & kPinghanPulserMask) continue;
    iEvent = gatChain->GetTree()->GetReadEntry();

    // copy the event-level info to the output fields
    run = int(*runIn);
    startTime = *startTimeIn;
    stopTime = *stopTimeIn;
    mH = 0;
    mL = 0;

    // clear all hit-level info fields
    iHit.resize(0);
    muType.resize(0);
    badScaler.resize(0);
    dtmu_s.resize(0);
    trapECal.resize(0);
    trapENFDBSGCal.resize(0);
    channel.resize(0);
    tloc_s.resize(0);
    time_s.resize(0);
    timeMT.resize(0);
    dateMT.resize(0);
    pos.resize(0);
    det.resize(0);
    cryo.resize(0);
    gain.resize(0);
    detID.resize(0);
    isEnr.resize(0);
    isNat.resize(0);
    isGood.resize(0);
    muVeto.resize(0);

    // loop over hits
    bool skipMe = false;
    size_t nHits = trapENFDBSGCalIn->size();
    for(size_t i=0; i<nHits; i++)
    {
      // original note:
      // skip all hits with E_H < 2 keV, E_L < 10 keV in -both- trapE and trapENF
      // For small skim files, skip all hits with E_H and E_L < 200 keV in trapE and trapENF
      // clint's note:
      // skip all hits with E_H and E_L < 200 keV in trapENFDBSGCal and trapECal
      double hitENFDBSGCal = (*trapENFDBSGCalIn)[i];
      double hitECal = (*trapECalIn)[i];
      int hitCh = (*channelIn)[i];
      if(hitCh%2 == 0 && hitENFDBSGCal <  5000. && hitECal < 5000.) continue;
      if(hitCh%2 == 1 && hitENFDBSGCal < 5000. && hitECal < 5000.) continue;

      // skip hits from totally "bad" detectors (not biased, etc), or from
      // use-for-veto-only detectors if E < 10 keV
      int hitDetID = (*detIDIn)[i];
      if(detIDIsBad[hitDetID] || (detIDIsVetoOnly[hitDetID] && hitECal < 10.)) continue;
      // copy over hit info
      iHit.push_back(i);
      trapECal.push_back(hitECal);
      trapENFDBSGCal.push_back(hitENFDBSGCal);
      channel.push_back(hitCh);
      cryo.push_back((*cryoIn)[i]);
      pos.push_back((*posIn)[i]);
      det.push_back((*detIn)[i]);
      gain.push_back(hitCh % 2);
      detID.push_back(hitDetID);
      isEnr.push_back((*detNameIn)[i][0] == 'P');
      isNat.push_back((*detNameIn)[i][0] == 'B');
      isGood.push_back(!detIDIsVetoOnly[hitDetID]);
      // sum energies and multiplicities
      if(hitCh%2 == 0) mH++;
      else mL++;

      // copy over time info
      double hitT_s = (*timestampIn)[i]*1.e-8;
      tloc_s.push_back(hitT_s);
      time_s.push_back( (startTime - startTime0) + hitT_s ); //Need to figure out what to do with continuous running, Clara 10/10/16
      timeMT.push_back((*timeMTIn)[i]);
      dateMT.push_back((*dateMTIn)[i]);

      double hitT_s_adj = hitT_s - firstGretinaTS*1.e-8;

      // Calculate mu vetos
      double dtGoodMu_s = 0.0002;
      double dtBadMu_s = 8;
      while(1) {
        if(iMu >= nMu-1) break;
        double dt = muRunTStarts_s[iMu+1]-startTime;
        dt += muTimes_s[iMu+1] - hitT_s_adj;
        if(badScalers[iMu+1]) dt -= dtBadMu_s;
        else dt -= dtGoodMu_s;
        if(dt > 0) break;
        iMu++;
      }
      double dtmu = (startTime -muRunTStarts_s[iMu]) + hitT_s_adj - muTimes_s[iMu];
      double tmu = muTimes_s[iMu];
      dtmu_s.push_back(dtmu);
      muType.push_back(muTypes[iMu]);
      badScaler.push_back(badScalers[iMu]);

      if(badScalers[iMu] && fabs(dtmu < dtBadMu_s)) muVeto.push_back(true);

      else if (dsNumber==4) muVeto.push_back(dtmu > -2. && dtmu < 2.);  // DS4 runs

      else muVeto.push_back(dtmu > -0.2e-3 && dtmu < 1);

      if (hitCh%2==0) continue;
      printf("%i  raw %-8.0f  adj %-8.1f  tmu %-5.0f  dtmu %-8.1f  mL %-3i  e %-6.0f  iMu %-3lu\n", run,hitT_s,hitT_s_adj,tmu,dtmu,mL,hitENFDBSGCal,iMu);


    }

    // If no good hits in the event or skipped for some other reason, don't
    // write this event to the output tree.
    if(trapENFDBSGCal.size() == 0 || skipMe) continue;

    // finally, fill the tree for this event
    skimTree->Fill();

  }


  // write output tree to output file
  cout << "Closing out skim file..." << endl;
  skimTree->Write("", TObject::kOverwrite);

  fOut->Close();
  return 0;
}