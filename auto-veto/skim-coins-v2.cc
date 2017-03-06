// another massively paired down skim file generator
// for checking veto coins.
// C Wiseman, 2/2/2017

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
#include <vector>
#include <algorithm>
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
void LoadActiveMasses(map<int,double>& activeMassForDetID_g, int dsNumber);

int main(int argc, const char** argv)
{
  if(argc < 3 || argc > 9) {
    cout << "Usage for data sets: " << argv[0] << " [dataset number] [runseq] (output path)" << endl;
    return 1;
  }

  GATDataSet ds;
  TChain* gatChain = NULL;
  TChain* vetoChain = NULL;
  string outputPath = "";
  int noDS = -999999;
  int dsNumber = noDS;
  int runSeq;
  double energyThresh = 2.0;
  bool smallOutput = true;
  bool simulatedInput = false;
  vector<string> args;
  for(int iArg=0; iArg<argc; ++iArg) args.push_back(argv[iArg]);

  // This version only runs over datasets
  dsNumber = stoi(args[1]);
  runSeq = stoi(args[2]);
  args.erase(args.begin()+1, args.begin()+3);
  cout << "loading dataset " << dsNumber << " run sequence " << runSeq << endl;
  LoadDataSet(ds, dsNumber, runSeq);

  // set up dataset
  cout << " getting chain\n";
  if(gatChain==NULL) gatChain = ds.GetGatifiedChain(false);
  cout << "got chain\n";
  TTreeReader gatReader(gatChain);

  // set up input chain value readers

  // run level variables and indices
  TTreeReaderValue<unsigned int> gatrevIn(gatReader, "gatrev");
  TTreeReaderValue<double> runIn(gatReader, "run");
  double runSave = -1;

  // ID variables
  TTreeReaderValue< vector<double> > channelIn(gatReader, "channel");
  TTreeReaderValue< vector<int> > detIDIn(gatReader, "detID");
  TTreeReaderValue< vector<int> > posIn(gatReader, "P");
  TTreeReaderValue< vector<int> > detIn(gatReader, "D");
  TTreeReaderValue< vector<int> > cryoIn(gatReader, "C");
  TTreeReaderValue< vector<int> > mageIDIn(gatReader, "mageID");
  TTreeReaderValue< vector<string> > detNameIn(gatReader, "detName");
  TTreeReaderValue< vector<bool> > isEnrIn(gatReader, "isEnr");
  TTreeReaderValue< vector<bool> > isNatIn(gatReader, "isNat");

  // time variables
  TTreeReaderValue<double> startTimeIn(gatReader, "startTime");
  TTreeReaderValue<double> stopTimeIn(gatReader, "stopTime");
  TTreeReaderValue< vector<double> > timestampIn(gatReader, "timestamp");
    // We have to use pointers here because these variables will not be in
    // simulated data, and the only way to set the TTreeReader and branch name
    // is in the constructor of the object. Thanks ROOT!
  TTreeReaderValue< vector<double> >* timeMTIn;
  TTreeReaderValue< vector<int> >* dateMTIn;
  if(gatChain->GetListOfBranches()->FindObject("timeMT"))
  {
    timeMTIn = 0;
    dateMTIn = 0;
    timeMTIn = new TTreeReaderValue< vector<double> >(gatReader, "timeMT");
    dateMTIn = new TTreeReaderValue< vector<int> >(gatReader, "dateMT");
  }
  else
  {
    simulatedInput = true;
    timeMTIn = 0;
    dateMTIn = 0;
  }

  cout << " Setting up energy variables\n";
  // energy variables
  TTreeReaderValue< vector<double> > trapENFIn(gatReader, "trapENF");
  TTreeReaderValue< vector<double> > trapENFCalIn(gatReader, "trapENFCal");
  TTreeReaderValue< vector<double> > trapENMCalIn(gatReader, "trapENMCal");
  TTreeReaderValue< vector<double> > trapECalIn(gatReader, "trapECal");
  TTreeReaderValue< vector<double> > energyIn(gatReader, "energy");

  // analysis cut variables
//  TTreeReaderValue< vector<double> > blrwfFMR1In(gatReader, "blrwfFMR1"); //no longer in gatified tree
  TTreeReaderValue< vector<double> > blrwfFMR50In(gatReader, "blrwfFMR50");
  TTreeReaderValue< vector<double> > tsCurrent50nsMaxIn(gatReader, "TSCurrent50nsMax");
  TTreeReaderValue< vector<double> > tsCurrent100nsMaxIn(gatReader, "TSCurrent100nsMax");
  TTreeReaderValue< vector<double> > tsCurrent200nsMaxIn(gatReader, "TSCurrent200nsMax");
  TTreeReaderValue< vector<double> > triTrapMaxIn(gatReader, "triTrapMax");
//  TTreeReaderValue< vector<double> > toeIn(gatReader, "toe"); //no longer in gatified tree
  TTreeReaderValue< vector<double> > dcrSlopeIn(gatReader, "nlcblrwfSlope");

  // data cleaning variables
  TTreeReaderValue<unsigned int> eventDC1BitsIn(gatReader, "EventDC1Bits");
  // pulser tag channel seems to miss a lot of pulsers! So let's just use
  // Pinghan's (only).
  //const int kDoublePulserMask = (0x1 << 1) + 0x1; // pinghan pulsers + pulser tag channels
  const int kPinghanPulserMask = 0x1 << 1; // pinghan pulsers
  TTreeReaderValue< vector<unsigned int> > wfDCBitsIn(gatReader, "wfDCBits");
  TTreeReaderValue< vector<double> > trapETailMinIn(gatReader, "trapETailMin");
  TTreeReaderValue< vector<double> > nRisingXIn(gatReader, "fastTrapNLCWFsnRisingX");
  //TTreeReaderValue< vector<double> > d2wfnoiseTagNormIn(gatReader, "d2wfnoiseTagNorm");


  // Load muon data
  if(vetoChain==NULL) vetoChain = ds.GetVetoChain();
  cout << "Found " << vetoChain->GetEntries() << " veto entries.  Creating muon list ...\n";
  vector<int> muRuns;
  vector<int> muTypes;
  vector<double> muRunTStarts;
  vector<double> muTimes;
  vector<double> muUncert;
  if (dsNumber != 4 && !simulatedInput)
  {
    cout << " Creating veto list .. " << endl;
    TTreeReader vetoReader(vetoChain);
    TTreeReaderValue<MJVetoEvent> vetoEventIn(vetoReader,"vetoEvent");
    TTreeReaderValue<int> vetoRunIn(vetoReader,"run");
  	TTreeReaderValue<Long64_t> vetoStart(vetoReader,"start");
  	TTreeReaderValue<Long64_t> vetoStop(vetoReader,"stop");
  	TTreeReaderValue<double> xTime(vetoReader,"xTime");
    TTreeReaderValue<double> timeUncert(vetoReader,"timeUncert");
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
        muRunTStarts.push_back(*vetoStart);
        muTypes.push_back(type);
        if (type!=3) muTimes.push_back(*xTime);
        else muTimes.push_back(*xTime); // time of the first veto entry in the run
        if (!veto.GetBadScaler()) muUncert.push_back(*timeUncert);
        else muUncert.push_back(8.0); // uncertainty for corrupted scalers
      }
  		prevStop = *vetoStop;  // end of entry, save the run and stop time
  		prevRun = run;
  	}
  }
  else if (dsNumber==4 && !simulatedInput) LoadDS4MuonList(muRuns,muRunTStarts,muTimes,muTypes,muUncert);
  size_t iMu = 0;
  size_t nMu = muTimes.size();
  if(nMu == 0 && !simulatedInput) {
    cout << "couldn't load mu data" << endl;
    return 0;
  }
  cout << "Muon list has " << muRuns.size() << " entries.\n";
  // for (int i = 0; i < (int)muRuns.size(); i++)
    // printf("%i  %i  %i  %.0f  %.3f +/- %.3f\n",i,muRuns[i],muTypes[i],muRunTStarts[i],muTimes[i],muUncert[i]);


  // set up output file and tree
  string filename = TString::Format("skimDS%d_", dsNumber).Data();
  char runStr[10];
  sprintf(runStr, "%d", runSeq);
  filename += runStr;
  filename += ".root";
  if(outputPath != "") filename = outputPath + "/" + filename;
  TFile *fOut = TFile::Open(filename.c_str(), "recreate");
  TTree* skimTree = new TTree("skimTree", "skimTree");

  // set up output chain branches

  // run level variables and indices
  unsigned int skimgatrev = strtol(GATUtils::GetGATRevision(), NULL, 16);
  skimTree->Branch("skimgatrev", &skimgatrev, "skimgatrev/i");
  cout << "skimgatrev_" << skimgatrev << endl;
  unsigned int gatrev = 0;
  skimTree->Branch("gatrev", &gatrev, "gatrev/i");
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
  vector<int> cryo;
  skimTree->Branch("C", &cryo);
  vector<int> gain;
  skimTree->Branch("gain", &gain);
  vector<int> mageID;
  skimTree->Branch("mageID", &mageID);
  vector<int> detID;
  skimTree->Branch("detID", &detID);
  vector<string> detName;
  skimTree->Branch("detName", &detName);
  vector<bool> isEnr;
  skimTree->Branch("isEnr", &isEnr);
  vector<bool> isNat;
  skimTree->Branch("isNat", &isNat);
  map<int, double> actM4Det_g;
  LoadActiveMasses(actM4Det_g, dsNumber);
  vector<double> mAct_g;
  skimTree->Branch("mAct_g", &mAct_g);
  vector<bool> isGood;
  skimTree->Branch("isGood", &isGood);

  // total mass variables
  double mAct_M1Total_kg = 0;
  double mAct_M1enr_kg = 0;
  double mAct_M1nat_kg = 0;
  double mAct_M2Total_kg = 0;
  double mAct_M2enr_kg = 0;
  double mAct_M2nat_kg = 0;
  if(dsNumber == 0) {
    mAct_M1Total_kg = 14.60;
    mAct_M1enr_kg = 10.69;
    mAct_M1nat_kg = 3.91;
  }
  else if(dsNumber == 1) {
    mAct_M1Total_kg = 12.43;
    mAct_M1enr_kg = 11.31;
    mAct_M1nat_kg = 1.12;
  }
  else if(dsNumber == 3) {
    mAct_M1Total_kg = 15.412;
    mAct_M1enr_kg = 12.631;
    mAct_M1nat_kg = 2.781;
  }
  else if(dsNumber == 4) {
    mAct_M2Total_kg = 9.4212;
    mAct_M2enr_kg = 5.4712;
    mAct_M2nat_kg = 3.9500;
  }
  else if(dsNumber == 5) {
    mAct_M1Total_kg = 15.412;
    mAct_M1enr_kg = 12.631;
    mAct_M1nat_kg = 2.781;
    mAct_M2Total_kg = 9.4212;
    mAct_M2enr_kg = 5.4712;
    mAct_M2nat_kg = 3.9500;
  }
  skimTree->Branch("mAct_M1Total_kg", &mAct_M1Total_kg, "mAct_M1Total_kg/D");
  skimTree->Branch("mAct_M1enr_kg", &mAct_M1enr_kg, "mAct_M1enr_kg/D");
  skimTree->Branch("mAct_M1nat_kg", &mAct_M1nat_kg, "mAct_M1nat_kg/D");
  skimTree->Branch("mAct_M2Total_kg", &mAct_M2Total_kg, "mAct_M2Total_kg/D");
  skimTree->Branch("mAct_M2enr_kg", &mAct_M2enr_kg, "mAct_M2enr_kg/D");
  skimTree->Branch("mAct_M2nat_kg", &mAct_M2nat_kg, "mAct_M2nat_kg/D");

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
    runTime_s = 5197670; //was 4728790; before blinding
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
  else if (dsNumber == 5) {
    runTime_s = 7129520;
    startTime0 = 1476396800; // start time of run 18623
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
  vector<int> dateMT;
  if(!simulatedInput)
  {
    skimTree->Branch("timeMT", &timeMT);
    skimTree->Branch("dateMT", &dateMT);
  }

  // energy variables
  vector<double> trapECal;
  vector<double> onBoardE;
  if(!smallOutput){
    skimTree->Branch("trapECal", &trapECal);
    skimTree->Branch("onBoardE", &onBoardE);
  }
  vector<double> trapENFCal;
  skimTree->Branch("trapENFCal", &trapENFCal);
  vector<double> trapENMCal;
  skimTree->Branch("trapENMCal", &trapENMCal);
  double sumEH = 0;
  skimTree->Branch("sumEH", &sumEH, "sumEH/D");
  double sumEL = 0;
  skimTree->Branch("sumEL", &sumEL, "sumEL/D");
  double sumEHClean = 0;
  skimTree->Branch("sumEHClean", &sumEHClean, "sumEHClean/D");
  double sumELClean = 0;
  skimTree->Branch("sumELClean", &sumELClean, "sumELClean/D");

  // analysis cut variables
  int mH = 0;
  skimTree->Branch("mH", &mH, "mH/I");
  int mL = 0;
  skimTree->Branch("mL", &mL, "mL/I");
  int mHClean = 0;
  skimTree->Branch("mHClean", &mHClean, "mHClean/I");
  int mLClean = 0;
  skimTree->Branch("mLClean", &mLClean, "mLClean/I");
  vector<double> aenorm;
  skimTree->Branch("aenorm", &aenorm);
  vector<double> avse;
  skimTree->Branch("avse", &avse);
//  vector<double> t150;
  vector<double> kvorrT;
//  vector<double> toe;
//  vector<double> d2wfnoiseTagNorm;
  vector<double> aenorm85;
  if(!smallOutput){
 //   skimTree->Branch("t150", &t150);
    skimTree->Branch("kvorrT", &kvorrT);
  }
  vector<double> rawDCR;
  vector<double> dcrSlope85;
  vector<double> dcr90;
  vector<double> dcrSlope95;
  vector<double> dcrSlope98;
  vector<double> dcrSlope99;
  vector<double> dcrctc90;
  vector<double> nlcblrwfSlope;

 // data cleaning variables
  unsigned int eventDC1Bits = 0;
  skimTree->Branch("EventDC1Bits", &eventDC1Bits, "eventDC1Bits/i");
  vector<unsigned int> wfDCBits;
  skimTree->Branch("wfDCBits", &wfDCBits);
  vector<double> lnFillTimes1;
  vector<double> lnFillTimes2;
  vector<bool> isLNFill1;
  skimTree->Branch("isLNFill1", &isLNFill1);
  vector<bool> isLNFill2;
  skimTree->Branch("isLNFill2", &isLNFill2);
  vector<int> nX;
  skimTree->Branch("nX", &nX);
  vector<double> trapETailMin;
  if(!smallOutput) skimTree->Branch("trapETailMin", &trapETailMin);

  // veto variables
  vector<double> dtmu_s;
  skimTree->Branch("dtmu_s", &dtmu_s);
  vector<int> muType;
  skimTree->Branch("muType", &muType);
  vector<bool> muTUnc;
  skimTree->Branch("muTUnc", &muTUnc);
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
  if (dsNumber == 2) {
    detIDIsBad[1426981] = true;
    detIDIsBad[1426622] = true;
    detIDIsBad[28455] = true;
    detIDIsBad[28470] = true;
    detIDIsBad[28463] = true;
    detIDIsBad[28465] = true;
    detIDIsBad[28469] = true;
    detIDIsBad[28477] = true;
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
  if(dsNumber == 5) {
    detIDIsBad[1426981] = true;
    detIDIsBad[1426622] = true;
    detIDIsBad[28477] = true;
    detIDIsBad[1425731] = true;
    detIDIsBad[1426611] = true;
    detIDIsBad[28595] = true;
    detIDIsBad[28461] = true;
    detIDIsBad[1428530] = true;
    detIDIsBad[28621] = true;
    detIDIsBad[28473] = true;
    detIDIsBad[1426651] = true;
    detIDIsBad[1429092] = true;
    detIDIsBad[1426652] = true;
    detIDIsBad[28619] = true;
    detIDIsBad[1427121] = true;
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
  if (dsNumber == 2) {
    detIDIsVetoOnly[28480] = true;
    detIDIsVetoOnly[1425751] = true;
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
 if(dsNumber == 5) {
    detIDIsVetoOnly[28480] = true;
    detIDIsVetoOnly[28470] = true;
    detIDIsVetoOnly[28463] = true;
    detIDIsVetoOnly[28459] = true;
    detIDIsVetoOnly[1426641] = true;
    detIDIsVetoOnly[1427481] = true;
    detIDIsVetoOnly[28456] = true;
    detIDIsVetoOnly[1427120] = true;
  }

  // start loop over all events
  while(gatReader.Next()) {

    // stuff to do on run boundaries
    if(runSave != *runIn) {
      runSave = *runIn;
      cout << "Processing run " << *runIn << ", "
           << skimTree->GetEntries() << " entries saved so far"
           << endl;
      skimTree->Write("", TObject::kOverwrite);
    }

    // Skip this event if it is a pulser event as identified by Pinghan
    if(*eventDC1BitsIn & kPinghanPulserMask) continue;
    iEvent = gatChain->GetTree()->GetReadEntry();

    // copy the event-level info to the output fields
    gatrev = *gatrevIn;
    run = int(*runIn);
    startTime = *startTimeIn;
    stopTime = *stopTimeIn;
    eventDC1Bits = *eventDC1BitsIn;
    sumEH = 0;
    sumEL = 0;
    mH = 0;
    mL = 0;
    sumEHClean = 0;
    sumELClean = 0;
    mHClean = 0;
    mLClean = 0;

    // clear all hit-level info fields
    iHit.resize(0);
    if(!smallOutput){
      trapECal.resize(0);
      onBoardE.resize(0);
   //   t150.resize(0);
      kvorrT.resize(0);
      aenorm85.resize(0);
   //   toe.resize(0);
   //   d2wfnoiseTagNorm.resize(0);
      trapETailMin.resize(0);
    }
    trapENFCal.resize(0);
    trapENMCal.resize(0);
    channel.resize(0);
    tloc_s.resize(0);
    time_s.resize(0);
    if(!simulatedInput)
    {
      timeMT.resize(0);
      dateMT.resize(0);
    }
    pos.resize(0);
    det.resize(0);
    cryo.resize(0);
    gain.resize(0);
    mageID.resize(0);
    detID.resize(0);
    detName.resize(0);
    isEnr.resize(0);
    isNat.resize(0);
    mAct_g.resize(0);
    isGood.resize(0);
    wfDCBits.resize(0);
    aenorm.resize(0);
    avse.resize(0);
    isLNFill1.resize(0);
    isLNFill2.resize(0);
    nX.resize(0);
    if (!simulatedInput)
    {
      muVeto.resize(0);
      muType.resize(0);
      muTUnc.resize(0);
      dtmu_s.resize(0);
    }

    // loop over hits
    bool skipMe = false;
    size_t nHits = trapENFCalIn->size();
    for(size_t i=0; i<nHits; i++) {

      // skip all hits with E_H < 2 keV, E_L < 10 keV in -both- trapE and trapENF
      // For small skim files, skip all hits with E_H and E_L < 200 keV in trapE and trapENF
      double hitENFCal = (*trapENFCalIn)[i];
      double hitENMCal = (*trapENMCalIn)[i];
      // double hitENF = (*trapENFIn)[i];
      double hitEMax = (*trapECalIn)[i];
      // double hitTrapMax = hitENMCal;
      int hitCh = (*channelIn)[i];
      if(!smallOutput && hitCh%2 == 0 && hitENFCal < 2000 && hitEMax < 2000) continue;  // set for veto
      if(!smallOutput && hitCh%2 == 1 && hitENFCal < 2000 && hitEMax < 2000) continue;
      if(smallOutput && hitCh%2 == 0 && hitENFCal <  200. && hitEMax < 200.) continue;
      if(smallOutput && hitCh%2 == 1 && hitENFCal < 200. && hitEMax < 200.) continue;
      // skip hits from totally "bad" detectors (not biased, etc), or from
      // use-for-veto-only detectors if E < 10 keV
      int hitDetID = (*detIDIn)[i];
      if(!simulatedInput)
      {
        if(detIDIsBad[hitDetID] || (detIDIsVetoOnly[hitDetID] && hitEMax < 10.)) continue;
      }
      // copy over hit info
      iHit.push_back(i);
      trapENFCal.push_back(hitENFCal);
      trapENMCal.push_back(hitENMCal);
      channel.push_back(hitCh);
      double hitT_s = (*timestampIn)[i]*1.e-8;
      tloc_s.push_back(hitT_s);
      cryo.push_back((*cryoIn)[i]);
      time_s.push_back( (startTime - startTime0) + hitT_s ); //Need to figure out what to do with continuous running, Clara 10/10/16
      if(!simulatedInput)
      {
        timeMT.push_back((*(*timeMTIn))[i]);
        dateMT.push_back((*(*dateMTIn))[i]);
      }
      pos.push_back((*posIn)[i]);
      det.push_back((*detIn)[i]);
      gain.push_back(hitCh % 2);

      mageID.push_back((*mageIDIn)[i]);
      detID.push_back(hitDetID);
      detName.push_back((*detNameIn)[i]);
      isEnr.push_back((*detNameIn)[i][0] == 'P');
      isNat.push_back((*detNameIn)[i][0] == 'B');
      mAct_g.push_back(actM4Det_g[hitDetID]);
      isGood.push_back(!detIDIsVetoOnly[hitDetID]);
      wfDCBits.push_back((*wfDCBitsIn)[i]);
     // d2wfnoiseTagNorm.push_back((*d2wfnoiseTagNormIn)[i]);
      nX.push_back((*nRisingXIn)[i]);
      if(!smallOutput){
        trapECal.push_back(hitEMax);
        onBoardE.push_back((*energyIn)[i]);
//        double t1 = (*blrwfFMR1In)[i];
//        double t50 = (*blrwfFMR50In)[i];
//        t150.push_back(t50-t1);
        kvorrT.push_back((*triTrapMaxIn)[i]);
//        toe.push_back((*toeIn)[i] / hitEMax);
        trapETailMin.push_back((*trapETailMinIn)[i]);
      }

      // sum energies and multiplicities
      if(hitCh%2 == 0) {
        mH++;
        if(!detIDIsVetoOnly[hitDetID]) sumEH += hitENFCal;
      }
      else {
        mL++;
        if(!detIDIsVetoOnly[hitDetID]) sumEL += hitENFCal;
      }
      if(hitCh%2 == 0 && ~((~0x010) & wfDCBits[wfDCBits.size()-1])) {
        mHClean++;
        if(!detIDIsVetoOnly[hitDetID]) sumEHClean += hitENFCal;
      }
      if (hitCh%2 == 1 && ~((~0x010) & wfDCBits[wfDCBits.size()-1])) {
        mLClean++;
        if(!detIDIsVetoOnly[hitDetID]) sumELClean += hitENFCal;
      }

      if(!simulatedInput)
      {
        // Find the most recent muon to this event
        while(1)
        {
          if(iMu >= nMu-1) break;
          double tmuUnc = 1.e-8; // normally 10ns uncertainty
          if (muUncert[iMu+1] > tmuUnc) tmuUnc = muUncert[iMu+1];
          if (muRuns[iMu+1] > run) break;
          else if (muRuns[iMu+1]==run && (muTimes[iMu+1]-tmuUnc) > hitT_s) break;
          // printf("Inc:iMu+1 %-4lu  gRun %-4i  mRun %-4i  tGe %-8.3f  tMu %-8.3f  dtRun %-8.0f  dtEvent %-8.0f\n" ,iMu+1, run, muRuns[iMu+1], hitT_s, muTimes[iMu], muRunTStarts[iMu+1]-startTime, (muTimes[iMu+1] - tmuUnc) - hitT_s);
          iMu++;
        }

        // Calculate time since last muon.
        // NOTE: If there has been a clock reset since the last muon hit, this delta-t will be incorrect.
        double dtmu = 0;
        if (dsNumber==0)
          dtmu = (startTime - muRunTStarts[iMu]) + (hitT_s - muTimes[iMu]);
        else
          dtmu = (hitT_s - muTimes[iMu]);

        bool vetoThisHit = (dtmu > -1.*(muUncert[iMu]) && dtmu < (1. + muUncert[iMu]));

        // DS-4 requires a larger window due to synchronization issues.
        if (dsNumber==4)
          vetoThisHit = (dtmu > -3.*(muUncert[iMu]) && dtmu < (4. + muUncert[iMu]));

        // if (vetoThisHit) printf("Coin: iMu %-4lu  det %i  gRun %-4i  mRun %-5i  tGe %-7.3f  tMu %-7.3f  ene %-6.0f  veto? %i  dtmu %.2f +/- %.2f\n", iMu,hitCh,run,muRuns[iMu],hitT_s,muTimes[iMu],hitENFCal,vetoThisHit,dtmu,muUncert[iMu]);

        dtmu_s.push_back(dtmu);
        muType.push_back(muTypes[iMu]);
        muTUnc.push_back(muUncert[iMu]);
        muVeto.push_back(vetoThisHit);


      }
    }

    if(!simulatedInput)
    {
      // If no good hits in the event or skipped for some other reason, don't
      // write this event to the output tree.
      if(trapENFCal.size() == 0 || skipMe) continue;
    }

    // finally, fill the tree for this event
    skimTree->Fill();
  }

  // write output tree to output file
  cout << "Closing out skim file..." << endl;
  skimTree->Write("", TObject::kOverwrite);
  fOut->Close();
  return 0;
}

void LoadRun(GATDataSet& ds, size_t i)
{
  ds.AddRunNumber(i);
}

void LoadActiveMasses(map<int,double>& activeMassForDetID_g, int dsNumber)
{
  if(dsNumber == 0 || dsNumber == 1 || dsNumber == 3) {
    activeMassForDetID_g[1426981] = 509.9;
    activeMassForDetID_g[1425750] = 978.8;
    activeMassForDetID_g[1426612] = 811.3;
    activeMassForDetID_g[1425380] = 967.9;
    activeMassForDetID_g[28474] = 587.0;
    activeMassForDetID_g[1426640] = 722.9;
    activeMassForDetID_g[1426650] = 659.1;
    activeMassForDetID_g[1426622] = 688.6;
    activeMassForDetID_g[28480] = 577.6;
    activeMassForDetID_g[1426980] = 886.3;
    activeMassForDetID_g[1425381] = 949.0;
    activeMassForDetID_g[1425730] = 1023.8;
    activeMassForDetID_g[28455] = 584.2;
    activeMassForDetID_g[28470] = 590.8;
    activeMassForDetID_g[28463] = 594.5;
    activeMassForDetID_g[28465] = 571.1;
    activeMassForDetID_g[28469] = 593.3;
    activeMassForDetID_g[28477] = 579.5;
    activeMassForDetID_g[1425751] = 730.6;
    activeMassForDetID_g[1426610] = 632.0;
    activeMassForDetID_g[1425731] = 982.0;
    activeMassForDetID_g[1425742] = 731.6;
    activeMassForDetID_g[1426611] = 675.5;
    activeMassForDetID_g[1425740] = 701.4;
    activeMassForDetID_g[1426620] = 572.3;
    activeMassForDetID_g[28482] = 588.0;
    activeMassForDetID_g[1425741] = 709.9;
    activeMassForDetID_g[1426621] = 590.9;
    activeMassForDetID_g[1425370] = 964.3;
  }
 else if(dsNumber == 4) {
    activeMassForDetID_g[28459] = 556.0;
    activeMassForDetID_g[1426641] = 576.0;
    activeMassForDetID_g[1427481] = 903.0;
    activeMassForDetID_g[1427480] = 917.0;
    activeMassForDetID_g[28481] = 581.0;
    activeMassForDetID_g[28576] = 562.0;
    activeMassForDetID_g[28594] = 559.0;
    activeMassForDetID_g[28595] = 558.0;
    activeMassForDetID_g[28461] = 557.0;
    activeMassForDetID_g[1427490] = 872.0;
    activeMassForDetID_g[1427491] = 852.0;
    activeMassForDetID_g[1428530] = 996.0;
    activeMassForDetID_g[28607] = 558.0;
    activeMassForDetID_g[28456] = 579.0;
    activeMassForDetID_g[28621] = 565.0;
    activeMassForDetID_g[28466] = 566.0;
    activeMassForDetID_g[28473] = 562.0;
    activeMassForDetID_g[28487] = 557.0;
    activeMassForDetID_g[1426651] = 591.0;
    activeMassForDetID_g[1428531] = 1031.0;
    activeMassForDetID_g[1427120] = 802.0;
    activeMassForDetID_g[1235170] = 462.2;
    activeMassForDetID_g[1429091] = 775.0;
    activeMassForDetID_g[1429092] = 821.0;
    activeMassForDetID_g[1426652] = 778.0;
    activeMassForDetID_g[28619] = 566.0;
    activeMassForDetID_g[1427121] = 968.0;
    activeMassForDetID_g[1429090] = 562.0;
    activeMassForDetID_g[28717] = 567.0;
}
  else if(dsNumber == 5) {
    activeMassForDetID_g[1426981] = 509.9;
    activeMassForDetID_g[1425750] = 978.8;
    activeMassForDetID_g[1426612] = 811.3;
    activeMassForDetID_g[1425380] = 967.9;
    activeMassForDetID_g[28474] = 587.0;
    activeMassForDetID_g[1426640] = 722.9;
    activeMassForDetID_g[1426650] = 659.1;
    activeMassForDetID_g[1426622] = 688.6;
    activeMassForDetID_g[28480] = 577.6;
    activeMassForDetID_g[1426980] = 886.3;
    activeMassForDetID_g[1425381] = 949.0;
    activeMassForDetID_g[1425730] = 1023.8;
    activeMassForDetID_g[28455] = 584.2;
    activeMassForDetID_g[28470] = 590.8;
    activeMassForDetID_g[28463] = 594.5;
    activeMassForDetID_g[28465] = 571.1;
    activeMassForDetID_g[28469] = 593.3;
    activeMassForDetID_g[28477] = 579.5;
    activeMassForDetID_g[1425751] = 730.6;
    activeMassForDetID_g[1426610] = 632.0;
    activeMassForDetID_g[1425731] = 982.0;
    activeMassForDetID_g[1425742] = 731.6;
    activeMassForDetID_g[1426611] = 675.5;
    activeMassForDetID_g[1425740] = 701.4;
    activeMassForDetID_g[1426620] = 572.3;
    activeMassForDetID_g[28482] = 588.0;
    activeMassForDetID_g[1425741] = 709.9;
    activeMassForDetID_g[1426621] = 590.9;
    activeMassForDetID_g[1425370] = 964.3;
    activeMassForDetID_g[28459] = 556.0;
    activeMassForDetID_g[1426641] = 576.0;
    activeMassForDetID_g[1427481] = 903.0;
    activeMassForDetID_g[1427480] = 917.0;
    activeMassForDetID_g[28481] = 581.0;
    activeMassForDetID_g[28576] = 562.0;
    activeMassForDetID_g[28594] = 559.0;
    activeMassForDetID_g[28595] = 558.0;
    activeMassForDetID_g[28461] = 557.0;
    activeMassForDetID_g[1427490] = 872.0;
    activeMassForDetID_g[1427491] = 852.0;
    activeMassForDetID_g[1428530] = 996.0;
    activeMassForDetID_g[28607] = 558.0;
    activeMassForDetID_g[28456] = 579.0;
    activeMassForDetID_g[28621] = 565.0;
    activeMassForDetID_g[28466] = 566.0;
    activeMassForDetID_g[28473] = 562.0;
    activeMassForDetID_g[28487] = 557.0;
    activeMassForDetID_g[1426651] = 591.0;
    activeMassForDetID_g[1428531] = 1031.0;
    activeMassForDetID_g[1427120] = 802.0;
    activeMassForDetID_g[1235170] = 462.2;
    activeMassForDetID_g[1429091] = 775.0;
    activeMassForDetID_g[1429092] = 821.0;
    activeMassForDetID_g[1426652] = 778.0;
    activeMassForDetID_g[28619] = 566.0;
    activeMassForDetID_g[1427121] = 968.0;
    activeMassForDetID_g[1429090] = 562.0;
    activeMassForDetID_g[28717] = 567.0;
}

  else cout << "LoadActiveMasses(): unknown dataset number DS" << dsNumber << endl;
}
