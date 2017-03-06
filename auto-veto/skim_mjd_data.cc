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
double GetAvsE(int channel, double TSCurrent50nsMax, double TSCurrent100nsMax, double TSCurrent200nsMax, double trapENF, double trapENFCal, int dsNumber);
void LoadLNFillTimes1(vector<double>& lnFillTimes1, int dsNumber);
void LoadLNFillTimes2(vector<double>& lnFillTimes2, int dsNumber);
double GetDCRraw(int channel, double nlcblrwfSlope, double trapMax, int dsNumber);
double GetDCR85(int channel, double nlcblrwfSlope, double trapMax, int dsNumber);
double GetDCR90(int channel, double nlcblrwfSlope, double trapMax, int dsNumber);
double GetDCR95(int channel, double nlcblrwfSlope, double trapMax, int dsNumber);
double GetDCR98(int channel, double nlcblrwfSlope, double trapMax, int dsNumber);
double GetDCR99(int channel, double nlcblrwfSlope, double trapMax, int dsNumber);
double GetDCRCTC(int channel, double nlcblrwfSlope, double trapE, double trapMax, int dsNumber);

int main(int argc, const char** argv)
{
  if(argc < 3 || argc > 9) {
    cout << "To include tail slope add flag -s. For raw DCR add flag -r " << endl;
    cout << "For minimal skim file add flag -m " << endl;
    cout << "For extensive skim file (multiple DCR and aenorm) add flag -e " << endl;
    cout << "For custom energy threshold: -t [number (default is 2 keV)]" << endl;
    cout << "Usage for single run: " << argv[0] << " -f [runNum] (output path)" << endl;
    cout << "Usage for custom file: " << argv[0] << " --filename [filename] [runNum] (output path)" << endl;
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
  bool singleFile = false;
  bool writeRawDCR = false;
  bool writeSlope = false;
  bool smallOutput = false;
  bool extendedOutput = false;
  bool simulatedInput = false;
  vector<string> args;
  for(int iArg=0; iArg<argc; ++iArg) args.push_back(argv[iArg]);

  auto slopeArg = find(args.begin(), args.end(), "-s");
  if(slopeArg!=args.end()){
    writeSlope = true;
    cout<<"Tail slope option selected."<<endl;
    args.erase(slopeArg);
  }
  auto rawDCRArg = find(args.begin(), args.end(), "-r");
  if(rawDCRArg!=args.end()){
    writeRawDCR = true;
    cout<<"Raw DCR option selected."<<endl;
    args.erase(rawDCRArg);
  }
  auto minimalSkimArg = find(args.begin(), args.end(), "-m");
  if(minimalSkimArg!=args.end()){
    smallOutput = true;
    cout<<"Minimal skim file option selected."<<endl;
    args.erase(minimalSkimArg);
  }
  auto extendedOutputArg = find(args.begin(), args.end(), "-e");
  if(extendedOutputArg!=args.end()){
    extendedOutput = true;
    cout<<"Extended skim file option selected."<<endl;
    args.erase(extendedOutputArg);
  }
  auto energyThreshArg = find(args.begin(), args.end(), "-t");
  if(energyThreshArg!=args.end()){
    int pos = find(args.begin(), args.end(), "-t") - args.begin();
    energyThresh = stod(args[pos+1]);
    cout << "Set HG energy threshold to " << energyThresh << " keV\n";
  }

  auto fileArg = find(args.begin(), args.end(), "-f");
  auto fileNameArg = find(args.begin(), args.end(), "--filename");
  if(fileArg != args.end())
  {
    try
    {
      singleFile = true;
      runSeq = atoi((fileArg+1)->c_str());
      args.erase(fileArg, fileArg+2);
      if(runSeq >= 2335 && runSeq <= 8183) dsNumber = 0;
        else if(runSeq >= 8722 && runSeq < 14502) dsNumber = 1;
        else if(runSeq >= 14503 && runSeq < 15892) dsNumber = 2;
        else if(runSeq >= 16289 && runSeq < 18622) dsNumber = 3;
        else if(runSeq >= 60000549 && runSeq < 70000000) dsNumber = 4;
	else if(runSeq >= 18623 && runSeq <= 30000) dsNumber = 5;
        else {
          cout << "Error: I don't know what dataset run " << runSeq << " is from." << endl;
          return 1;
        }
      cout << "loading run " << runSeq << endl;
      LoadRun(ds, runSeq);
    }
    catch(exception& e)
    {
      cerr << "Invalid command-line argument." << endl;
      cerr << "Terminating skim_mjd_data." << endl;
      return 1;
    }
  }
  else if(fileNameArg != args.end())
  {
    try
    {
      singleFile = true;
      runSeq = atoi((fileNameArg+2)->c_str());
      string fileName = *(fileNameArg+1);
      args.erase(fileNameArg, fileNameArg+3);
      if(runSeq >= 2335 && runSeq <= 8183) dsNumber = 0;
        else if(runSeq >= 8722 && runSeq < 14502) dsNumber = 1;
        else if(runSeq >= 14503 && runSeq < 15892) dsNumber = 2;
        else if(runSeq >= 16289 && runSeq < 18622) dsNumber = 3;
        else if(runSeq >= 60000549 && runSeq < 70000000) dsNumber = 4;
	else if(runSeq >= 18623 && runSeq <= 30000) dsNumber = 5;
        else {
          cout << "Error: I don't know what dataset run " << runSeq << " is from." << endl;
          return 1;
        }
      cout << "Reading file " << fileName << " as run " << runSeq << endl;
      gatChain = new TChain("mjdTree","mjdTree");
      gatChain->AddFile(fileName.c_str());
      vetoChain = new TChain("vetoTree","vetoTree");
      string vetoPath = ds.GetPathToRun(runSeq,GATDataSet::kVeto);
      vetoChain->Add(vetoPath.c_str());
    }
    catch(exception& e)
    {
      cerr << "Invalid command-line argument." << endl;
      cerr << "Terminating skim_mjd_data." << endl;
      return 1;
    }
  }
  else
  {
    try
    {
      dsNumber = stoi(args[1]);
      runSeq = stoi(args[2]);
      args.erase(args.begin()+1, args.begin()+3);
      cout << "loading dataset " << dsNumber << " run sequence " << runSeq << endl;
      LoadDataSet(ds, dsNumber, runSeq);
    }
    catch(exception& e)
    {
      cerr << "Invalid command-line argument." << endl;
      cerr << "Terminating skim_mjd_data." << endl;
      return 1;
    }
  }

  if(args.size() > 1) outputPath += args[1];

  // set up dataset
  if(gatChain==NULL) gatChain = ds.GetGatifiedChain(false);
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
  if(singleFile) filename += "run";
  filename += runStr;
  if(writeSlope) filename += "_slope";
  if(writeRawDCR) filename += "_rawDCR";
  if(smallOutput) filename += "_small";
  if(extendedOutput) filename += "_ext";
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
 //   skimTree->Branch("toe", &toe);
 //   skimTree->Branch("d2wfnoiseTagNorm", &d2wfnoiseTagNorm);
    if(extendedOutput) skimTree->Branch("aenorm85", &aenorm85);
  }
  vector<double> rawDCR;
  vector<double> dcrSlope85;
  vector<double> dcr90;
  vector<double> dcrSlope95;
  vector<double> dcrSlope98;
  vector<double> dcrSlope99;
  vector<double> dcrctc90;
  vector<double> nlcblrwfSlope;
  if(writeSlope) skimTree->Branch("nlcblrwfSlope", &nlcblrwfSlope);
  if(writeRawDCR) skimTree->Branch("rawDCR", &rawDCR);
  else{
    if(!smallOutput && extendedOutput) {
      skimTree->Branch("dcrSlope85", &dcrSlope85);
      skimTree->Branch("dcrSlope95", &dcrSlope95);
      skimTree->Branch("dcrSlope98", &dcrSlope98);
      skimTree->Branch("dcrSlope99", &dcrSlope99);
    }
    skimTree->Branch("dcr90", &dcr90);
    skimTree->Branch("dcrctc90", &dcrctc90);
  }
 // data cleaning variables
  unsigned int eventDC1Bits = 0;
  skimTree->Branch("EventDC1Bits", &eventDC1Bits, "eventDC1Bits/i");
  vector<unsigned int> wfDCBits;
  skimTree->Branch("wfDCBits", &wfDCBits);
  vector<double> lnFillTimes1;
  LoadLNFillTimes1(lnFillTimes1, dsNumber);
  vector<double> lnFillTimes2;
  LoadLNFillTimes2(lnFillTimes2, dsNumber);
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
    if(writeSlope) nlcblrwfSlope.resize(0);
    if(writeRawDCR) rawDCR.resize(0);
    else {
      if(!smallOutput){
        dcrSlope85.resize(0);
        dcrSlope95.resize(0);
        dcrSlope98.resize(0);
      }
      dcr90.resize(0);
      dcrctc90.resize(0);
      dcrSlope99.resize(0);
    }
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
      double hitENF = (*trapENFIn)[i];
      double hitEMax = (*trapECalIn)[i];
      double hitTrapMax = hitENMCal;
      int hitCh = (*channelIn)[i];
      if(!smallOutput && hitCh%2 == 0 && hitENFCal < energyThresh && hitEMax < energyThresh) continue;
      if(!smallOutput && hitCh%2 == 1 && hitENFCal < 10. && hitEMax < 10.) continue;
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
      double a50 = (*tsCurrent50nsMaxIn)[i];
      double a100 = (*tsCurrent100nsMaxIn)[i];
      double a200 = (*tsCurrent200nsMaxIn)[i];
      avse.push_back(GetAvsE(hitCh, a50, a100, a200, hitENF, hitENFCal, dsNumber));
      if(writeSlope) nlcblrwfSlope.push_back((*dcrSlopeIn)[i]);
      if(writeRawDCR) rawDCR.push_back(GetDCRraw(hitCh, (*dcrSlopeIn)[i], hitTrapMax, dsNumber));
      else {
        dcr90.push_back(GetDCR90(hitCh, (*dcrSlopeIn)[i], hitTrapMax, dsNumber));
        dcrctc90.push_back(GetDCRCTC(hitCh, (*dcrSlopeIn)[i], hitENFCal, hitTrapMax, dsNumber));
        if(!smallOutput){
          dcrSlope85.push_back(GetDCR85(hitCh, (*dcrSlopeIn)[i], hitTrapMax, dsNumber));
          dcrSlope95.push_back(GetDCR95(hitCh, (*dcrSlopeIn)[i], hitTrapMax, dsNumber));
          dcrSlope98.push_back(GetDCR98(hitCh, (*dcrSlopeIn)[i], hitTrapMax, dsNumber));
        //dcrSlope99.push_back(GetDCR99(hitCh, (*dcrSlopeIn)[i], hitTrapMax, dsNumber));
        }
      }

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

        // tag LN fills
        double utctime = startTime + hitT_s;
        isLNFill1.push_back(false);
        isLNFill2.push_back(false);
        for(size_t i=0; i<lnFillTimes1.size(); i++) {
          if(lnFillTimes1[i]+300 < utctime) continue;
          if(lnFillTimes1[i]-900 > utctime) break;
          isLNFill1[isLNFill1.size()-1] = true;
          break;
        for(size_t i=0; i<lnFillTimes2.size(); i++) {
          if(lnFillTimes2[i]+300 < utctime) continue;
          if(lnFillTimes2[i]-900 > utctime) break;
          isLNFill2[isLNFill2.size()-1] = true;
          break;
        }
        }
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

double GetAvsE(int channel, double TSCurrent50nsMax, double TSCurrent100nsMax, double TSCurrent200nsMax, double trapENF, double trapENFCal, int dsNumber)
{
  if(dsNumber == 1) {

    if    (channel==582) return (-1*((TSCurrent200nsMax*(trapENFCal)/(trapENF))-(0.01854576347230)-(0.00468276622549*(trapENFCal))-(-0.00000001407301*(trapENFCal)*(trapENFCal)))/-0.0702263);
    if    (channel==583) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.16369855142045)-(0.00628861548784*(trapENFCal))-(0.00000009715465*(trapENFCal)*(trapENFCal)))/-0.208138);
    if    (channel==580) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(0.05674199580730)-(0.00425035340190*(trapENFCal))-(0.00000000171536*(trapENFCal)*(trapENFCal)))/-0.118208);
    if    (channel==581) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.32060502828949)-(0.00423312877746*(trapENFCal))-(0.00000014226922*(trapENFCal)*(trapENFCal)))/-0.222847);
    if    (channel==578) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.08927516275763)-(0.00703937089571*(trapENFCal))-(0.00000003604893*(trapENFCal)*(trapENFCal)))/-0.248755);
    if    (channel==579) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.29254857541595)-(0.00689197923341*(trapENFCal))-(0.00000015350572*(trapENFCal)*(trapENFCal)))/-0.216434);
    if    (channel==692) return (-1*((TSCurrent200nsMax*(trapENFCal)/(trapENF))-(0.06259010350846)-(0.00483586165925*(trapENFCal))-(0.00000001627416*(trapENFCal)*(trapENFCal)))/-0.0820596);
    if    (channel==693) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.12134467443495)-(0.00691895832301*(trapENFCal))-(-0.00000000470471*(trapENFCal)*(trapENFCal)))/-0.224861);
    if    (channel==648) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.06437389210040)-(0.00719101008714*(trapENFCal))-(-0.00000005202873*(trapENFCal)*(trapENFCal)))/-0.504072);
    if    (channel==649) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.19367715373649)-(0.00704134129392*(trapENFCal))-(0.00000000395276*(trapENFCal)*(trapENFCal)))/-0.47942);
    if    (channel==640) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.10070555974961)-(0.00742514512250*(trapENFCal))-(0.00000002369359*(trapENFCal)*(trapENFCal)))/-0.387218);
    if    (channel==641) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(0.07526220050901)-(0.00667997978625*(trapENFCal))-(0.00000002407041*(trapENFCal)*(trapENFCal)))/-0.270941);
    if    (channel==616) return (-1*((TSCurrent200nsMax*(trapENFCal)/(trapENF))-(0.02486114215672)-(0.00455738079802*(trapENFCal))-(-0.00000010811773*(trapENFCal)*(trapENFCal)))/-0.129789);
    if    (channel==617) return (-1*((TSCurrent200nsMax*(trapENFCal)/(trapENF))-(0.02787931415456)-(0.00458580492241*(trapENFCal))-(-0.00000012117977*(trapENFCal)*(trapENFCal)))/-0.124577);
    if    (channel==610) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.05649949694814)-(0.00716511550872*(trapENFCal))-(-0.00000003055894*(trapENFCal)*(trapENFCal)))/-0.199234);
    if    (channel==611) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.14399987715688)-(0.00714778413265*(trapENFCal))-(0.00000008102384*(trapENFCal)*(trapENFCal)))/-0.310912);
    if    (channel==608) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.05524811515552)-(0.00743711815084*(trapENFCal))-(-0.00000003756238*(trapENFCal)*(trapENFCal)))/-0.281594);
    if    (channel==609) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.15876502286021)-(0.00732122138428*(trapENFCal))-(0.00000004729350*(trapENFCal)*(trapENFCal)))/-0.430407);
    if    (channel==664) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.07779551311681)-(0.00607567539278*(trapENFCal))-(-0.00000004657088*(trapENFCal)*(trapENFCal)))/-0.250382);
    if    (channel==665) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(0.07304968321561)-(0.00564943039142*(trapENFCal))-(-0.00000008192775*(trapENFCal)*(trapENFCal)))/-0.182482);
    if    (channel==672) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.06848460163308)-(0.00725365163400*(trapENFCal))-(-0.00000002128673*(trapENFCal)*(trapENFCal)))/-0.290611);
    if    (channel==673) return (-1*((TSCurrent200nsMax*(trapENFCal)/(trapENF))-(0.05831439575329)-(0.00500538585273*(trapENFCal))-(0.00000000830948*(trapENFCal)*(trapENFCal)))/-0.137654);
    if    (channel==632) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.07917762952496)-(0.00694792043183*(trapENFCal))-(-0.00000000758383*(trapENFCal)*(trapENFCal)))/-0.181945);
    if    (channel==633) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(0.06954722595344)-(0.00626199652534*(trapENFCal))-(-0.00000001155030*(trapENFCal)*(trapENFCal)))/-0.133663);
    if    (channel==626) return (-1*((TSCurrent200nsMax*(trapENFCal)/(trapENF))-(0.01810678070787)-(0.00488771705012*(trapENFCal))-(-0.00000000804786*(trapENFCal)*(trapENFCal)))/-0.086283);
    if    (channel==627) return (-1*((TSCurrent200nsMax*(trapENFCal)/(trapENF))-(0.02984410185282)-(0.00496024053732*(trapENFCal))-(-0.00000000937571*(trapENFCal)*(trapENFCal)))/-0.102878);
    if    (channel==690) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(0.02821723031971)-(0.00621139001811*(trapENFCal))-(0.00000000902572*(trapENFCal)*(trapENFCal)))/-0.154616);
    if    (channel==691) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(0.06022352707069)-(0.00629275630451*(trapENFCal))-(0.00000005072433*(trapENFCal)*(trapENFCal)))/-0.152819);
    if    (channel==600) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.10055600376558)-(0.00656624540158*(trapENFCal))-(0.00000001755970*(trapENFCal)*(trapENFCal)))/-0.195797);
    if    (channel==601) return (-1*((TSCurrent200nsMax*(trapENFCal)/(trapENF))-(0.07060566139558)-(0.00462456482179*(trapENFCal))-(0.00000003567736*(trapENFCal)*(trapENFCal)))/-0.059421);
    if    (channel==598) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.08363810522446)-(0.00697336661343*(trapENFCal))-(-0.00000004853185*(trapENFCal)*(trapENFCal)))/-0.174753);
    if    (channel==599) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(0.10050448277932)-(0.00620078423930*(trapENFCal))-(-0.00000000963453*(trapENFCal)*(trapENFCal)))/-0.132879);
    if    (channel==594) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(0.03637253831151)-(0.00614898547205*(trapENFCal))-(0.00000002174726*(trapENFCal)*(trapENFCal)))/-0.16588);
    if    (channel==595) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.13382021226886)-(0.00684308040220*(trapENFCal))-(0.00000009806745*(trapENFCal)*(trapENFCal)))/-0.2612);
    if    (channel==592) return (-1*((TSCurrent200nsMax*(trapENFCal)/(trapENF))-(0.00308216986982)-(0.00513834129293*(trapENFCal))-(-0.00000006730542*(trapENFCal)*(trapENFCal)))/-0.231685);
    if    (channel==593) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.10854259868427)-(0.00695707196683*(trapENFCal))-(0.00000002953986*(trapENFCal)*(trapENFCal)))/-0.961633);
}

  else if(dsNumber == 3) {
    if    (channel==582) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.03291013810328)-(0.00592019288408*(trapENFCal))-(-0.00000003872612*(trapENFCal)*(trapENFCal)))/-0.0519063);
    if    (channel==580) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.04157771946995)-(0.00436286247927*(trapENFCal))-(-0.00000003049836*(trapENFCal)*(trapENFCal)))/-0.0226936);
    if    (channel==581) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.00395765626708)-(0.00433733122891*(trapENFCal))-(-0.00000001661893*(trapENFCal)*(trapENFCal)))/-0.0362425);
    if    (channel==578) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.05751673378761)-(0.00639392239166*(trapENFCal))-(-0.00000003439059*(trapENFCal)*(trapENFCal)))/-0.0607365);
    if    (channel==579) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.05670651469833)-(0.00652829078356*(trapENFCal))-(-0.00000003014225*(trapENFCal)*(trapENFCal)))/-0.0544917);
    if    (channel==692) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.05738334053388)-(0.00631201007769*(trapENFCal))-(-0.00000003718222*(trapENFCal)*(trapENFCal)))/-0.0552909);
    if    (channel==693) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.03127857380050)-(0.00630596961348*(trapENFCal))-(-0.00000001610488*(trapENFCal)*(trapENFCal)))/-0.0796882);
    if    (channel==648) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.04437940443903)-(0.00644791547964*(trapENFCal))-(-0.00000005788124*(trapENFCal)*(trapENFCal)))/-0.176498);
    if    (channel==649) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(0.00653881959029)-(0.00639808779261*(trapENFCal))-(-0.00000001324761*(trapENFCal)*(trapENFCal)))/-0.205588);
    if    (channel==640) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.04654487223408)-(0.00668445021629*(trapENFCal))-(-0.00000005487387*(trapENFCal)*(trapENFCal)))/-0.165336);
    if    (channel==641) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.08479370822675)-(0.00681362022520*(trapENFCal))-(-0.00000009047305*(trapENFCal)*(trapENFCal)))/-0.184052);
    if    (channel==610) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.06949704676841)-(0.00640494467551*(trapENFCal))-(-0.00000002681229*(trapENFCal)*(trapENFCal)))/-0.0581112);
    if    (channel==611) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.05027380433695)-(0.00655200480122*(trapENFCal))-(-0.00000002506488*(trapENFCal)*(trapENFCal)))/-0.0651642);
    if    (channel==608) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.08762236936258)-(0.00673752674219*(trapENFCal))-(-0.00000006087619*(trapENFCal)*(trapENFCal)))/-0.0491428);
    if    (channel==609) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.07787163252312)-(0.00677962842063*(trapENFCal))-(-0.00000005729099*(trapENFCal)*(trapENFCal)))/-0.0535326);
    if    (channel==664) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(0.00168443667268)-(0.00545865678715*(trapENFCal))-(-0.00000004582835*(trapENFCal)*(trapENFCal)))/-0.072863);
    if    (channel==665) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(0.01185450788664)-(0.00547756401790*(trapENFCal))-(-0.00000005037805*(trapENFCal)*(trapENFCal)))/-0.0777594);
    if    (channel==624) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.03196232221245)-(0.00669630305457*(trapENFCal))-(-0.00000006856325*(trapENFCal)*(trapENFCal)))/-0.0220731);
    if    (channel==625) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.03172665453744)-(0.00679070274489*(trapENFCal))-(-0.00000009629184*(trapENFCal)*(trapENFCal)))/-0.048417);
    if    (channel==694) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.07371634328114)-(0.00619574151085*(trapENFCal))-(-0.00000002946204*(trapENFCal)*(trapENFCal)))/-0.0296665);
    if    (channel==695) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.07251210743669)-(0.00625933953258*(trapENFCal))-(-0.00000004110540*(trapENFCal)*(trapENFCal)))/-0.0530002);
    if    (channel==614) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.04306365346654)-(0.00548092183053*(trapENFCal))-(-0.00000003122997*(trapENFCal)*(trapENFCal)))/-0.0348515);
    if    (channel==615) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.00943573405492)-(0.00543688230843*(trapENFCal))-(-0.00000000242063*(trapENFCal)*(trapENFCal)))/-0.0168583);
    if    (channel==678) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.03266810971048)-(0.00625155169960*(trapENFCal))-(-0.00000002525494*(trapENFCal)*(trapENFCal)))/-0.127189);
    if    (channel==679) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.00815352577182)-(0.00626419919635*(trapENFCal))-(-0.00000003489606*(trapENFCal)*(trapENFCal)))/-0.10936);
    if    (channel==672) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.06694992954694)-(0.00649629295173*(trapENFCal))-(-0.00000005858563*(trapENFCal)*(trapENFCal)))/-0.104179);
    if    (channel==673) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.02420738047400)-(0.00643537603280*(trapENFCal))-(-0.00000001162148*(trapENFCal)*(trapENFCal)))/-0.115225);
    if    (channel==632) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.01034719451947)-(0.00616626919649*(trapENFCal))-(-0.00000000601213*(trapENFCal)*(trapENFCal)))/-0.0523791);
    if    (channel==633) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.01583930828042)-(0.00623585291784*(trapENFCal))-(-0.00000003613045*(trapENFCal)*(trapENFCal)))/-0.065838);
    if    (channel==626) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.07853390754111)-(0.00634137282673*(trapENFCal))-(-0.00000003428436*(trapENFCal)*(trapENFCal)))/-0.0730604);
    if    (channel==627) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.06225130747970)-(0.00647374783928*(trapENFCal))-(-0.00000003735658*(trapENFCal)*(trapENFCal)))/-0.0660374);
    if    (channel==690) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.04078165666878)-(0.00619888558137*(trapENFCal))-(-0.00000002544554*(trapENFCal)*(trapENFCal)))/-0.0438099);
    if    (channel==691) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.05340367688251)-(0.00638822824024*(trapENFCal))-(-0.00000003475414*(trapENFCal)*(trapENFCal)))/-0.0509051);
    if    (channel==600) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.07582324092566)-(0.00603959595659*(trapENFCal))-(-0.00000004553222*(trapENFCal)*(trapENFCal)))/-0.0666385);
    if    (channel==601) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.02325229400866)-(0.00600601204756*(trapENFCal))-(-0.00000003416987*(trapENFCal)*(trapENFCal)))/-0.0431727);
    if    (channel==598) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.04159923205737)-(0.00623997890391*(trapENFCal))-(-0.00000003500470*(trapENFCal)*(trapENFCal)))/-0.0659015);
    if    (channel==599) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(0.03348663613225)-(0.00610943051544*(trapENFCal))-(0.00000000849683*(trapENFCal)*(trapENFCal)))/-0.059581);
    if    (channel==594) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.12943193566487)-(0.00628691811907*(trapENFCal))-(-0.00000007468196*(trapENFCal)*(trapENFCal)))/-0.0482716);
    if    (channel==592) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(0.02180669907036)-(0.00511006227965*(trapENFCal))-(-0.00000001470141*(trapENFCal)*(trapENFCal)))/-0.241891);
    if    (channel==593) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.03063844341763)-(0.00525325302212*(trapENFCal))-(-0.00000006921246*(trapENFCal)*(trapENFCal)))/-0.189573);
  }

  else if(dsNumber == 4) {
    if    (channel==1204) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.04492676133258)-(0.00440276781199*(trapENFCal))-(-0.00000003937151*(trapENFCal)*(trapENFCal)))/-0.0219138);
    if    (channel==1205) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.03214880397544)-(0.00442416880221*(trapENFCal))-(-0.00000003438751*(trapENFCal)*(trapENFCal)))/-0.0308027);
    if    (channel==1174) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.05627726461000)-(0.00564090989246*(trapENFCal))-(-0.00000002930835*(trapENFCal)*(trapENFCal)))/-0.0243243);
    if    (channel==1144) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.06332310559740)-(0.00599764085195*(trapENFCal))-(-0.00000003918864*(trapENFCal)*(trapENFCal)))/-0.0446303);
    if    (channel==1145) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.06266688638121)-(0.00604620690710*(trapENFCal))-(-0.00000005354455*(trapENFCal)*(trapENFCal)))/-0.068454);
    if    (channel==1106) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.04529584551927)-(0.00651802637111*(trapENFCal))-(-0.00000005349529*(trapENFCal)*(trapENFCal)))/-0.0322532);
    if    (channel==1107) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.02315328560714)-(0.00666508533271*(trapENFCal))-(-0.00000005222197*(trapENFCal)*(trapENFCal)))/-0.035906);
    if    (channel==1176) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.06393963024283)-(0.00455079156776*(trapENFCal))-(-0.00000005155169*(trapENFCal)*(trapENFCal)))/-0.0607681);
    if    (channel==1177) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.03132838326082)-(0.00454712193412*(trapENFCal))-(-0.00000005810435*(trapENFCal)*(trapENFCal)))/-0.0598209);
    if    (channel==1172) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.04708588875779)-(0.00469099720240*(trapENFCal))-(-0.00000002590479*(trapENFCal)*(trapENFCal)))/-0.100638);
    if    (channel==1173) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.02211835105220)-(0.00471355069226*(trapENFCal))-(-0.00000002806763*(trapENFCal)*(trapENFCal)))/-0.108748);
    if    (channel==1170) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.03926728269556)-(0.00476128290556*(trapENFCal))-(-0.00000001885519*(trapENFCal)*(trapENFCal)))/-0.0155709);
    if    (channel==1171) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.02591239821003)-(0.00487749086294*(trapENFCal))-(-0.00000001925330*(trapENFCal)*(trapENFCal)))/-0.0437917);
    if    (channel==1330) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.06581667873247)-(0.00564316330196*(trapENFCal))-(-0.00000002618512*(trapENFCal)*(trapENFCal)))/-0.021023);
    if    (channel==1331) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.03341117475806)-(0.00572173918090*(trapENFCal))-(-0.00000001203680*(trapENFCal)*(trapENFCal)))/-0.0282779);
    if    (channel==1136) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.06183800184960)-(0.00444858206971*(trapENFCal))-(-0.00000003570746*(trapENFCal)*(trapENFCal)))/-0.0170052);
    if    (channel==1137) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.04514103255668)-(0.00445998785948*(trapENFCal))-(-0.00000004218729*(trapENFCal)*(trapENFCal)))/-0.037604);
    if    (channel==1332) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.05020837432235)-(0.00646895691849*(trapENFCal))-(-0.00000005016942*(trapENFCal)*(trapENFCal)))/-0.111793);
    if    (channel==1333) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.06522379230084)-(0.00655720014415*(trapENFCal))-(-0.00000005509542*(trapENFCal)*(trapENFCal)))/-0.129923);
    if    (channel==1296) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.02905902234780)-(0.00635690240881*(trapENFCal))-(-0.00000005804185*(trapENFCal)*(trapENFCal)))/-1.21662);
    if    (channel==1297) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.10080734537872)-(0.00659471659076*(trapENFCal))-(-0.00000011353877*(trapENFCal)*(trapENFCal)))/-1.3617);
    if    (channel==1298) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.02602283663240)-(0.00392855907799*(trapENFCal))-(-0.00000002441196*(trapENFCal)*(trapENFCal)))/-0.0237908);
    if    (channel==1299) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.03819461764877)-(0.00406795064026*(trapENFCal))-(-0.00000004724361*(trapENFCal)*(trapENFCal)))/-0.0423676);
    if    (channel==1236) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.06964458207570)-(0.00610095710446*(trapENFCal))-(-0.00000004577157*(trapENFCal)*(trapENFCal)))/-0.0533531);
    if    (channel==1237) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.02098993091207)-(0.00611008983295*(trapENFCal))-(-0.00000004433716*(trapENFCal)*(trapENFCal)))/-0.0703124);
    if    (channel==1232) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.07091244235699)-(0.00492158632910*(trapENFCal))-(-0.00000003514492*(trapENFCal)*(trapENFCal)))/-0.024136);
    if    (channel==1233) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(-0.04787974448220)-(0.00494105584215*(trapENFCal))-(-0.00000003379360*(trapENFCal)*(trapENFCal)))/-0.0379984);
  }
  //else cout << "GetAvsE(): unknown dataset number DS" << dsNumber << endl;
  return 0.0;
}

void LoadLNFillTimes1(vector<double>& lnFillTimes1, int dsNumber)
{
  // we don't really need to make DS-specific lists, but look-up
  // time is shorter if the lists are shorter.

  if(dsNumber == 0) {
    lnFillTimes1.push_back(1435870160);
    lnFillTimes1.push_back(1436020533);
    lnFillTimes1.push_back(1436168622);
    lnFillTimes1.push_back(1436316362);
    lnFillTimes1.push_back(1436463827);
    lnFillTimes1.push_back(1436614832);
    lnFillTimes1.push_back(1436759402);
    lnFillTimes1.push_back(1436908966);
    lnFillTimes1.push_back(1437058884);
    lnFillTimes1.push_back(1437060687);
    lnFillTimes1.push_back(1437204895);
    lnFillTimes1.push_back(1437350621);
    lnFillTimes1.push_back(1437638022);
    lnFillTimes1.push_back(1437781652);
    lnFillTimes1.push_back(1437926383);
    lnFillTimes1.push_back(1438070819);
    lnFillTimes1.push_back(1438121812);
    lnFillTimes1.push_back(1438212611);
    lnFillTimes1.push_back(1438352837);
    lnFillTimes1.push_back(1438496004);
    lnFillTimes1.push_back(1438639835);
    lnFillTimes1.push_back(1438782282);
    lnFillTimes1.push_back(1438932169);
    lnFillTimes1.push_back(1439080850);
    lnFillTimes1.push_back(1439226624);
    lnFillTimes1.push_back(1439371591);
    lnFillTimes1.push_back(1439514291);
    lnFillTimes1.push_back(1439657247);
    lnFillTimes1.push_back(1439801649);
    lnFillTimes1.push_back(1439944622);
    lnFillTimes1.push_back(1440080091);
    lnFillTimes1.push_back(1440192295);
    lnFillTimes1.push_back(1440335863);
    lnFillTimes1.push_back(1440477294);
    lnFillTimes1.push_back(1440618411);
    lnFillTimes1.push_back(1440759782);
    lnFillTimes1.push_back(1440902658);
    lnFillTimes1.push_back(1441046730);
    lnFillTimes1.push_back(1441187019);
    lnFillTimes1.push_back(1441325878);
    lnFillTimes1.push_back(1441462000);
    lnFillTimes1.push_back(1441600116);
    lnFillTimes1.push_back(1441741779);
    lnFillTimes1.push_back(1441883940);
    lnFillTimes1.push_back(1442027368);
    lnFillTimes1.push_back(1442169713);
    lnFillTimes1.push_back(1442312599);
    lnFillTimes1.push_back(1442453920);
    lnFillTimes1.push_back(1442595578);
    lnFillTimes1.push_back(1442737259);
    lnFillTimes1.push_back(1442879000);
    lnFillTimes1.push_back(1443021647);
  }

  if(dsNumber == 1) {
    lnFillTimes1.push_back(1452519860);
    lnFillTimes1.push_back(1452655867);
    lnFillTimes1.push_back(1452791415);
    lnFillTimes1.push_back(1453541032);
    lnFillTimes1.push_back(1453670314);
    lnFillTimes1.push_back(1453800137);
    lnFillTimes1.push_back(1453929779);
    lnFillTimes1.push_back(1453937178);
    lnFillTimes1.push_back(1453998125);
    lnFillTimes1.push_back(1454000591);
    lnFillTimes1.push_back(1454002456);
    lnFillTimes1.push_back(1454014971);
    lnFillTimes1.push_back(1454107981);
    lnFillTimes1.push_back(1454219392);
    lnFillTimes1.push_back(1454332307);
    lnFillTimes1.push_back(1454447212);
    lnFillTimes1.push_back(1454559953);
    lnFillTimes1.push_back(1454679496);
    lnFillTimes1.push_back(1454769079);
    lnFillTimes1.push_back(1454882301);
    lnFillTimes1.push_back(1454946492);
    lnFillTimes1.push_back(1454951907);
    lnFillTimes1.push_back(1454954012);
    lnFillTimes1.push_back(1454955958);
    lnFillTimes1.push_back(1455225161);
    lnFillTimes1.push_back(1455228289);
    lnFillTimes1.push_back(1455440690);
    lnFillTimes1.push_back(1455568585);
    lnFillTimes1.push_back(1455696789);
    lnFillTimes1.push_back(1455822149);
    lnFillTimes1.push_back(1455952573);
    lnFillTimes1.push_back(1456082166);
    lnFillTimes1.push_back(1456206057);
    lnFillTimes1.push_back(1456333236);
    lnFillTimes1.push_back(1456460237);
    lnFillTimes1.push_back(1456588495);
    lnFillTimes1.push_back(1456717776);
    lnFillTimes1.push_back(1456846882);
    lnFillTimes1.push_back(1456943983);
    lnFillTimes1.push_back(1456981934);
    lnFillTimes1.push_back(1457110918);
    lnFillTimes1.push_back(1457238098);
    lnFillTimes1.push_back(1457365179);
    lnFillTimes1.push_back(1457491997);
    lnFillTimes1.push_back(1457619662);
    lnFillTimes1.push_back(1457747884);
    lnFillTimes1.push_back(1457874684);
    lnFillTimes1.push_back(1458001143);
    lnFillTimes1.push_back(1458130495);
    lnFillTimes1.push_back(1458259402);
    lnFillTimes1.push_back(1458387930);
    lnFillTimes1.push_back(1458515657);
    lnFillTimes1.push_back(1458639685);
    lnFillTimes1.push_back(1458767518);
    lnFillTimes1.push_back(1458896260);
    lnFillTimes1.push_back(1459023862);
    lnFillTimes1.push_back(1459150863);
    lnFillTimes1.push_back(1459275989);
    lnFillTimes1.push_back(1459402548);
    lnFillTimes1.push_back(1459529663);
    lnFillTimes1.push_back(1459654930);
    lnFillTimes1.push_back(1459785212);
    lnFillTimes1.push_back(1459912507);
    lnFillTimes1.push_back(1460042708);
    lnFillTimes1.push_back(1460169429);
    lnFillTimes1.push_back(1460297779);
    lnFillTimes1.push_back(1460426527);
    lnFillTimes1.push_back(1460552742);
    lnFillTimes1.push_back(1460678641);
    lnFillTimes1.push_back(1460808916);
    lnFillTimes1.push_back(1460939150);
    lnFillTimes1.push_back(1461064619);
    lnFillTimes1.push_back(1461189868);
    lnFillTimes1.push_back(1461318559);
    lnFillTimes1.push_back(1461443967);
    lnFillTimes1.push_back(1461569658);
    lnFillTimes1.push_back(1461695992);
    lnFillTimes1.push_back(1461824752);
    lnFillTimes1.push_back(1461952713);
    lnFillTimes1.push_back(1462258909);
    lnFillTimes1.push_back(1462314666);
    lnFillTimes1.push_back(1462371538);
    lnFillTimes1.push_back(1462426019);
    lnFillTimes1.push_back(1462462316);
    lnFillTimes1.push_back(1462588708);
    lnFillTimes1.push_back(1462715136);
    lnFillTimes1.push_back(1462840823);
    lnFillTimes1.push_back(1462966744);
    lnFillTimes1.push_back(1463093301);
    lnFillTimes1.push_back(1463220332);
    lnFillTimes1.push_back(1463346237);
    lnFillTimes1.push_back(1463472892);
    lnFillTimes1.push_back(1463599111);
    lnFillTimes1.push_back(1463694217);
  }

   if(dsNumber == 3) {
   //M1
    lnFillTimes1.push_back(1472129474);
    lnFillTimes1.push_back(1472268007);
    lnFillTimes1.push_back(1472404916);
    lnFillTimes1.push_back(1472545253);
    lnFillTimes1.push_back(1472683456);
    lnFillTimes1.push_back(1472822134);
    lnFillTimes1.push_back(1472961379);
    lnFillTimes1.push_back(1473096880);
    lnFillTimes1.push_back(1473236517);
    lnFillTimes1.push_back(1473369324);
    lnFillTimes1.push_back(1473510316);
    lnFillTimes1.push_back(1473635684);
    lnFillTimes1.push_back(1473771858);
    lnFillTimes1.push_back(1473909355);
    lnFillTimes1.push_back(1474045268);
    lnFillTimes1.push_back(1474182801);
    lnFillTimes1.push_back(1474319605);
    lnFillTimes1.push_back(1474457662);
    lnFillTimes1.push_back(1474569273);
    lnFillTimes1.push_back(1474705001);
    lnFillTimes1.push_back(1474842466);
    lnFillTimes1.push_back(1474977413);
   }
}
   void LoadLNFillTimes2(vector<double>& lnFillTimes2, int dsNumber)
   {
    if(dsNumber == 4) {
    //M2
    lnFillTimes2.push_back(1472057495);
    lnFillTimes2.push_back(1472169344);
    lnFillTimes2.push_back(1472283139);
    lnFillTimes2.push_back(1472392397);
    lnFillTimes2.push_back(1472508671);
    lnFillTimes2.push_back(1472626749);
    lnFillTimes2.push_back(1472749055);
    lnFillTimes2.push_back(1472848201);
    lnFillTimes2.push_back(1472959270);
    lnFillTimes2.push_back(1473069251);
    lnFillTimes2.push_back(1473180721);
    lnFillTimes2.push_back(1473291137);
    lnFillTimes2.push_back(1473369625);
    lnFillTimes2.push_back(1473485389);
    lnFillTimes2.push_back(1473594355);
    lnFillTimes2.push_back(1473702993);
    lnFillTimes2.push_back(1473815615);
    lnFillTimes2.push_back(1473926807);
    lnFillTimes2.push_back(1474040572);
    lnFillTimes2.push_back(1474150140);
    lnFillTimes2.push_back(1474260642);
    lnFillTimes2.push_back(1474370902);
    lnFillTimes2.push_back(1474482449);
    lnFillTimes2.push_back(1474569695);
    lnFillTimes2.push_back(1474679736);
    lnFillTimes2.push_back(1474789636);
    lnFillTimes2.push_back(1474901605);
    lnFillTimes2.push_back(1475011276);
  }
}

double GetDCRraw(int channel, double nlcblrwfSlope, double trapMax, int dsNumber)
{
  if(dsNumber == 0) {
    if(channel == 576) return nlcblrwfSlope-( 3.673385E-05+trapMax*-1.279E-05);
    if(channel == 577) return nlcblrwfSlope-( 3.122161E-05+trapMax*-1.278E-05);
    if(channel == 592) return nlcblrwfSlope-( 1.371969E-04+trapMax*-1.287E-05);
    if(channel == 593) return nlcblrwfSlope-( 4.661555E-05+trapMax*-1.285E-05);
    if(channel == 594) return nlcblrwfSlope-( 7.684589E-05+trapMax*-1.294E-05);
    if(channel == 595) return nlcblrwfSlope-( 5.906154E-05+trapMax*-1.294E-05);
    if(channel == 598) return nlcblrwfSlope-( 1.842168E-04+trapMax*-1.315E-05);
    if(channel == 599) return nlcblrwfSlope-( 7.434566E-05+trapMax*-1.315E-05);
    if(channel == 600) return nlcblrwfSlope-( 5.645349E-05+trapMax*-1.264E-05);
    if(channel == 601) return nlcblrwfSlope-( 1.644013E-05+trapMax*-1.260E-05);
    if(channel == 608) return nlcblrwfSlope-( 9.091576E-05+trapMax*-1.298E-05);
    if(channel == 609) return nlcblrwfSlope-( 1.319922E-05+trapMax*-1.296E-05);
    if(channel == 610) return nlcblrwfSlope-( 3.708218E-05+trapMax*-1.331E-05);
    if(channel == 611) return nlcblrwfSlope-( 1.691204E-05+trapMax*-1.331E-05);
    if(channel == 614) return nlcblrwfSlope-( 7.716150E-05+trapMax*-1.272E-05);
    if(channel == 615) return nlcblrwfSlope-( 1.479029E-05+trapMax*-1.271E-05);
    if(channel == 624) return nlcblrwfSlope-( 3.891078E-05+trapMax*-1.262E-05);
    if(channel == 625) return nlcblrwfSlope-( 4.344515E-06+trapMax*-1.262E-05);
    if(channel == 626) return nlcblrwfSlope-( 2.578407E-05+trapMax*-1.298E-05);
    if(channel == 627) return nlcblrwfSlope-( 1.913080E-05+trapMax*-1.296E-05);
    if(channel == 628) return nlcblrwfSlope-( 3.266857E-05+trapMax*-1.297E-05);
    if(channel == 629) return nlcblrwfSlope-( 2.171356E-05+trapMax*-1.295E-05);
    if(channel == 640) return nlcblrwfSlope-( 7.818519E-05+trapMax*-1.296E-05);
    if(channel == 641) return nlcblrwfSlope-( 3.668705E-05+trapMax*-1.296E-05);
    if(channel == 642) return nlcblrwfSlope-( 5.146852E-05+trapMax*-1.373E-05);
    if(channel == 643) return nlcblrwfSlope-( 2.370897E-05+trapMax*-1.375E-05);
    if(channel == 644) return nlcblrwfSlope-( 6.066896E-05+trapMax*-1.308E-05);
    if(channel == 645) return nlcblrwfSlope-(-1.188733E-05+trapMax*-1.305E-05);
    if(channel == 646) return nlcblrwfSlope-( 7.972405E-05+trapMax*-1.291E-05);
    if(channel == 647) return nlcblrwfSlope-( 2.540927E-05+trapMax*-1.286E-05);
    if(channel == 656) return nlcblrwfSlope-( 1.262340E-05+trapMax*-1.337E-05);
    if(channel == 657) return nlcblrwfSlope-( 1.086301E-05+trapMax*-1.342E-05);
    if(channel == 662) return nlcblrwfSlope-( 1.485632E-04+trapMax*-1.460E-05);
    if(channel == 663) return nlcblrwfSlope-( 2.742351E-06+trapMax*-1.457E-05);
    if(channel == 664) return nlcblrwfSlope-( 5.110392E-05+trapMax*-1.321E-05);
    if(channel == 665) return nlcblrwfSlope-( 6.499560E-06+trapMax*-1.319E-05);
    if(channel == 674) return nlcblrwfSlope-( 9.983398E-05+trapMax*-1.310E-05);
    if(channel == 675) return nlcblrwfSlope-( 3.757473E-05+trapMax*-1.309E-05);
    if(channel == 688) return nlcblrwfSlope-( 2.461617E-05+trapMax*-1.287E-05);
    if(channel == 689) return nlcblrwfSlope-( 3.018505E-05+trapMax*-1.287E-05);
    if(channel == 690) return nlcblrwfSlope-( 9.944427E-05+trapMax*-1.300E-05);
    if(channel == 691) return nlcblrwfSlope-(-4.737995E-06+trapMax*-1.300E-05);
    if(channel == 692) return nlcblrwfSlope-( 8.146529E-05+trapMax*-1.354E-05);
    if(channel == 693) return nlcblrwfSlope-( 1.839544E-05+trapMax*-1.356E-05);
    if(channel == 696) return nlcblrwfSlope-( 2.093490E-05+trapMax*-1.308E-05);
    if(channel == 697) return nlcblrwfSlope-( 8.938457E-06+trapMax*-1.310E-05);
    if(trapMax == 0) return 0;
    return nlcblrwfSlope/trapMax;
  }

  if(dsNumber == 1) {//updated 17 Aug 2016 using runs 11507-11592
    if(channel == 578) return nlcblrwfSlope-(4.620083E-05+trapMax*-1.306E-05) ;
    if(channel == 579) return nlcblrwfSlope-(1.163801E-05+trapMax*-1.304E-05) ;
    if(channel == 580) return nlcblrwfSlope-(2.604215E-05+trapMax*-1.291E-05) ;
    if(channel == 581) return nlcblrwfSlope-(2.994274E-05+trapMax*-1.290E-05) ;
    if(channel == 582) return nlcblrwfSlope-(5.175215E-05+trapMax*-1.304E-05) ;
    if(channel == 583) return nlcblrwfSlope-(1.902729E-05+trapMax*-1.302E-05) ;
    if(channel == 592) return nlcblrwfSlope-(5.145000E-05+trapMax*-1.297E-05) ;
    if(channel == 593) return nlcblrwfSlope-(4.298575E-05+trapMax*-1.297E-05) ;
    if(channel == 594) return nlcblrwfSlope-(7.298870E-05+trapMax*-1.353E-05) ;
    if(channel == 595) return nlcblrwfSlope-(6.178016E-06+trapMax*-1.351E-05) ;
    if(channel == 598) return nlcblrwfSlope-(3.497582E-05+trapMax*-1.298E-05) ;
    if(channel == 599) return nlcblrwfSlope-(9.553672E-06+trapMax*-1.298E-05) ;
    if(channel == 600) return nlcblrwfSlope-(5.649999E-05+trapMax*-1.287E-05) ;
    if(channel == 601) return nlcblrwfSlope-(5.482757E-05+trapMax*-1.288E-05) ;
    if(channel == 608) return nlcblrwfSlope-(4.586047E-05+trapMax*-1.284E-05) ;
    if(channel == 609) return nlcblrwfSlope-(1.016108E-05+trapMax*-1.282E-05) ;
    if(channel == 610) return nlcblrwfSlope-(2.209244E-05+trapMax*-1.298E-05) ;
    if(channel == 611) return nlcblrwfSlope-(3.159581E-05+trapMax*-1.299E-05) ;
    if(channel == 616) return nlcblrwfSlope-(7.722590E-05+trapMax*-1.316E-05) ;
    if(channel == 617) return nlcblrwfSlope-(6.568340E-05+trapMax*-1.309E-05) ;
    if(channel == 626) return nlcblrwfSlope-(1.786671E-05+trapMax*-1.266E-05) ;
    if(channel == 627) return nlcblrwfSlope-(8.580088E-06+trapMax*-1.264E-05) ;
    if(channel == 632) return nlcblrwfSlope-(9.481784E-05+trapMax*-1.307E-05) ;
    if(channel == 633) return nlcblrwfSlope-(6.909768E-06+trapMax*-1.304E-05) ;
    if(channel == 640) return nlcblrwfSlope-(7.977354E-05+trapMax*-1.287E-05) ;
    if(channel == 641) return nlcblrwfSlope-(3.680967E-05+trapMax*-1.287E-05) ;
    if(channel == 648) return nlcblrwfSlope-(9.035880E-05+trapMax*-1.312E-05) ;
    if(channel == 649) return nlcblrwfSlope-(6.199340E-06+trapMax*-1.308E-05) ;
    if(channel == 664) return nlcblrwfSlope-(5.873798E-05+trapMax*-1.321E-05) ;
    if(channel == 665) return nlcblrwfSlope-(2.416351E-05+trapMax*-1.320E-05) ;
    if(channel == 672) return nlcblrwfSlope-(2.815562E-05+trapMax*-1.328E-05) ;
    if(channel == 673) return nlcblrwfSlope-(4.809700E-05+trapMax*-1.334E-05) ;
    if(channel == 690) return nlcblrwfSlope-(4.148262E-05+trapMax*-1.308E-05) ;
    if(channel == 691) return nlcblrwfSlope-(1.779754E-05+trapMax*-1.310E-05) ;
    if(channel == 692) return nlcblrwfSlope-(4.737009E-05+trapMax*-1.309E-05) ;
    if(channel == 693) return nlcblrwfSlope-(2.395867E-05+trapMax*-1.312E-05) ;

    if(trapMax == 0) return 0;
    return nlcblrwfSlope/trapMax;
  }

  if(dsNumber == 3) {
    //M1 updated 14 Nov 2016 using runs 17183-17302, with DS 3 AvsE cut
    if (channel == 578) return (nlcblrwfSlope-(5.561101E-05+trapMax*-3.209E-05));
    if (channel == 579) return (nlcblrwfSlope-(2.685407E-05+trapMax*-9.546E-06));
    if (channel == 580) return (nlcblrwfSlope-(7.351231E-05+trapMax*-3.186E-05));
    if (channel == 581) return (nlcblrwfSlope-(3.017863E-05+trapMax*-9.503E-06));
    if (channel == 582) return (nlcblrwfSlope-(9.632288E-05+trapMax*-3.431E-05));
    if (channel == 592) return (nlcblrwfSlope-(1.567474E-04+trapMax*-3.107E-05));
    if (channel == 593) return (nlcblrwfSlope-(7.869555E-05+trapMax*-9.440E-06));
    if (channel == 594) return (nlcblrwfSlope-(1.822680E-04+trapMax*-3.382E-05));
    if (channel == 598) return (nlcblrwfSlope-(4.848408E-05+trapMax*-3.198E-05));
    if (channel == 599) return (nlcblrwfSlope-(3.762619E-05+trapMax*-9.507E-06));
    if (channel == 600) return (nlcblrwfSlope-(6.206060E-05+trapMax*-3.107E-05));
    if (channel == 601) return (nlcblrwfSlope-(6.057838E-05+trapMax*-9.191E-06));
    if (channel == 608) return (nlcblrwfSlope-(6.751434E-05+trapMax*-3.200E-05));
    if (channel == 609) return (nlcblrwfSlope-(2.672019E-05+trapMax*-9.396E-06));
    if (channel == 610) return (nlcblrwfSlope-(3.915655E-05+trapMax*-3.179E-05));
    if (channel == 611) return (nlcblrwfSlope-(4.076262E-05+trapMax*-9.499E-06));
    if (channel == 614) return (nlcblrwfSlope-(6.801647E-05+trapMax*-1.276E-05));
    if (channel == 615) return (nlcblrwfSlope-(3.130201E-05+trapMax*-1.278E-05));
    if (channel == 624) return (nlcblrwfSlope-(1.028492E-04+trapMax*-1.304E-05));
    if (channel == 625) return (nlcblrwfSlope-(-1.356258E-05+trapMax*-1.302E-05));
    if (channel == 626) return (nlcblrwfSlope-(3.598488E-05+trapMax*-3.202E-05));
    if (channel == 627) return (nlcblrwfSlope-(1.723288E-05+trapMax*-9.479E-06));
    if (channel == 632) return (nlcblrwfSlope-(1.160388E-04+trapMax*-3.170E-05));
    if (channel == 633) return (nlcblrwfSlope-(1.645620E-05+trapMax*-9.467E-06));
    if (channel == 640) return (nlcblrwfSlope-(9.364262E-05+trapMax*-3.211E-05));
    if (channel == 641) return (nlcblrwfSlope-(3.945656E-05+trapMax*-9.469E-06));
    if (channel == 648) return (nlcblrwfSlope-(1.547883E-04+trapMax*-3.272E-05));
    if (channel == 649) return (nlcblrwfSlope-(4.500279E-05+trapMax*-9.744E-06));
    if (channel == 664) return (nlcblrwfSlope-(8.303413E-05+trapMax*-3.208E-05));
    if (channel == 665) return (nlcblrwfSlope-(2.918286E-05+trapMax*-9.616E-06));
    if (channel == 672) return (nlcblrwfSlope-(5.707301E-05+trapMax*-3.276E-05));
    if (channel == 673) return (nlcblrwfSlope-(3.122337E-05+trapMax*-9.688E-06));
    if (channel == 678) return (nlcblrwfSlope-(1.379480E-04+trapMax*-1.461E-05));
    if (channel == 679) return (nlcblrwfSlope-(8.178871E-05+trapMax*-1.462E-05));
    if (channel == 690) return (nlcblrwfSlope-(4.751138E-05+trapMax*-3.223E-05));
    if (channel == 691) return (nlcblrwfSlope-(2.206048E-05+trapMax*-9.579E-06));
    if (channel == 692) return (nlcblrwfSlope-(7.984028E-05+trapMax*-3.270E-05));
    if (channel == 693) return (nlcblrwfSlope-(4.228087E-05+trapMax*-9.851E-06));
    if (channel == 694) return (nlcblrwfSlope-(7.351794E-05+trapMax*-1.288E-05));
    if (channel == 695) return (nlcblrwfSlope-(6.237732E-05+trapMax*-1.288E-05));
	//undepleted channels, use for veto only
	if (channel == 616) return (nlcblrwfSlope-(4.166850E-04+trapMax*-3.260E-05));
	if (channel == 617) return (nlcblrwfSlope-(2.040554E-04+trapMax*-9.619E-06));
	if (channel == 628) return (nlcblrwfSlope-(1.243452E-04+trapMax*-1.319E-05));
	if (channel == 629) return (nlcblrwfSlope-(3.542933E-05+trapMax*-1.315E-05));
	if (channel == 688) return (nlcblrwfSlope-(1.328691E-04+trapMax*-1.268E-05));
	if (channel == 689) return (nlcblrwfSlope-(4.139136E-05+trapMax*-1.267E-05));
	}

   if(dsNumber == 4) {
	//M2 updated 28 Sept 2016 using runs 60001207-60001306, with AvsE cut where available
     if (channel == 1106) return nlcblrwfSlope-(1.573309E-04+trapMax*-1.273E-05);
     if (channel == 1107) return nlcblrwfSlope-(6.297092E-05+trapMax*-1.272E-05);
     if (channel == 1110) return nlcblrwfSlope-(6.987162E-05+trapMax*-1.259E-05);
     if (channel == 1111) return nlcblrwfSlope-(-3.050914E-06+trapMax*-1.257E-05);
     if (channel == 1136) return nlcblrwfSlope-(7.232481E-05+trapMax*-1.265E-05);
     if (channel == 1137) return nlcblrwfSlope-(2.011832E-05+trapMax*-1.263E-05);
     if (channel == 1140) return nlcblrwfSlope-(9.959902E-05+trapMax*-1.249E-05);
     if (channel == 1141) return nlcblrwfSlope-(3.367713E-05+trapMax*-1.248E-05);
     if (channel == 1142) return nlcblrwfSlope-(1.673876E-04+trapMax*-1.206E-05);//no AvsE
     if (channel == 1143) return nlcblrwfSlope-(9.746344E-05+trapMax*-1.208E-05);//no AvsE
     if (channel == 1144) return nlcblrwfSlope-(2.081933E-04+trapMax*-1.268E-05);
     if (channel == 1145) return nlcblrwfSlope-(-2.144171E-05+trapMax*-1.260E-05);
     if (channel == 1170) return nlcblrwfSlope-(1.659349E-04+trapMax*-1.282E-05);
     if (channel == 1171) return nlcblrwfSlope-(1.156556E-05+trapMax*-1.279E-05);
     if (channel == 1172) return nlcblrwfSlope-(1.015104E-04+trapMax*-1.277E-05);
     if (channel == 1173) return nlcblrwfSlope-(6.846852E-05+trapMax*-1.278E-05);
     if (channel == 1174) return nlcblrwfSlope-(9.796817E-05+trapMax*-1.291E-05);
     if (channel == 1175) return nlcblrwfSlope-(1.164940E-04+trapMax*-1.296E-05);
     if (channel == 1176) return nlcblrwfSlope-(1.041141E-04+trapMax*-1.297E-05);
     if (channel == 1177) return nlcblrwfSlope-(3.407831E-06+trapMax*-1.290E-05);
     if (channel == 1204) return nlcblrwfSlope-(8.106601E-05+trapMax*-1.295E-05);
     if (channel == 1205) return nlcblrwfSlope-(6.240128E-05+trapMax*-1.297E-05);
     if (channel == 1208) return nlcblrwfSlope-(4.656507E-05+trapMax*-1.278E-05);
     if (channel == 1209) return nlcblrwfSlope-(1.959101E-05+trapMax*-1.277E-05);
     if (channel == 1232) return nlcblrwfSlope-(1.405523E-04+trapMax*-1.301E-05);
     if (channel == 1233) return nlcblrwfSlope-(-1.258021E-05+trapMax*-1.294E-05);
     if (channel == 1236) return nlcblrwfSlope-(3.628867E-05+trapMax*-1.318E-05);
     if (channel == 1237) return nlcblrwfSlope-(8.980344E-05+trapMax*-1.323E-05);
     if (channel == 1238) return nlcblrwfSlope-(2.571538E-04+trapMax*-1.294E-05);//no AvsE
     if (channel == 1239) return nlcblrwfSlope-(-1.758461E-04+trapMax*-1.269E-05);
     if (channel == 1296) return nlcblrwfSlope-(4.452692E-04+trapMax*-1.350E-05);//noAvsE
     if (channel == 1297) return nlcblrwfSlope-(1.846946E-04+trapMax*-1.354E-05);
     if (channel == 1298) return nlcblrwfSlope-(6.628176E-05+trapMax*-1.277E-05);
     if (channel == 1299) return nlcblrwfSlope-(1.287783E-05+trapMax*-1.276E-05);
     if (channel == 1302) return nlcblrwfSlope-(8.618478E-05+trapMax*-1.277E-05);
     if (channel == 1303) return nlcblrwfSlope-(7.087915E-05+trapMax*-1.281E-05);
     if (channel == 1330) return nlcblrwfSlope-(1.790647E-04+trapMax*-1.277E-05);
     if (channel == 1331) return nlcblrwfSlope-(7.959402E-05+trapMax*-1.275E-05);
     if (channel == 1332) return nlcblrwfSlope-(1.120067E-04+trapMax*-1.306E-05);
     if (channel == 1333) return nlcblrwfSlope-(3.282432E-05+trapMax*-1.304E-05);

    if(trapMax == 0) return 0;
    return nlcblrwfSlope/trapMax;
  }


  else //cout << "GetDCRraw(): unknown dataset number DS" << dsNumber << endl;
  return 0;
}

double GetDCR85(int channel, double nlcblrwfSlope, double trapMax, int dsNumber)
{
  if(dsNumber == 0) {
    if(channel == 674) return nlcblrwfSlope-( 9.983398E-05+trapMax*-1.310E-05)-1.404733E-04;
    if(channel == 675) return nlcblrwfSlope-( 3.757473E-05+trapMax*-1.309E-05)-5.853623E-05;
    if(channel == 688) return nlcblrwfSlope-( 2.461617E-05+trapMax*-1.287E-05)-1.969379E-04;
    if(channel == 689) return nlcblrwfSlope-( 3.018505E-05+trapMax*-1.287E-05)-6.711741E-05;
    if(channel == 690) return nlcblrwfSlope-( 9.944427E-05+trapMax*-1.300E-05)-1.273623E-04;
    if(channel == 691) return nlcblrwfSlope-(-4.737995E-06+trapMax*-1.300E-05)-5.400800E-05;
    if(channel == 692) return nlcblrwfSlope-( 8.146529E-05+trapMax*-1.354E-05)-7.336924E-05;
    if(channel == 693) return nlcblrwfSlope-( 1.839544E-05+trapMax*-1.356E-05)-4.684072E-05;
    if(channel == 696) return nlcblrwfSlope-( 2.093490E-05+trapMax*-1.308E-05)-1.193772E-04;
    if(channel == 697) return nlcblrwfSlope-( 8.938457E-06+trapMax*-1.310E-05)-5.622187E-05;
    if(channel == 576) return nlcblrwfSlope-( 3.673385E-05+trapMax*-1.279E-05)-7.634619E-05;
    if(channel == 577) return nlcblrwfSlope-( 3.122161E-05+trapMax*-1.278E-05)-3.847147E-05;
    if(channel == 592) return nlcblrwfSlope-( 1.371969E-04+trapMax*-1.287E-05)-9.726243E-05;
    if(channel == 593) return nlcblrwfSlope-( 4.661555E-05+trapMax*-1.285E-05)-4.177841E-05;
    if(channel == 594) return nlcblrwfSlope-( 7.684589E-05+trapMax*-1.294E-05)-8.493529E-05;
    if(channel == 595) return nlcblrwfSlope-( 5.906154E-05+trapMax*-1.294E-05)-3.811032E-05;
    if(channel == 598) return nlcblrwfSlope-( 1.842168E-04+trapMax*-1.315E-05)-4.101792E-05;
    if(channel == 599) return nlcblrwfSlope-( 7.434566E-05+trapMax*-1.315E-05)-2.216033E-05;
    if(channel == 600) return nlcblrwfSlope-( 5.645349E-05+trapMax*-1.264E-05)-1.110490E-04;
    if(channel == 601) return nlcblrwfSlope-( 1.644013E-05+trapMax*-1.260E-05)-5.854669E-05;
    if(channel == 608) return nlcblrwfSlope-( 9.091576E-05+trapMax*-1.298E-05)-8.393719E-05;
    if(channel == 609) return nlcblrwfSlope-( 1.319922E-05+trapMax*-1.296E-05)-3.169238E-05;
    if(channel == 610) return nlcblrwfSlope-( 3.708218E-05+trapMax*-1.331E-05)-1.104404E-04;
    if(channel == 611) return nlcblrwfSlope-( 1.691204E-05+trapMax*-1.331E-05)-5.060284E-05;
    if(channel == 614) return nlcblrwfSlope-( 7.716150E-05+trapMax*-1.272E-05)-1.857343E-04;
    if(channel == 615) return nlcblrwfSlope-( 1.479029E-05+trapMax*-1.271E-05)-7.794238E-05;
    if(channel == 624) return nlcblrwfSlope-( 3.891078E-05+trapMax*-1.262E-05)-9.827360E-05;
    if(channel == 625) return nlcblrwfSlope-( 4.344515E-06+trapMax*-1.262E-05)-4.216143E-05;
    if(channel == 626) return nlcblrwfSlope-( 2.578407E-05+trapMax*-1.298E-05)-8.512630E-05;
    if(channel == 627) return nlcblrwfSlope-( 1.913080E-05+trapMax*-1.296E-05)-4.005612E-05;
    if(channel == 628) return nlcblrwfSlope-( 3.266857E-05+trapMax*-1.297E-05)-1.161883E-04;
    if(channel == 629) return nlcblrwfSlope-( 2.171356E-05+trapMax*-1.295E-05)-5.772356E-05;
    if(channel == 640) return nlcblrwfSlope-( 7.818519E-05+trapMax*-1.296E-05)-1.079951E-04;
    if(channel == 641) return nlcblrwfSlope-( 3.668705E-05+trapMax*-1.296E-05)-5.054795E-05;
    if(channel == 642) return nlcblrwfSlope-( 5.146852E-05+trapMax*-1.373E-05)-8.726088E-05;
    if(channel == 643) return nlcblrwfSlope-( 2.370897E-05+trapMax*-1.375E-05)-4.125616E-05;
    if(channel == 644) return nlcblrwfSlope-( 6.066896E-05+trapMax*-1.308E-05)-1.252195E-04;
    if(channel == 645) return nlcblrwfSlope-(-1.188733E-05+trapMax*-1.305E-05)-4.269257E-05;
    if(channel == 646) return nlcblrwfSlope-( 7.972405E-05+trapMax*-1.291E-05)-7.771228E-05;
    if(channel == 647) return nlcblrwfSlope-( 2.540927E-05+trapMax*-1.286E-05)-2.347411E-05;
    if(channel == 656) return nlcblrwfSlope-( 1.262340E-05+trapMax*-1.337E-05)-7.504308E-05;
    if(channel == 657) return nlcblrwfSlope-( 1.086301E-05+trapMax*-1.342E-05)-3.671951E-05;
    if(channel == 662) return nlcblrwfSlope-( 1.485632E-04+trapMax*-1.460E-05)-2.332299E-04;
    if(channel == 663) return nlcblrwfSlope-( 2.742351E-06+trapMax*-1.457E-05)-4.664937E-05;
    if(channel == 664) return nlcblrwfSlope-( 5.110392E-05+trapMax*-1.321E-05)-7.504187E-05;
    if(channel == 665) return nlcblrwfSlope-( 6.499560E-06+trapMax*-1.319E-05)-3.670170E-05;

    if(trapMax == 0) return 0;
    return nlcblrwfSlope/trapMax;
  }


  if(dsNumber == 1) {
    //params from 11507-11592
    if(channel == 672) return nlcblrwfSlope-(2.833847E-05+trapMax*-1.328E-05)-8.966599E-05;
    if(channel == 673) return nlcblrwfSlope-(4.734834E-05+trapMax*-1.334E-05)-4.815561E-05;
    if(channel == 690) return nlcblrwfSlope-(3.690947E-05+trapMax*-1.308E-05)-6.733372E-05;
    if(channel == 691) return nlcblrwfSlope-(1.605341E-05+trapMax*-1.310E-05)-3.688012E-05;
    if(channel == 692) return nlcblrwfSlope-(5.681803E-05+trapMax*-1.309E-05)-8.342106E-05;
    if(channel == 693) return nlcblrwfSlope-(2.453714E-05+trapMax*-1.312E-05)-4.434964E-05;
    if(channel == 578) return nlcblrwfSlope-(4.417989E-05+trapMax*-1.306E-05)-1.178722E-04;
    if(channel == 579) return nlcblrwfSlope-(1.106159E-05+trapMax*-1.304E-05)-4.880046E-05;
    if(channel == 580) return nlcblrwfSlope-(2.577260E-05+trapMax*-1.291E-05)-2.066446E-04;
    if(channel == 581) return nlcblrwfSlope-(3.041374E-05+trapMax*-1.290E-05)-6.389357E-05;
    if(channel == 582) return nlcblrwfSlope-(5.396620E-05+trapMax*-1.304E-05)-8.078366E-05;
    if(channel == 583) return nlcblrwfSlope-(1.977025E-05+trapMax*-1.302E-05)-4.854558E-05; //average for all LG channels
    if(channel == 592) return nlcblrwfSlope-(6.803550E-05+trapMax*-1.297E-05)-1.025187E-04;
    if(channel == 593) return nlcblrwfSlope-(4.309706E-05+trapMax*-1.296E-05)-5.053970E-05;
    if(channel == 594) return nlcblrwfSlope-(4.643130E-05+trapMax*-1.352E-05)-8.236275E-05;
    if(channel == 595) return nlcblrwfSlope-(1.550955E-06+trapMax*-1.351E-05)-4.854558E-05; //average for all LG channels
    if(channel == 598) return nlcblrwfSlope-(2.897476E-05+trapMax*-1.298E-05)-1.014487E-04;
    if(channel == 599) return nlcblrwfSlope-(9.612110E-06+trapMax*-1.298E-05)-4.279017E-05;
    if(channel == 600) return nlcblrwfSlope-(5.842046E-05+trapMax*-1.287E-05)-6.607836E-05;
    if(channel == 601) return nlcblrwfSlope-(5.372889E-05+trapMax*-1.288E-05)-3.617333E-05;
    if(channel == 608) return nlcblrwfSlope-(4.693324E-05+trapMax*-1.284E-05)-1.825068E-04;
    if(channel == 609) return nlcblrwfSlope-(1.055481E-05+trapMax*-1.282E-05)-6.448547E-05;
    if(channel == 610) return nlcblrwfSlope-(1.216405E-05+trapMax*-1.298E-05)-9.937699E-05;
    if(channel == 611) return nlcblrwfSlope-(2.884038E-05+trapMax*-1.299E-05)-4.632804E-05;
    if(channel == 616) return nlcblrwfSlope-(8.295462E-05+trapMax*-1.316E-05)-8.402205E-05;
    if(channel == 617) return nlcblrwfSlope-(7.294466E-05+trapMax*-1.310E-05)-1.073512E-04;
    if(channel == 626) return nlcblrwfSlope-(1.860959E-06+trapMax*-1.265E-05)-8.080291E-05;
    if(channel == 627) return nlcblrwfSlope-(7.000642E-06+trapMax*-1.264E-05)-4.500772E-05;
    if(channel == 632) return nlcblrwfSlope-(9.473908E-05+trapMax*-1.307E-05)-1.300685E-04;
    if(channel == 633) return nlcblrwfSlope-(6.874599E-06+trapMax*-1.304E-05)-6.070530E-05;
    if(channel == 640) return nlcblrwfSlope-(8.088126E-05+trapMax*-1.287E-05)-7.355431E-05;
    if(channel == 641) return nlcblrwfSlope-(3.667507E-05+trapMax*-1.287E-05)-4.635975E-05;
    if(channel == 648) return nlcblrwfSlope-(8.317390E-05+trapMax*-1.311E-05)-6.962783E-05;
    if(channel == 649) return nlcblrwfSlope-(4.360668E-06+trapMax*-1.308E-05)-5.785846E-05;
    if(channel == 664) return nlcblrwfSlope-(5.916704E-05+trapMax*-1.321E-05)-1.088147E-04;
    if(channel == 665) return nlcblrwfSlope-(2.373893E-05+trapMax*-1.320E-05)-4.178785E-05;


    if(trapMax == 0) return 0;
    return nlcblrwfSlope/trapMax;
  }

  else //cout << "GetDCR85(): unknown dataset number DS" << dsNumber << endl;
  return 0;
}

double GetDCR90(int channel, double nlcblrwfSlope, double trapMax, int dsNumber)
{
  if(dsNumber == 0) {
    if(channel == 576) return nlcblrwfSlope-( 3.673385E-05+trapMax*-1.279E-05)-1.019405E-04;
    if(channel == 577) return nlcblrwfSlope-( 3.122161E-05+trapMax*-1.278E-05)-4.913911E-05;
    if(channel == 592) return nlcblrwfSlope-( 1.371969E-04+trapMax*-1.287E-05)-1.279785E-04;
    if(channel == 593) return nlcblrwfSlope-( 4.661555E-05+trapMax*-1.285E-05)-5.424818E-05;
    if(channel == 594) return nlcblrwfSlope-( 7.684589E-05+trapMax*-1.294E-05)-1.092631E-04;
    if(channel == 595) return nlcblrwfSlope-( 5.906154E-05+trapMax*-1.294E-05)-4.862555E-05;
    if(channel == 598) return nlcblrwfSlope-( 1.842168E-04+trapMax*-1.315E-05)-7.087017E-05;
    if(channel == 599) return nlcblrwfSlope-( 7.434566E-05+trapMax*-1.315E-05)-3.625546E-05;
    if(channel == 600) return nlcblrwfSlope-( 5.645349E-05+trapMax*-1.264E-05)-1.345640E-04;
    if(channel == 601) return nlcblrwfSlope-( 1.644013E-05+trapMax*-1.260E-05)-7.199881E-05;
    if(channel == 608) return nlcblrwfSlope-( 9.091576E-05+trapMax*-1.298E-05)-1.075847E-04;
    if(channel == 609) return nlcblrwfSlope-( 1.319922E-05+trapMax*-1.296E-05)-4.179749E-05;
    if(channel == 610) return nlcblrwfSlope-( 3.708218E-05+trapMax*-1.331E-05)-1.458138E-04;
    if(channel == 611) return nlcblrwfSlope-( 1.691204E-05+trapMax*-1.331E-05)-6.535173E-05;
    if(channel == 614) return nlcblrwfSlope-( 7.716150E-05+trapMax*-1.272E-05)-2.365925E-04;
    if(channel == 615) return nlcblrwfSlope-( 1.479029E-05+trapMax*-1.271E-05)-9.628056E-05;
    if(channel == 624) return nlcblrwfSlope-( 3.891078E-05+trapMax*-1.262E-05)-1.227120E-04;
    if(channel == 625) return nlcblrwfSlope-( 4.344515E-06+trapMax*-1.262E-05)-5.426306E-05;
    if(channel == 626) return nlcblrwfSlope-( 2.578407E-05+trapMax*-1.298E-05)-1.125461E-04;
    if(channel == 627) return nlcblrwfSlope-( 1.913080E-05+trapMax*-1.296E-05)-5.166456E-05;
    if(channel == 628) return nlcblrwfSlope-( 3.266857E-05+trapMax*-1.297E-05)-1.513889E-04;
    if(channel == 629) return nlcblrwfSlope-( 2.171356E-05+trapMax*-1.295E-05)-7.207791E-05;
    if(channel == 640) return nlcblrwfSlope-( 7.818519E-05+trapMax*-1.296E-05)-1.396431E-04;
    if(channel == 641) return nlcblrwfSlope-( 3.668705E-05+trapMax*-1.296E-05)-6.263415E-05;
    if(channel == 642) return nlcblrwfSlope-( 5.146852E-05+trapMax*-1.373E-05)-1.099235E-04;
    if(channel == 643) return nlcblrwfSlope-( 2.370897E-05+trapMax*-1.375E-05)-5.128737E-05;
    if(channel == 644) return nlcblrwfSlope-( 6.066896E-05+trapMax*-1.308E-05)-1.573303E-04;
    if(channel == 645) return nlcblrwfSlope-(-1.188733E-05+trapMax*-1.305E-05)-5.677784E-05;
    if(channel == 646) return nlcblrwfSlope-( 7.972405E-05+trapMax*-1.291E-05)-1.042533E-04;
    if(channel == 647) return nlcblrwfSlope-( 2.540927E-05+trapMax*-1.286E-05)-5.084706E-05;
    if(channel == 656) return nlcblrwfSlope-( 1.262340E-05+trapMax*-1.337E-05)-9.601873E-05;
    if(channel == 657) return nlcblrwfSlope-( 1.086301E-05+trapMax*-1.342E-05)-4.681626E-05;
    if(channel == 662) return nlcblrwfSlope-( 1.485632E-04+trapMax*-1.460E-05)-2.948478E-04;
    if(channel == 663) return nlcblrwfSlope-( 2.742351E-06+trapMax*-1.457E-05)-6.062255E-05;
    if(channel == 664) return nlcblrwfSlope-( 5.110392E-05+trapMax*-1.321E-05)-9.675150E-05;
    if(channel == 665) return nlcblrwfSlope-( 6.499560E-06+trapMax*-1.319E-05)-4.822259E-05;
    if(channel == 674) return nlcblrwfSlope-( 9.983398E-05+trapMax*-1.310E-05)-2.173713E-04;
    if(channel == 675) return nlcblrwfSlope-( 3.757473E-05+trapMax*-1.309E-05)-8.310702E-05;
    if(channel == 688) return nlcblrwfSlope-( 2.461617E-05+trapMax*-1.287E-05)-2.461283E-04;
    if(channel == 689) return nlcblrwfSlope-( 3.018505E-05+trapMax*-1.287E-05)-8.321697E-05;
    if(channel == 690) return nlcblrwfSlope-( 9.944427E-05+trapMax*-1.300E-05)-1.720172E-04;
    if(channel == 691) return nlcblrwfSlope-(-4.737995E-06+trapMax*-1.300E-05)-6.872953E-05;
    if(channel == 692) return nlcblrwfSlope-( 8.146529E-05+trapMax*-1.354E-05)-1.009845E-04;
    if(channel == 693) return nlcblrwfSlope-( 1.839544E-05+trapMax*-1.356E-05)-5.797197E-05;
    if(channel == 696) return nlcblrwfSlope-( 2.093490E-05+trapMax*-1.308E-05)-1.537598E-04;
    if(channel == 697) return nlcblrwfSlope-( 8.938457E-06+trapMax*-1.310E-05)-7.108055E-05;
    if(trapMax == 0) return 0;
    return nlcblrwfSlope/trapMax;
  }


  if(dsNumber == 1) {
    //params from 11507-11592 Updated 29 Aug 2016
	if(channel == 578) return nlcblrwfSlope-(4.620083E-05+trapMax*-1.306E-05)-1.414597E-04 ;
	if(channel == 579) return nlcblrwfSlope-(1.163801E-05+trapMax*-1.304E-05)-6.024276E-05 ;
	if(channel == 580) return nlcblrwfSlope-(2.604215E-05+trapMax*-1.291E-05)-2.600445E-04 ;
	if(channel == 581) return nlcblrwfSlope-(2.994274E-05+trapMax*-1.290E-05)-7.895203E-05 ;
	if(channel == 582) return nlcblrwfSlope-(5.175215E-05+trapMax*-1.304E-05)-8.452177E-05 ;
	if(channel == 583) return nlcblrwfSlope-(1.902729E-05+trapMax*-1.302E-05)-5.807816E-05 ;
	if(channel == 592) return nlcblrwfSlope-(5.145000E-05+trapMax*-1.297E-05)-1.158769E-04 ;
	if(channel == 593) return nlcblrwfSlope-(4.298575E-05+trapMax*-1.297E-05)-6.062252E-05 ;
	if(channel == 594) return nlcblrwfSlope-(7.298870E-05+trapMax*-1.353E-05)-1.216112E-04 ;
	if(channel == 595) return nlcblrwfSlope-(6.178016E-06+trapMax*-1.351E-05)-6.158457E-05 ;
	if(channel == 598) return nlcblrwfSlope-(3.497582E-05+trapMax*-1.298E-05)-1.198058E-04 ;
	if(channel == 599) return nlcblrwfSlope-(9.553672E-06+trapMax*-1.298E-05)-5.338186E-05 ;
	if(channel == 600) return nlcblrwfSlope-(5.649999E-05+trapMax*-1.287E-05)-9.150926E-05 ;
	if(channel == 601) return nlcblrwfSlope-(5.482757E-05+trapMax*-1.288E-05)-4.738173E-05 ;
	if(channel == 608) return nlcblrwfSlope-(4.586047E-05+trapMax*-1.284E-05)-2.255628E-04 ;
	if(channel == 609) return nlcblrwfSlope-(1.016108E-05+trapMax*-1.282E-05)-8.085080E-05 ;
	if(channel == 610) return nlcblrwfSlope-(2.209244E-05+trapMax*-1.298E-05)-1.144211E-04 ;
	if(channel == 611) return nlcblrwfSlope-(3.159581E-05+trapMax*-1.299E-05)-5.614616E-05 ;
	if(channel == 616) return nlcblrwfSlope-(7.722590E-05+trapMax*-1.316E-05)-1.175098E-04 ;
	if(channel == 617) return nlcblrwfSlope-(6.568340E-05+trapMax*-1.309E-05)-3.279914E-05 ;
	if(channel == 626) return nlcblrwfSlope-(1.786671E-05+trapMax*-1.266E-05)-1.285227E-04 ;
	if(channel == 627) return nlcblrwfSlope-(8.580088E-06+trapMax*-1.264E-05)-5.442935E-05 ;
	if(channel == 632) return nlcblrwfSlope-(9.481784E-05+trapMax*-1.307E-05)-1.709287E-04 ;
	if(channel == 633) return nlcblrwfSlope-(6.909768E-06+trapMax*-1.304E-05)-8.019256E-05 ;
	if(channel == 640) return nlcblrwfSlope-(7.977354E-05+trapMax*-1.287E-05)-9.866131E-05 ;
	if(channel == 641) return nlcblrwfSlope-(3.680967E-05+trapMax*-1.287E-05)-5.716358E-05 ;
	if(channel == 648) return nlcblrwfSlope-(9.035880E-05+trapMax*-1.312E-05)-1.874037E-04 ;
	if(channel == 649) return nlcblrwfSlope-(6.199340E-06+trapMax*-1.308E-05)-8.963580E-05 ;
	if(channel == 664) return nlcblrwfSlope-(5.873798E-05+trapMax*-1.321E-05)-1.320418E-04 ;
	if(channel == 665) return nlcblrwfSlope-(2.416351E-05+trapMax*-1.320E-05)-5.207803E-05 ;
	if(channel == 672) return nlcblrwfSlope-(2.815562E-05+trapMax*-1.328E-05)-1.103071E-04 ;
	if(channel == 673) return nlcblrwfSlope-(4.809700E-05+trapMax*-1.334E-05)-5.966680E-05 ;
	if(channel == 690) return nlcblrwfSlope-(4.148262E-05+trapMax*-1.308E-05)-8.145841E-05 ;
	if(channel == 691) return nlcblrwfSlope-(1.779754E-05+trapMax*-1.310E-05)-4.578676E-05 ;
	if(channel == 692) return nlcblrwfSlope-(4.737009E-05+trapMax*-1.309E-05)-1.106426E-04 ;
	if(channel == 693) return nlcblrwfSlope-(2.395867E-05+trapMax*-1.312E-05)-5.510341E-05 ;
    if(trapMax == 0) return 0;
    return nlcblrwfSlope/trapMax;
  }
  if(dsNumber == 3) {
    //params updated 26 Jan 2017, using updated DS 3 AvsE
    if (channel == 578) return (nlcblrwfSlope-(4.969765E-05+trapMax*-3.202149E-05)-1.305706E-04);
    if (channel == 579) return (nlcblrwfSlope-(2.112294E-05+trapMax*-9.527327E-06)-5.421336E-05);
    if (channel == 580) return (nlcblrwfSlope-(3.216003E-05+trapMax*-3.177636E-05)-2.416331E-04);
    if (channel == 581) return (nlcblrwfSlope-(1.605000E-05+trapMax*-9.465058E-06)-8.618267E-05);
    if (channel == 582) return (nlcblrwfSlope-(8.352656E-05+trapMax*-3.424562E-05)-1.445981E-04);
    if (channel == 592) return (nlcblrwfSlope-(1.551853E-04+trapMax*-3.083398E-05)-2.846104E-04);
    if (channel == 593) return (nlcblrwfSlope-(7.158773E-05+trapMax*-9.365700E-06)-9.634583E-05);
    if (channel == 594) return (nlcblrwfSlope-(1.697975E-04+trapMax*-3.371482E-05)-9.725605E-05);
    if (channel == 598) return (nlcblrwfSlope-(4.847808E-05+trapMax*-3.193175E-05)-1.172522E-04);
    if (channel == 599) return (nlcblrwfSlope-(2.946223E-05+trapMax*-9.478967E-06)-5.511747E-05);
    if (channel == 600) return (nlcblrwfSlope-(6.744091E-05+trapMax*-3.101750E-05)-1.006691E-04);
    if (channel == 601) return (nlcblrwfSlope-(5.674194E-05+trapMax*-9.176644E-06)-5.797139E-05);
    if (channel == 608) return (nlcblrwfSlope-(6.195887E-05+trapMax*-3.193785E-05)-2.350848E-04);
    if (channel == 609) return (nlcblrwfSlope-(1.194460E-05+trapMax*-9.356596E-06)-8.063594E-05);
    if (channel == 610) return (nlcblrwfSlope-(4.300218E-05+trapMax*-3.172544E-05)-1.282914E-04);
    if (channel == 611) return (nlcblrwfSlope-(3.435305E-05+trapMax*-9.473748E-06)-5.709263E-05);
    if (channel == 614) return (nlcblrwfSlope-(6.739594E-05+trapMax*-3.347227E-05)-1.735536E-04);
    if (channel == 615) return (nlcblrwfSlope-(2.487351E-05+trapMax*-9.882128E-06)-5.818863E-05);
    if (channel == 624) return (nlcblrwfSlope-(1.037784E-04+trapMax*-3.371539E-05)-1.151629E-04);
    if (channel == 625) return (nlcblrwfSlope-(-1.369093E-05+trapMax*-9.824005E-06)-5.430976E-05);
    if (channel == 626) return (nlcblrwfSlope-(4.360833E-05+trapMax*-3.177044E-05)-1.047656E-04);
    if (channel == 627) return (nlcblrwfSlope-(1.045450E-05+trapMax*-9.387932E-06)-4.929234E-05);
    if (channel == 632) return (nlcblrwfSlope-(1.064173E-04+trapMax*-3.163814E-05)-1.761009E-04);
    if (channel == 633) return (nlcblrwfSlope-(4.261676E-06+trapMax*-9.440554E-06)-7.412162E-05);
    if (channel == 640) return (nlcblrwfSlope-(8.553435E-05+trapMax*-3.203242E-05)-1.081273E-04);
    if (channel == 641) return (nlcblrwfSlope-(3.859470E-05+trapMax*-9.443690E-06)-5.198290E-05);
    if (channel == 648) return (nlcblrwfSlope-(1.531563E-04+trapMax*-3.265891E-05)-1.732467E-04);
    if (channel == 649) return (nlcblrwfSlope-(2.349817E-05+trapMax*-9.703813E-06)-8.717817E-05);
    if (channel == 664) return (nlcblrwfSlope-(7.833005E-05+trapMax*-3.202713E-05)-1.282682E-04);
    if (channel == 665) return (nlcblrwfSlope-(3.095469E-05+trapMax*-9.595590E-06)-5.746954E-05);
    if (channel == 672) return (nlcblrwfSlope-(3.802716E-05+trapMax*-3.267933E-05)-9.456284E-05);
    if (channel == 673) return (nlcblrwfSlope-(2.906844E-05+trapMax*-9.668007E-06)-4.921254E-05);
    if (channel == 678) return (nlcblrwfSlope-(1.237485E-04+trapMax*-3.675035E-05)-2.355340E-04);
    if (channel == 679) return (nlcblrwfSlope-(7.193220E-05+trapMax*-1.097365E-05)-9.160413E-05);
    if (channel == 690) return (nlcblrwfSlope-(3.250624E-05+trapMax*-3.214250E-05)-9.273793E-05);
    if (channel == 691) return (nlcblrwfSlope-(1.172497E-05+trapMax*-9.542865E-06)-4.658490E-05);
    if (channel == 692) return (nlcblrwfSlope-(6.658605E-05+trapMax*-3.263892E-05)-1.908410E-04);
    if (channel == 693) return (nlcblrwfSlope-(3.125473E-05+trapMax*-9.809244E-06)-6.439295E-05);
    if (channel == 694) return (nlcblrwfSlope-(6.968617E-05+trapMax*-3.409197E-05)-1.264238E-04);
    if (channel == 695) return (nlcblrwfSlope-(5.918364E-05+trapMax*-1.009706E-05)-4.901754E-05);

   if(trapMax == 0) return 0;
    return nlcblrwfSlope/trapMax;
  }
  if(dsNumber == 4) {
     if (channel == 1106) return (nlcblrwfSlope-(1.043332E-04+trapMax*-3.117063E-05)-1.611705E-04 );
     if (channel == 1107) return (nlcblrwfSlope-(6.943577E-05+trapMax*-9.261796E-06)-5.523397E-05 );
     if (channel == 1136) return (nlcblrwfSlope-(1.023021E-04+trapMax*-3.118825E-05)-9.249370E-05 );
     if (channel == 1137) return (nlcblrwfSlope-(4.532656E-05+trapMax*-9.245771E-06)-3.782950E-05 );
     if (channel == 1144) return (nlcblrwfSlope-(1.731251E-04+trapMax*-3.102642E-05)-7.360783E-05 );
     if (channel == 1145) return (nlcblrwfSlope-(7.414058E-05+trapMax*-9.344486E-06)-3.217599E-05 );
     if (channel == 1170) return (nlcblrwfSlope-(1.208446E-04+trapMax*-3.024511E-05)-9.963974E-05 );
     if (channel == 1171) return (nlcblrwfSlope-(5.247239E-05+trapMax*-8.995409E-06)-5.767731E-05 );
     if (channel == 1172) return (nlcblrwfSlope-(8.341973E-05+trapMax*-3.119219E-05)-2.706324E-04 );
     if (channel == 1173) return (nlcblrwfSlope-(4.997503E-05+trapMax*-9.378188E-06)-9.747403E-05 );
     if (channel == 1174) return (nlcblrwfSlope-(8.231800E-05+trapMax*-3.233127E-05)-1.153748E-04 );
     if (channel == 1176) return (nlcblrwfSlope-(7.421930E-05+trapMax*-3.133153E-05)-1.194989E-04 );
     if (channel == 1177) return (nlcblrwfSlope-(4.489232E-05+trapMax*-9.351887E-06)-4.731424E-05 );
     if (channel == 1204) return (nlcblrwfSlope-(3.585310E-05+trapMax*-3.245930E-05)-1.151355E-04 );
     if (channel == 1205) return (nlcblrwfSlope-(5.514260E-05+trapMax*-9.588399E-06)-5.293414E-05 );
     if (channel == 1232) return (nlcblrwfSlope-(8.035306E-05+trapMax*-3.229609E-05)-8.604362E-05 );
     if (channel == 1233) return (nlcblrwfSlope-(6.434726E-05+trapMax*-9.593176E-06)-4.550998E-05 );
     if (channel == 1236) return (nlcblrwfSlope-(3.372169E-05+trapMax*-3.287308E-05)-8.864796E-05 );
     if (channel == 1237) return (nlcblrwfSlope-(5.888586E-05+trapMax*-9.776967E-06)-5.218772E-05 );
     if (channel == 1296) return (nlcblrwfSlope-(2.948720E-04+trapMax*-3.244349E-05)-4.899846E-05 );
     if (channel == 1297) return (nlcblrwfSlope-(2.062004E-04+trapMax*-9.777431E-06)-2.899245E-05 );
     if (channel == 1298) return (nlcblrwfSlope-(1.125306E-05+trapMax*-3.216418E-05)-9.907714E-05 );
     if (channel == 1299) return (nlcblrwfSlope-(4.456495E-05+trapMax*-9.406765E-06)-5.186341E-05 );
     if (channel == 1330) return (nlcblrwfSlope-(1.488386E-04+trapMax*-3.124473E-05)-9.610995E-05 );
     if (channel == 1331) return (nlcblrwfSlope-(8.532411E-05+trapMax*-9.369605E-06)-5.126190E-05 );
     if (channel == 1332) return (nlcblrwfSlope-(5.693212E-05+trapMax*-3.189530E-05)-1.762632E-04 );
     if (channel == 1333) return (nlcblrwfSlope-(6.527037E-05+trapMax*-9.668432E-06)-7.727936E-05 );


if(trapMax == 0) return 0;
    return nlcblrwfSlope/trapMax;
  }


  else //cout << "GetDCR90(): unknown dataset number DS" << dsNumber << endl;
  return 0;
}


double GetDCRCTC(int channel, double nlcblrwfSlope, double trapE, double trapMax, int dsNumber) //DCR with charge trapping correction applied
{
  if(dsNumber == 3) {
    //params updated 26 Jan 2017, using updated DS 3 AvsE
    if (channel == 578) return (nlcblrwfSlope-((1.127394E-04*exp(-2.433297E-04*trapMax))*(trapE-trapMax))-(6.153314E-05+trapMax*-3.202111E-05)-1.150849E-04);
    if (channel == 579) return (nlcblrwfSlope-((7.235269E-05*exp(-5.314217E-04*trapMax))*(trapE-trapMax))-(1.867161E-05+trapMax*-9.518036E-06)-4.930772E-05);
    if (channel == 580) return (nlcblrwfSlope-((1.738967E-04*exp(-1.606743E-04*trapMax))*(trapE-trapMax))-(4.458943E-05+trapMax*-3.179558E-05)-1.454872E-04);
    if (channel == 581) return (nlcblrwfSlope-((6.522808E-05*exp(-1.579528E-04*trapMax))*(trapE-trapMax))-(4.755118E-05+trapMax*-9.476431E-06)-5.590125E-05);
    if (channel == 582) return (nlcblrwfSlope-((1.235261E-04*exp(-3.158980E-04*trapMax))*(trapE-trapMax))-(4.344448E-05+trapMax*-3.423551E-05)-1.325392E-04);
    if (channel == 592) return (nlcblrwfSlope-((2.108174E-04*exp(-1.846621E-04*trapMax))*(trapE-trapMax))-(2.330795E-05+trapMax*-3.082218E-05)-1.835619E-04);
    if (channel == 593) return (nlcblrwfSlope-((9.618042E-05*exp(-2.638546E-04*trapMax))*(trapE-trapMax))-(3.677377E-05+trapMax*-9.338063E-06)-7.914180E-05);
    if (channel == 594) return (nlcblrwfSlope-((1.766883E-04*exp(-6.883060E-04*trapMax))*(trapE-trapMax))-(1.569557E-04+trapMax*-3.370489E-05)-1.086235E-04);
    if (channel == 598) return (nlcblrwfSlope-((1.538365E-04*exp(-3.788839E-04*trapMax))*(trapE-trapMax))-(4.309318E-05+trapMax*-3.192628E-05)-1.158028E-04);
    if (channel == 599) return (nlcblrwfSlope-((5.331744E-05*exp(-2.472075E-04*trapMax))*(trapE-trapMax))-(3.662790E-05+trapMax*-9.482648E-06)-5.251377E-05);
    if (channel == 600) return (nlcblrwfSlope-((3.970803E-04*exp(-6.189204E-04*trapMax))*(trapE-trapMax))-(9.320193E-05+trapMax*-3.101883E-05)-9.693572E-05);
    if (channel == 601) return (nlcblrwfSlope-((1.355323E-04*exp(-4.784889E-04*trapMax))*(trapE-trapMax))-(-1.547546E-05+trapMax*-9.212459E-06)-5.789839E-05);
    if (channel == 608) return (nlcblrwfSlope-((1.513094E-04*exp(-2.056450E-04*trapMax))*(trapE-trapMax))-(4.419241E-05+trapMax*-3.198043E-05)-1.642101E-04);
    if (channel == 609) return (nlcblrwfSlope-((5.864419E-05*exp(-2.491374E-04*trapMax))*(trapE-trapMax))-(1.282546E-05+trapMax*-9.367739E-06)-5.935112E-05);
    if (channel == 610) return (nlcblrwfSlope-((1.635264E-04*exp(-5.584308E-04*trapMax))*(trapE-trapMax))-(2.390021E-06+trapMax*-3.172120E-05)-1.147398E-04);
    if (channel == 611) return (nlcblrwfSlope-((5.995284E-05*exp(-4.999357E-04*trapMax))*(trapE-trapMax))-(3.654757E-05+trapMax*-9.470468E-06)-5.311655E-05);
    if (channel == 614) return (nlcblrwfSlope-((1.571440E-04*exp(-2.201333E-04*trapMax))*(trapE-trapMax))-(9.170683E-05+trapMax*-3.349689E-05)-1.526042E-04);
    if (channel == 615) return (nlcblrwfSlope-((5.917700E-05*exp(-2.414154E-04*trapMax))*(trapE-trapMax))-(3.863961E-05+trapMax*-9.895673E-06)-5.389743E-05);
    if (channel == 624) return (nlcblrwfSlope-((2.698744E-04*exp(-5.246425E-04*trapMax))*(trapE-trapMax))-(1.073635E-04+trapMax*-3.372771E-05)-1.188775E-04);
    if (channel == 625) return (nlcblrwfSlope-((1.145199E-04*exp(-4.571575E-04*trapMax))*(trapE-trapMax))-(2.245688E-05+trapMax*-9.844269E-06)-5.938712E-05);
    if (channel == 626) return (nlcblrwfSlope-((9.943684E-05*exp(-4.822379E-04*trapMax))*(trapE-trapMax))-(6.128326E-05+trapMax*-3.176414E-05)-1.075445E-04);
    if (channel == 627) return (nlcblrwfSlope-((3.446802E-05*exp(-4.026681E-04*trapMax))*(trapE-trapMax))-(2.297720E-05+trapMax*-9.381861E-06)-5.024683E-05);
    if (channel == 632) return (nlcblrwfSlope-((1.729608E-04*exp(-3.541229E-04*trapMax))*(trapE-trapMax))-(-3.290756E-05+trapMax*-3.161975E-05)-1.421823E-04);
    if (channel == 633) return (nlcblrwfSlope-((5.123398E-05*exp(-2.273723E-04*trapMax))*(trapE-trapMax))-(1.023045E-06+trapMax*-9.457910E-06)-6.580567E-05);
    if (channel == 640) return (nlcblrwfSlope-((2.069584E-04*exp(-4.945272E-04*trapMax))*(trapE-trapMax))-(7.722267E-05+trapMax*-3.203165E-05)-1.161178E-04);
    if (channel == 641) return (nlcblrwfSlope-((6.346157E-05*exp(-3.367259E-04*trapMax))*(trapE-trapMax))-(3.714462E-05+trapMax*-9.444573E-06)-5.274827E-05);
    if (channel == 648) return (nlcblrwfSlope-((1.725229E-04*exp(-3.497286E-04*trapMax))*(trapE-trapMax))-(-1.395188E-05+trapMax*-3.264036E-05)-1.254127E-04);
    if (channel == 649) return (nlcblrwfSlope-((5.430675E-05*exp(-2.314091E-04*trapMax))*(trapE-trapMax))-(-8.944687E-06+trapMax*-9.703637E-06)-7.539341E-05);
    if (channel == 664) return (nlcblrwfSlope-((8.734885E-05*exp(-3.937745E-04*trapMax))*(trapE-trapMax))-(8.975668E-05+trapMax*-3.202976E-05)-1.243938E-04);
    if (channel == 665) return (nlcblrwfSlope-((2.433784E-05*exp(-2.586932E-04*trapMax))*(trapE-trapMax))-(3.979345E-05+trapMax*-9.595988E-06)-5.574878E-05);
    if (channel == 672) return (nlcblrwfSlope-((2.563171E-04*exp(-6.349953E-04*trapMax))*(trapE-trapMax))-(6.899191E-05+trapMax*-3.267645E-05)-1.020221E-04);
    if (channel == 673) return (nlcblrwfSlope-((4.690967E-05*exp(-1.717128E-04*trapMax))*(trapE-trapMax))-(2.536466E-05+trapMax*-9.654297E-06)-4.957713E-05);
    if (channel == 678) return (nlcblrwfSlope-((1.411855E-04*exp(-2.147832E-04*trapMax))*(trapE-trapMax))-(3.559417E-05+trapMax*-3.674342E-05)-1.590019E-04);
    if (channel == 679) return (nlcblrwfSlope-((6.143846E-05*exp(-2.302520E-04*trapMax))*(trapE-trapMax))-(2.409803E-05+trapMax*-1.095705E-05)-5.802272E-05);
    if (channel == 690) return (nlcblrwfSlope-((2.192544E-04*exp(-8.526074E-04*trapMax))*(trapE-trapMax))-(8.128862E-05+trapMax*-3.215197E-05)-1.009718E-04);
    if (channel == 691) return (nlcblrwfSlope-((1.295922E-04*exp(-1.067913E-03*trapMax))*(trapE-trapMax))-(6.163567E-05+trapMax*-9.559764E-06)-4.852831E-05);
    if (channel == 692) return (nlcblrwfSlope-((2.372333E-04*exp(-3.010077E-04*trapMax))*(trapE-trapMax))-(8.859344E-05+trapMax*-3.262989E-05)-1.782932E-04);
    if (channel == 693) return (nlcblrwfSlope-((1.596905E-04*exp(-5.558878E-04*trapMax))*(trapE-trapMax))-(5.573115E-05+trapMax*-9.830307E-06)-6.304652E-05);
    if (channel == 694) return (nlcblrwfSlope-((2.179297E-04*exp(-3.221275E-04*trapMax))*(trapE-trapMax))-(6.854911E-05+trapMax*-3.408329E-05)-1.182403E-04);
    if (channel == 695) return (nlcblrwfSlope-((8.287052E-05*exp(-3.267436E-04*trapMax))*(trapE-trapMax))-(5.095433E-05+trapMax*-1.010929E-05)-4.731396E-05);

  if(trapMax == 0) return 0;
    return nlcblrwfSlope/trapMax;
  }
  if(dsNumber == 4) {
     if (channel == 1106) return (nlcblrwfSlope-((2.280664E-04*exp(-3.005270E-04*trapMax))*(trapE-trapMax))-(1.585164E-04+trapMax*-3.110136E-05)-1.319380E-04);
     if (channel == 1107) return (nlcblrwfSlope-((6.735549E-05*exp(-2.409237E-04*trapMax))*(trapE-trapMax))-(8.310126E-05+trapMax*-9.257286E-06)-4.639398E-05);
     if (channel == 1136) return (nlcblrwfSlope-((2.432154E-04*exp(-3.340925E-04*trapMax))*(trapE-trapMax))-(1.373025E-04+trapMax*-3.117127E-05)-9.067373E-05);
     if (channel == 1137) return (nlcblrwfSlope-((1.155699E-04*exp(-4.071065E-04*trapMax))*(trapE-trapMax))-(5.520871E-05+trapMax*-9.248161E-06)-3.814742E-05);
     if (channel == 1144) return (nlcblrwfSlope-((3.996635E-04*exp(-5.550602E-04*trapMax))*(trapE-trapMax))-(2.612822E-04+trapMax*-3.102647E-05)-8.407006E-05);
     if (channel == 1145) return (nlcblrwfSlope-((2.510294E-04*exp(-6.922776E-04*trapMax))*(trapE-trapMax))-(1.742175E-04+trapMax*-9.358991E-06)-3.928366E-05);
     if (channel == 1170) return (nlcblrwfSlope-((2.634281E-04*exp(-3.208118E-04*trapMax))*(trapE-trapMax))-(1.680486E-04+trapMax*-3.024031E-05)-9.169512E-05);
     if (channel == 1171) return (nlcblrwfSlope-((1.558809E-04*exp(-5.037593E-04*trapMax))*(trapE-trapMax))-(1.080515E-04+trapMax*-9.002678E-06)-5.471471E-05);
     if (channel == 1172) return (nlcblrwfSlope-((2.294414E-04*exp(-2.090438E-04*trapMax))*(trapE-trapMax))-(1.313878E-04+trapMax*-3.115949E-05)-1.516341E-04);
     if (channel == 1173) return (nlcblrwfSlope-((8.541671E-05*exp(-2.124352E-04*trapMax))*(trapE-trapMax))-(8.131950E-05+trapMax*-9.374774E-06)-7.122980E-05);
     if (channel == 1174) return (nlcblrwfSlope-((2.338199E-04*exp(-1.781426E-04*trapMax))*(trapE-trapMax))-(7.496339E-05+trapMax*-3.228964E-05)-1.141873E-04);
     if (channel == 1176) return (nlcblrwfSlope-((-1.430449E-08*exp(3.443593E-03*trapMax))*(trapE-trapMax))-(8.169231E-05+trapMax*-3.133717E-05)-1.163027E-04);
     if (channel == 1177) return (nlcblrwfSlope-((1.582782E-04*exp(-1.103157E-03*trapMax))*(trapE-trapMax))-(5.857280E-05+trapMax*-9.356245E-06)-4.853332E-05);
     if (channel == 1204) return (nlcblrwfSlope-((1.290170E-04*exp(-1.577457E-04*trapMax))*(trapE-trapMax))-(5.987632E-05+trapMax*-3.245611E-05)-1.105412E-04);
     if (channel == 1205) return (nlcblrwfSlope-((5.155111E-05*exp(-1.882470E-04*trapMax))*(trapE-trapMax))-(6.155698E-05+trapMax*-9.583817E-06)-5.184528E-05);
     if (channel == 1232) return (nlcblrwfSlope-((5.862128E-04*exp(-8.824675E-04*trapMax))*(trapE-trapMax))-(1.358949E-04+trapMax*-3.231454E-05)-9.651307E-05);
     if (channel == 1233) return (nlcblrwfSlope-((2.021734E-04*exp(-6.922000E-04*trapMax))*(trapE-trapMax))-(7.511586E-05+trapMax*-9.600501E-06)-4.866403E-05);
     if (channel == 1236) return (nlcblrwfSlope-((5.169144E-04*exp(-1.796068E-03*trapMax))*(trapE-trapMax))-(3.933840E-05+trapMax*-3.287496E-05)-8.986392E-05);
     if (channel == 1237) return (nlcblrwfSlope-((5.053424E-05*exp(-6.511982E-04*trapMax))*(trapE-trapMax))-(5.594828E-05+trapMax*-9.776989E-06)-5.220814E-05);
     if (channel == 1296) return (nlcblrwfSlope-((7.507806E-04*exp(-4.953029E-04*trapMax))*(trapE-trapMax))-(2.479104E-04+trapMax*-3.237509E-05)-1.073754E-04);
     if (channel == 1297) return (nlcblrwfSlope-((5.060572E-04*exp(-6.827946E-04*trapMax))*(trapE-trapMax))-(4.391638E-05+trapMax*-9.749360E-06)-5.963386E-05);
     if (channel == 1298) return (nlcblrwfSlope-((1.517395E-04*exp(-5.723059E-04*trapMax))*(trapE-trapMax))-(1.302316E-05+trapMax*-3.216506E-05)-1.000615E-04);
     if (channel == 1299) return (nlcblrwfSlope-((4.137805E-05*exp(-4.922765E-04*trapMax))*(trapE-trapMax))-(5.215931E-05+trapMax*-9.408508E-06)-5.167422E-05);
     if (channel == 1330) return (nlcblrwfSlope-((1.008143E-03*exp(-3.838487E-04*trapMax))*(trapE-trapMax))-(1.565639E-04+trapMax*-3.112789E-05)-1.154411E-04);
     if (channel == 1331) return (nlcblrwfSlope-((3.287850E-04*exp(-3.059999E-04*trapMax))*(trapE-trapMax))-(1.051466E-05+trapMax*-9.332279E-06)-5.966420E-05);
     if (channel == 1332) return (nlcblrwfSlope-((1.720117E-04*exp(-2.109293E-04*trapMax))*(trapE-trapMax))-(7.220260E-05+trapMax*-3.190263E-05)-1.329770E-04);
     if (channel == 1333) return (nlcblrwfSlope-((6.012840E-05*exp(-1.833330E-04*trapMax))*(trapE-trapMax))-(6.935094E-05+trapMax*-9.671974E-06)-6.534142E-05);


if(trapMax == 0) return 0;
    return nlcblrwfSlope/trapMax;
  }


  else //cout << "GetDCRCTC(): unknown dataset number DS" << dsNumber << endl;
  return 0;
}

double GetDCR95(int channel, double nlcblrwfSlope, double trapMax, int dsNumber)
{
  if(dsNumber == 0) {
    if(channel == 576) return nlcblrwfSlope-( 3.673385E-05+trapMax*-1.279E-05)-1.438638E-04;
    if(channel == 577) return nlcblrwfSlope-( 3.122161E-05+trapMax*-1.278E-05)-6.627588E-05;
    if(channel == 592) return nlcblrwfSlope-( 1.371969E-04+trapMax*-1.287E-05)-1.777946E-04;
    if(channel == 593) return nlcblrwfSlope-( 4.661555E-05+trapMax*-1.285E-05)-7.428328E-05;
    if(channel == 594) return nlcblrwfSlope-( 7.684589E-05+trapMax*-1.294E-05)-1.487208E-04;
    if(channel == 595) return nlcblrwfSlope-( 5.906154E-05+trapMax*-1.294E-05)-6.560241E-05;
    if(channel == 598) return nlcblrwfSlope-( 1.842168E-04+trapMax*-1.315E-05)-1.310455E-04;
    if(channel == 599) return nlcblrwfSlope-( 7.434566E-05+trapMax*-1.315E-05)-6.411484E-05;
    if(channel == 600) return nlcblrwfSlope-( 5.645349E-05+trapMax*-1.264E-05)-1.709638E-04;
    if(channel == 601) return nlcblrwfSlope-( 1.644013E-05+trapMax*-1.260E-05)-9.198958E-05;
    if(channel == 608) return nlcblrwfSlope-( 9.091576E-05+trapMax*-1.298E-05)-1.461563E-04;
    if(channel == 609) return nlcblrwfSlope-( 1.319922E-05+trapMax*-1.296E-05)-5.850014E-05;
    if(channel == 610) return nlcblrwfSlope-( 3.708218E-05+trapMax*-1.331E-05)-2.135500E-04;
    if(channel == 611) return nlcblrwfSlope-( 1.691204E-05+trapMax*-1.331E-05)-9.026296E-05;
    if(channel == 614) return nlcblrwfSlope-( 7.716150E-05+trapMax*-1.272E-05)-3.157850E-04;
    if(channel == 615) return nlcblrwfSlope-( 1.479029E-05+trapMax*-1.271E-05)-1.252092E-04;
    if(channel == 624) return nlcblrwfSlope-( 3.891078E-05+trapMax*-1.262E-05)-1.609082E-04;
    if(channel == 625) return nlcblrwfSlope-( 4.344515E-06+trapMax*-1.262E-05)-7.434378E-05;
    if(channel == 626) return nlcblrwfSlope-( 2.578407E-05+trapMax*-1.298E-05)-1.518572E-04;
    if(channel == 627) return nlcblrwfSlope-( 1.913080E-05+trapMax*-1.296E-05)-6.920177E-05;
    if(channel == 628) return nlcblrwfSlope-( 3.266857E-05+trapMax*-1.297E-05)-2.094821E-04;
    if(channel == 629) return nlcblrwfSlope-( 2.171356E-05+trapMax*-1.295E-05)-9.346685E-05;
    if(channel == 640) return nlcblrwfSlope-( 7.818519E-05+trapMax*-1.296E-05)-1.868514E-04;
    if(channel == 641) return nlcblrwfSlope-( 3.668705E-05+trapMax*-1.296E-05)-8.219029E-05;
    if(channel == 642) return nlcblrwfSlope-( 5.146852E-05+trapMax*-1.373E-05)-1.453889E-04;
    if(channel == 643) return nlcblrwfSlope-( 2.370897E-05+trapMax*-1.375E-05)-6.690158E-05;
    if(channel == 644) return nlcblrwfSlope-( 6.066896E-05+trapMax*-1.308E-05)-2.107852E-04;
    if(channel == 645) return nlcblrwfSlope-(-1.188733E-05+trapMax*-1.305E-05)-8.136305E-05;
    if(channel == 646) return nlcblrwfSlope-( 7.972405E-05+trapMax*-1.291E-05)-1.473133E-04;
    if(channel == 647) return nlcblrwfSlope-( 2.540927E-05+trapMax*-1.286E-05)-3.482776E-04;
    if(channel == 656) return nlcblrwfSlope-( 1.262340E-05+trapMax*-1.337E-05)-1.287283E-04;
    if(channel == 657) return nlcblrwfSlope-( 1.086301E-05+trapMax*-1.342E-05)-6.244994E-05;
    if(channel == 662) return nlcblrwfSlope-( 1.485632E-04+trapMax*-1.460E-05)-3.979238E-04;
    if(channel == 663) return nlcblrwfSlope-( 2.742351E-06+trapMax*-1.457E-05)-8.256436E-05;
    if(channel == 664) return nlcblrwfSlope-( 5.110392E-05+trapMax*-1.321E-05)-1.323383E-04;
    if(channel == 665) return nlcblrwfSlope-( 6.499560E-06+trapMax*-1.319E-05)-6.696894E-05;
    if(channel == 674) return nlcblrwfSlope-( 9.983398E-05+trapMax*-1.310E-05)-3.899950E-04;
    if(channel == 675) return nlcblrwfSlope-( 3.757473E-05+trapMax*-1.309E-05)-1.327352E-04;
    if(channel == 688) return nlcblrwfSlope-( 2.461617E-05+trapMax*-1.287E-05)-3.258232E-04;
    if(channel == 689) return nlcblrwfSlope-( 3.018505E-05+trapMax*-1.287E-05)-1.088782E-04;
    if(channel == 690) return nlcblrwfSlope-( 9.944427E-05+trapMax*-1.300E-05)-2.461961E-04;
    if(channel == 691) return nlcblrwfSlope-(-4.737995E-06+trapMax*-1.300E-05)-9.352369E-05;
    if(channel == 692) return nlcblrwfSlope-( 8.146529E-05+trapMax*-1.354E-05)-1.456008E-04;
    if(channel == 693) return nlcblrwfSlope-( 1.839544E-05+trapMax*-1.356E-05)-7.560564E-05;
    if(channel == 696) return nlcblrwfSlope-( 2.093490E-05+trapMax*-1.308E-05)-2.075377E-04;
    if(channel == 697) return nlcblrwfSlope-( 8.938457E-06+trapMax*-1.310E-05)-9.389924E-05;
    if(trapMax == 0) return 0;
    return nlcblrwfSlope/trapMax;
  }

  if(dsNumber == 1) {
    //params from 11507-11592
    if(channel == 672) return nlcblrwfSlope-(2.833847E-05+trapMax*-1.328E-05)-1.418052E-04;
    if(channel == 673) return nlcblrwfSlope-(4.734834E-05+trapMax*-1.334E-05)-8.184476E-05;
    if(channel == 690) return nlcblrwfSlope-(3.690947E-05+trapMax*-1.308E-05)-1.221969E-04;
    if(channel == 691) return nlcblrwfSlope-(1.605341E-05+trapMax*-1.310E-05)-6.775107E-05;
    if(channel == 692) return nlcblrwfSlope-(5.681803E-05+trapMax*-1.309E-05)-1.519825E-04;
    if(channel == 693) return nlcblrwfSlope-(2.453714E-05+trapMax*-1.312E-05)-7.302743E-05;
    if(channel == 578) return nlcblrwfSlope-(4.417989E-05+trapMax*-1.306E-05)-1.928345E-04;
    if(channel == 579) return nlcblrwfSlope-(1.106159E-05+trapMax*-1.304E-05)-8.146071E-05;
    if(channel == 580) return nlcblrwfSlope-(2.577260E-05+trapMax*-1.291E-05)-3.260552E-04;
    if(channel == 581) return nlcblrwfSlope-(3.041374E-05+trapMax*-1.290E-05)-1.044327E-04;
    if(channel == 582) return nlcblrwfSlope-(5.396620E-05+trapMax*-1.304E-05)-1.747793E-04;
    if(channel == 583) return nlcblrwfSlope-(1.977025E-05+trapMax*-1.302E-05)-8.834911E-05; //average for all LG channels
    if(channel == 592) return nlcblrwfSlope-(6.803550E-05+trapMax*-1.297E-05)-2.718900E-04;
    if(channel == 593) return nlcblrwfSlope-(4.309706E-05+trapMax*-1.296E-05)-1.013371E-04;
    if(channel == 594) return nlcblrwfSlope-(4.643130E-05+trapMax*-1.352E-05)-1.443377E-04;
    if(channel == 595) return nlcblrwfSlope-(1.550955E-06+trapMax*-1.351E-05)-8.834911E-05; //average for all LG channels
    if(channel == 598) return nlcblrwfSlope-(2.897476E-05+trapMax*-1.298E-05)-1.715410E-04;
    if(channel == 599) return nlcblrwfSlope-(9.612110E-06+trapMax*-1.298E-05)-7.209174E-05;
    if(channel == 600) return nlcblrwfSlope-(5.842046E-05+trapMax*-1.287E-05)-1.289802E-04;
    if(channel == 601) return nlcblrwfSlope-(5.372889E-05+trapMax*-1.288E-05)-6.895361E-05;
    if(channel == 608) return nlcblrwfSlope-(4.693324E-05+trapMax*-1.284E-05)-3.084515E-04;
    if(channel == 609) return nlcblrwfSlope-(1.055481E-05+trapMax*-1.282E-05)-1.090966E-04;
    if(channel == 610) return nlcblrwfSlope-(1.216405E-05+trapMax*-1.298E-05)-1.758732E-04;
    if(channel == 611) return nlcblrwfSlope-(2.884038E-05+trapMax*-1.299E-05)-7.925371E-05;
    if(channel == 616) return nlcblrwfSlope-(8.295462E-05+trapMax*-1.316E-05)-1.587131E-04;
    if(channel == 617) return nlcblrwfSlope-(7.294466E-05+trapMax*-1.310E-05)-4.362600E-04;
    if(channel == 626) return nlcblrwfSlope-(1.860959E-06+trapMax*-1.265E-05)-1.412070E-04;
    if(channel == 627) return nlcblrwfSlope-(7.000642E-06+trapMax*-1.264E-05)-7.404711E-05;
    if(channel == 632) return nlcblrwfSlope-(9.473908E-05+trapMax*-1.307E-05)-2.461187E-04;
    if(channel == 633) return nlcblrwfSlope-(6.874599E-06+trapMax*-1.304E-05)-1.122029E-04;
    if(channel == 640) return nlcblrwfSlope-(8.088126E-05+trapMax*-1.287E-05)-1.427142E-04;
    if(channel == 641) return nlcblrwfSlope-(3.667507E-05+trapMax*-1.287E-05)-7.623847E-05;
    if(channel == 648) return nlcblrwfSlope-(8.317390E-05+trapMax*-1.311E-05)-1.863428E-04;
    if(channel == 649) return nlcblrwfSlope-(4.360668E-06+trapMax*-1.308E-05)-1.220285E-04;
    if(channel == 664) return nlcblrwfSlope-(5.916704E-05+trapMax*-1.321E-05)-1.822770E-04;
    if(channel == 665) return nlcblrwfSlope-(2.373893E-05+trapMax*-1.320E-05)-7.736052E-05;

    if(trapMax == 0) return 0;
    return nlcblrwfSlope/trapMax;
  }

  else //cout << "GetDCR95(): unknown dataset number DS" << dsNumber << endl;
  return 0;
}


double GetDCR98(int channel, double nlcblrwfSlope, double trapMax, int dsNumber)
{
  if(dsNumber == 0) {
    if(channel == 576) return nlcblrwfSlope-( 3.673385E-05+trapMax*-1.279E-05)-2.128491E-04;
    if(channel == 577) return nlcblrwfSlope-( 3.122161E-05+trapMax*-1.278E-05)-9.009598E-05;
    if(channel == 592) return nlcblrwfSlope-( 1.371969E-04+trapMax*-1.287E-05)-2.478023E-04;
    if(channel == 593) return nlcblrwfSlope-( 4.661555E-05+trapMax*-1.285E-05)-9.960543E-05;
    if(channel == 594) return nlcblrwfSlope-( 7.684589E-05+trapMax*-1.294E-05)-2.021753E-04;
    if(channel == 595) return nlcblrwfSlope-( 5.906154E-05+trapMax*-1.294E-05)-8.815796E-05;
    if(channel == 598) return nlcblrwfSlope-( 1.842168E-04+trapMax*-1.315E-05)-5.932375E-04;
    if(channel == 599) return nlcblrwfSlope-( 7.434566E-05+trapMax*-1.315E-05)-2.631672E-04;
    if(channel == 600) return nlcblrwfSlope-( 5.645349E-05+trapMax*-1.264E-05)-2.128397E-04;
    if(channel == 601) return nlcblrwfSlope-( 1.644013E-05+trapMax*-1.260E-05)-1.156453E-04;
    if(channel == 608) return nlcblrwfSlope-( 9.091576E-05+trapMax*-1.298E-05)-2.108129E-04;
    if(channel == 609) return nlcblrwfSlope-( 1.319922E-05+trapMax*-1.296E-05)-8.500168E-05;
    if(channel == 610) return nlcblrwfSlope-( 3.708218E-05+trapMax*-1.331E-05)-3.363783E-04;
    if(channel == 611) return nlcblrwfSlope-( 1.691204E-05+trapMax*-1.331E-05)-1.291137E-04;
    if(channel == 614) return nlcblrwfSlope-( 7.716150E-05+trapMax*-1.272E-05)-4.345738E-04;
    if(channel == 615) return nlcblrwfSlope-( 1.479029E-05+trapMax*-1.271E-05)-1.682398E-04;
    if(channel == 624) return nlcblrwfSlope-( 3.891078E-05+trapMax*-1.262E-05)-2.108711E-04;
    if(channel == 625) return nlcblrwfSlope-( 4.344515E-06+trapMax*-1.262E-05)-1.006748E-04;
    if(channel == 626) return nlcblrwfSlope-( 2.578407E-05+trapMax*-1.298E-05)-1.994223E-04;
    if(channel == 627) return nlcblrwfSlope-( 1.913080E-05+trapMax*-1.296E-05)-9.152186E-05;
    if(channel == 628) return nlcblrwfSlope-( 3.266857E-05+trapMax*-1.297E-05)-2.816064E-04;
    if(channel == 629) return nlcblrwfSlope-( 2.171356E-05+trapMax*-1.295E-05)-1.197040E-04;
    if(channel == 640) return nlcblrwfSlope-( 7.818519E-05+trapMax*-1.296E-05)-2.513343E-04;
    if(channel == 641) return nlcblrwfSlope-( 3.668705E-05+trapMax*-1.296E-05)-1.063627E-04;
    if(channel == 642) return nlcblrwfSlope-( 5.146852E-05+trapMax*-1.373E-05)-1.911555E-04;
    if(channel == 643) return nlcblrwfSlope-( 2.370897E-05+trapMax*-1.375E-05)-8.660089E-05;
    if(channel == 644) return nlcblrwfSlope-( 6.066896E-05+trapMax*-1.308E-05)-2.776511E-04;
    if(channel == 645) return nlcblrwfSlope-(-1.188733E-05+trapMax*-1.305E-05)-1.153726E-04;
    if(channel == 646) return nlcblrwfSlope-( 7.972405E-05+trapMax*-1.291E-05)-2.211304E-04;
    if(channel == 647) return nlcblrwfSlope-( 2.540927E-05+trapMax*-1.286E-05)-8.877251E-04;
    if(channel == 656) return nlcblrwfSlope-( 1.262340E-05+trapMax*-1.337E-05)-1.702642E-04;
    if(channel == 657) return nlcblrwfSlope-( 1.086301E-05+trapMax*-1.342E-05)-8.161982E-05;
    if(channel == 662) return nlcblrwfSlope-( 1.485632E-04+trapMax*-1.460E-05)-5.350548E-04;
    if(channel == 663) return nlcblrwfSlope-( 2.742351E-06+trapMax*-1.457E-05)-1.100618E-04;
    if(channel == 664) return nlcblrwfSlope-( 5.110392E-05+trapMax*-1.321E-05)-1.847819E-04;
    if(channel == 665) return nlcblrwfSlope-( 6.499560E-06+trapMax*-1.319E-05)-9.712239E-05;
    if(channel == 674) return nlcblrwfSlope-( 9.983398E-05+trapMax*-1.310E-05)-6.408777E-04;
    if(channel == 675) return nlcblrwfSlope-( 3.757473E-05+trapMax*-1.309E-05)-2.113617E-04;
    if(channel == 688) return nlcblrwfSlope-( 2.461617E-05+trapMax*-1.287E-05)-4.178157E-04;
    if(channel == 689) return nlcblrwfSlope-( 3.018505E-05+trapMax*-1.287E-05)-1.380808E-04;
    if(channel == 690) return nlcblrwfSlope-( 9.944427E-05+trapMax*-1.300E-05)-3.561727E-04;
    if(channel == 691) return nlcblrwfSlope-(-4.737995E-06+trapMax*-1.300E-05)-1.279477E-04;
    if(channel == 692) return nlcblrwfSlope-( 8.146529E-05+trapMax*-1.354E-05)-2.265576E-04;
    if(channel == 693) return nlcblrwfSlope-( 1.839544E-05+trapMax*-1.356E-05)-1.004911E-04;
    if(channel == 696) return nlcblrwfSlope-( 2.093490E-05+trapMax*-1.308E-05)-2.745397E-04;
    if(channel == 697) return nlcblrwfSlope-( 8.938457E-06+trapMax*-1.310E-05)-1.215824E-04;
    if(trapMax == 0) return 0;
    return nlcblrwfSlope/trapMax;
  }

  if(dsNumber == 1) {
    //params from 11507-11592
    if(channel == 672) return nlcblrwfSlope-(2.833847E-05+trapMax*-1.328E-05)-1.815976E-04 ;
    if(channel == 673) return nlcblrwfSlope-(4.734834E-05+trapMax*-1.334E-05)-1.102766E-04 ;
    if(channel == 690) return nlcblrwfSlope-(3.690947E-05+trapMax*-1.308E-05)-1.861127E-04 ;
    if(channel == 691) return nlcblrwfSlope-(1.605341E-05+trapMax*-1.310E-05)-1.024712E-04 ;
    if(channel == 692) return nlcblrwfSlope-(5.681803E-05+trapMax*-1.309E-05)-2.062761E-04 ;
    if(channel == 693) return nlcblrwfSlope-(2.453714E-05+trapMax*-1.312E-05)-9.459848E-05 ;
    if(channel == 578) return nlcblrwfSlope-(4.417989E-05+trapMax*-1.306E-05)-2.471341E-04 ;
    if(channel == 579) return nlcblrwfSlope-(1.106159E-05+trapMax*-1.304E-05)-1.087355E-04 ;
    if(channel == 580) return nlcblrwfSlope-(2.577260E-05+trapMax*-1.291E-05)-4.108245E-04 ;
    if(channel == 581) return nlcblrwfSlope-(3.041374E-05+trapMax*-1.290E-05)-1.306187E-04 ;
    if(channel == 582) return nlcblrwfSlope-(5.396620E-05+trapMax*-1.304E-05)-2.544365E-04 ;
    if(channel == 583) return nlcblrwfSlope-(1.977025E-05+trapMax*-1.302E-05)-1.316353E-04 ;//average for all LG channels
    if(channel == 592) return nlcblrwfSlope-(6.803550E-05+trapMax*-1.297E-05)-4.974235E-04 ;
    if(channel == 593) return nlcblrwfSlope-(4.309706E-05+trapMax*-1.296E-05)-1.586067E-04 ;
    if(channel == 594) return nlcblrwfSlope-(4.643130E-05+trapMax*-1.352E-05)-1.914537E-04 ;
    if(channel == 595) return nlcblrwfSlope-(1.550955E-06+trapMax*-1.351E-05)-1.316353E-04 ;//average for all LG channels
    if(channel == 598) return nlcblrwfSlope-(2.897476E-05+trapMax*-1.298E-05)-2.391778E-04 ;
    if(channel == 599) return nlcblrwfSlope-(9.612110E-06+trapMax*-1.298E-05)-9.659031E-05 ;
    if(channel == 600) return nlcblrwfSlope-(5.842046E-05+trapMax*-1.287E-05)-2.038157E-04 ;
    if(channel == 601) return nlcblrwfSlope-(5.372889E-05+trapMax*-1.288E-05)-1.089046E-04 ;
    if(channel == 608) return nlcblrwfSlope-(4.693324E-05+trapMax*-1.284E-05)-4.175150E-04 ;
    if(channel == 609) return nlcblrwfSlope-(1.055481E-05+trapMax*-1.282E-05)-1.439053E-04 ;
    if(channel == 610) return nlcblrwfSlope-(1.216405E-05+trapMax*-1.298E-05)-2.325430E-04 ;
    if(channel == 611) return nlcblrwfSlope-(2.884038E-05+trapMax*-1.299E-05)-1.038957E-04 ;
    if(channel == 616) return nlcblrwfSlope-(8.295462E-05+trapMax*-1.316E-05)-2.999593E-04 ;
    if(channel == 617) return nlcblrwfSlope-(7.294466E-05+trapMax*-1.310E-05)-6.829416E-04 ;
    if(channel == 626) return nlcblrwfSlope-(1.860959E-06+trapMax*-1.265E-05)-1.895566E-04 ;
    if(channel == 627) return nlcblrwfSlope-(7.000642E-06+trapMax*-1.264E-05)-9.793286E-05 ;
    if(channel == 632) return nlcblrwfSlope-(9.473908E-05+trapMax*-1.307E-05)-3.559619E-04 ;
    if(channel == 633) return nlcblrwfSlope-(6.874599E-06+trapMax*-1.304E-05)-1.610833E-04 ;
    if(channel == 640) return nlcblrwfSlope-(8.088126E-05+trapMax*-1.287E-05)-2.227665E-04 ;
    if(channel == 641) return nlcblrwfSlope-(3.667507E-05+trapMax*-1.287E-05)-1.031861E-04 ;
    if(channel == 648) return nlcblrwfSlope-(8.317390E-05+trapMax*-1.311E-05)-4.997305E-04 ;
    if(channel == 649) return nlcblrwfSlope-(4.360668E-06+trapMax*-1.308E-05)-1.948269E-04 ;
    if(channel == 664) return nlcblrwfSlope-(5.916704E-05+trapMax*-1.321E-05)-2.549099E-04 ;
    if(channel == 665) return nlcblrwfSlope-(2.373893E-05+trapMax*-1.320E-05)-1.133329E-04 ;


    if(trapMax == 0) return 0;
    return nlcblrwfSlope/trapMax;
  }

  else //cout << "GetDCR98(): unknown dataset number DS" << dsNumber << endl;
  return 0;
}


double GetDCR99(int channel, double nlcblrwfSlope, double trapMax, int dsNumber)
{
  if(dsNumber == 0) {
    if(channel == 576) return nlcblrwfSlope-( 3.673385E-05+trapMax*-1.279E-05)-2.943105E-04;
    if(channel == 577) return nlcblrwfSlope-( 3.122161E-05+trapMax*-1.278E-05)-1.116883E-04;
    if(channel == 592) return nlcblrwfSlope-( 1.371969E-04+trapMax*-1.287E-05)-3.338696E-04;
    if(channel == 593) return nlcblrwfSlope-( 4.661555E-05+trapMax*-1.285E-05)-1.257624E-04;
    if(channel == 594) return nlcblrwfSlope-( 7.684589E-05+trapMax*-1.294E-05)-2.709597E-04;
    if(channel == 595) return nlcblrwfSlope-( 5.906154E-05+trapMax*-1.294E-05)-1.215498E-04;
    if(channel == 598) return nlcblrwfSlope-( 1.842168E-04+trapMax*-1.315E-05)-1.895059E-03;
    if(channel == 599) return nlcblrwfSlope-( 7.434566E-05+trapMax*-1.315E-05)-7.290344E-04;
    if(channel == 600) return nlcblrwfSlope-( 5.645349E-05+trapMax*-1.264E-05)-2.532124E-04;
    if(channel == 601) return nlcblrwfSlope-( 1.644013E-05+trapMax*-1.260E-05)-1.343034E-04;
    if(channel == 608) return nlcblrwfSlope-( 9.091576E-05+trapMax*-1.298E-05)-3.866514E-04;
    if(channel == 609) return nlcblrwfSlope-( 1.319922E-05+trapMax*-1.296E-05)-1.362231E-04;
    if(channel == 610) return nlcblrwfSlope-( 3.708218E-05+trapMax*-1.331E-05)-4.881074E-04;
    if(channel == 611) return nlcblrwfSlope-( 1.691204E-05+trapMax*-1.331E-05)-1.643671E-04;
    if(channel == 614) return nlcblrwfSlope-( 7.716150E-05+trapMax*-1.272E-05)-5.700044E-04;
    if(channel == 615) return nlcblrwfSlope-( 1.479029E-05+trapMax*-1.271E-05)-2.195198E-04;
    if(channel == 624) return nlcblrwfSlope-( 3.891078E-05+trapMax*-1.262E-05)-2.539550E-04;
    if(channel == 625) return nlcblrwfSlope-( 4.344515E-06+trapMax*-1.262E-05)-1.206225E-04;
    if(channel == 626) return nlcblrwfSlope-( 2.578407E-05+trapMax*-1.298E-05)-2.380339E-04;
    if(channel == 627) return nlcblrwfSlope-( 1.913080E-05+trapMax*-1.296E-05)-1.086605E-04;
    if(channel == 628) return nlcblrwfSlope-( 3.266857E-05+trapMax*-1.297E-05)-3.433920E-04;
    if(channel == 629) return nlcblrwfSlope-( 2.171356E-05+trapMax*-1.295E-05)-1.375756E-04;
    if(channel == 640) return nlcblrwfSlope-( 7.818519E-05+trapMax*-1.296E-05)-3.114656E-04;
    if(channel == 641) return nlcblrwfSlope-( 3.668705E-05+trapMax*-1.296E-05)-1.278493E-04;
    if(channel == 642) return nlcblrwfSlope-( 5.146852E-05+trapMax*-1.373E-05)-2.291227E-04;
    if(channel == 643) return nlcblrwfSlope-( 2.370897E-05+trapMax*-1.375E-05)-1.034860E-04;
    if(channel == 644) return nlcblrwfSlope-( 6.066896E-05+trapMax*-1.308E-05)-3.411170E-04;
    if(channel == 645) return nlcblrwfSlope-(-1.188733E-05+trapMax*-1.305E-05)-1.401627E-04;
    if(channel == 646) return nlcblrwfSlope-( 7.972405E-05+trapMax*-1.291E-05)-7.431032E-04;
    if(channel == 647) return nlcblrwfSlope-( 2.540927E-05+trapMax*-1.286E-05)-1.122732E-03;
    if(channel == 656) return nlcblrwfSlope-( 1.262340E-05+trapMax*-1.337E-05)-2.086850E-04;
    if(channel == 657) return nlcblrwfSlope-( 1.086301E-05+trapMax*-1.342E-05)-9.781184E-05;
    if(channel == 662) return nlcblrwfSlope-( 1.485632E-04+trapMax*-1.460E-05)-6.434840E-04;
    if(channel == 663) return nlcblrwfSlope-( 2.742351E-06+trapMax*-1.457E-05)-1.311618E-04;
    if(channel == 664) return nlcblrwfSlope-( 5.110392E-05+trapMax*-1.321E-05)-2.597014E-04;
    if(channel == 665) return nlcblrwfSlope-( 6.499560E-06+trapMax*-1.319E-05)-1.339893E-04;
    if(channel == 674) return nlcblrwfSlope-( 9.983398E-05+trapMax*-1.310E-05)-8.513833E-04;
    if(channel == 675) return nlcblrwfSlope-( 3.757473E-05+trapMax*-1.309E-05)-2.744186E-04;
    if(channel == 688) return nlcblrwfSlope-( 2.461617E-05+trapMax*-1.287E-05)-4.855323E-04;
    if(channel == 689) return nlcblrwfSlope-( 3.018505E-05+trapMax*-1.287E-05)-1.607455E-04;
    if(channel == 690) return nlcblrwfSlope-( 9.944427E-05+trapMax*-1.300E-05)-4.550779E-04;
    if(channel == 691) return nlcblrwfSlope-(-4.737995E-06+trapMax*-1.300E-05)-1.582763E-04;
    if(channel == 692) return nlcblrwfSlope-( 8.146529E-05+trapMax*-1.354E-05)-3.431355E-04;
    if(channel == 693) return nlcblrwfSlope-( 1.839544E-05+trapMax*-1.356E-05)-1.244729E-04;
    if(channel == 696) return nlcblrwfSlope-( 2.093490E-05+trapMax*-1.308E-05)-3.268483E-04;
    if(channel == 697) return nlcblrwfSlope-( 8.938457E-06+trapMax*-1.310E-05)-1.413939E-04;
    if(trapMax == 0) return 0;
    return nlcblrwfSlope/trapMax;
  }

  if(dsNumber == 1) {
    //params from 11507-11592
    if(channel == 672) return nlcblrwfSlope-(2.833847E-05+trapMax*-1.328E-05)-2.196966E-04;
    if(channel == 673) return nlcblrwfSlope-(4.734834E-05+trapMax*-1.334E-05)-1.301283E-04;
    if(channel == 690) return nlcblrwfSlope-(3.690947E-05+trapMax*-1.308E-05)-3.892826E-04;
    if(channel == 691) return nlcblrwfSlope-(1.605341E-05+trapMax*-1.310E-05)-1.785248E-04;
    if(channel == 692) return nlcblrwfSlope-(5.681803E-05+trapMax*-1.309E-05)-2.678930E-04;
    if(channel == 693) return nlcblrwfSlope-(2.453714E-05+trapMax*-1.312E-05)-1.176745E-04;
    if(channel == 578) return nlcblrwfSlope-(4.417989E-05+trapMax*-1.306E-05)-2.913426E-04;
    if(channel == 579) return nlcblrwfSlope-(1.106159E-05+trapMax*-1.304E-05)-1.342151E-04;
    if(channel == 580) return nlcblrwfSlope-(2.577260E-05+trapMax*-1.291E-05)-4.716845E-04;
    if(channel == 581) return nlcblrwfSlope-(3.041374E-05+trapMax*-1.290E-05)-1.497926E-04;
    if(channel == 582) return nlcblrwfSlope-(5.396620E-05+trapMax*-1.304E-05)-3.266592E-04;
    if(channel == 583) return nlcblrwfSlope-(1.977025E-05+trapMax*-1.302E-05)-1.795207E-04;//average value for all LG channels
    if(channel == 592) return nlcblrwfSlope-(6.803550E-05+trapMax*-1.297E-05)-7.344548E-04;
    if(channel == 593) return nlcblrwfSlope-(4.309706E-05+trapMax*-1.296E-05)-2.255520E-04;
    if(channel == 594) return nlcblrwfSlope-(4.643130E-05+trapMax*-1.352E-05)-2.328330E-04;
    if(channel == 595) return nlcblrwfSlope-(1.550955E-06+trapMax*-1.351E-05)-1.795207E-04;//average value for all LG channels
    if(channel == 598) return nlcblrwfSlope-(2.897476E-05+trapMax*-1.298E-05)-3.320390E-04;
    if(channel == 599) return nlcblrwfSlope-(9.612110E-06+trapMax*-1.298E-05)-1.228940E-04;
    if(channel == 600) return nlcblrwfSlope-(5.842046E-05+trapMax*-1.287E-05)-8.579336E-04;
    if(channel == 601) return nlcblrwfSlope-(5.372889E-05+trapMax*-1.288E-05)-2.916546E-04;
    if(channel == 608) return nlcblrwfSlope-(4.693324E-05+trapMax*-1.284E-05)-5.116411E-04;
    if(channel == 609) return nlcblrwfSlope-(1.055481E-05+trapMax*-1.282E-05)-1.752132E-04;
    if(channel == 610) return nlcblrwfSlope-(1.216405E-05+trapMax*-1.298E-05)-2.729985E-04;
    if(channel == 611) return nlcblrwfSlope-(2.884038E-05+trapMax*-1.299E-05)-1.225513E-04;
    if(channel == 616) return nlcblrwfSlope-(8.295462E-05+trapMax*-1.316E-05)-6.543178E-04;
    if(channel == 617) return nlcblrwfSlope-(7.294466E-05+trapMax*-1.310E-05)-1.052964E-03;
    if(channel == 626) return nlcblrwfSlope-(1.860959E-06+trapMax*-1.265E-05)-2.306604E-04;
    if(channel == 627) return nlcblrwfSlope-(7.000642E-06+trapMax*-1.264E-05)-1.191646E-04;
    if(channel == 632) return nlcblrwfSlope-(9.473908E-05+trapMax*-1.307E-05)-4.958992E-04;
    if(channel == 633) return nlcblrwfSlope-(6.874599E-06+trapMax*-1.304E-05)-1.981503E-04;
    if(channel == 640) return nlcblrwfSlope-(8.088126E-05+trapMax*-1.287E-05)-3.634485E-04;
    if(channel == 641) return nlcblrwfSlope-(3.667507E-05+trapMax*-1.287E-05)-1.330648E-04;
    if(channel == 648) return nlcblrwfSlope-(8.317390E-05+trapMax*-1.311E-05)-9.261722E-04;
    if(channel == 649) return nlcblrwfSlope-(4.360668E-06+trapMax*-1.308E-05)-2.763007E-04;
    if(channel == 664) return nlcblrwfSlope-(5.916704E-05+trapMax*-1.321E-05)-4.001758E-04;
    if(channel == 665) return nlcblrwfSlope-(2.373893E-05+trapMax*-1.320E-05)-1.605581E-04;

    if(trapMax == 0) return 0;
    return nlcblrwfSlope/trapMax;
  }

  else //cout << "GetDCR99(): unknown dataset number DS" << dsNumber << endl;
  return 0;
}
