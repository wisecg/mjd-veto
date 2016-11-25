#include <iostream>
#include <string>
#include <map>
#include "TTreeReader.h"
#include "TROOT.h"
#include "TH1.h"
#include "TCanvas.h"
#include "GATDataSet.hh"
#include "MJTRun.hh"
#include "MGTBasicEvent.hh"
#include "MJVetoEvent.hh"
#include "MGTEvent.hh"
#include "MGVDigitizerData.hh"

#include "MJTChannelMap.hh"
#include "MJTChannelSettings.hh"

using namespace std;

void durationCheckGAT()
{
  GATDataSet ds(19070);
  TChain *c = ds.GetGatifiedChain();
  TTreeReader reader(c);
  TTreeReaderValue<double> startTime(reader, "startTime");
  TTreeReaderValue<double> stopTime(reader, "stopTime");
  TTreeReaderValue<vector<double>> timestamp(reader,"timestamp");
  reader.SetEntry(0);
  double start = *startTime;
  double stop = *stopTime;
  // duration = (double)(stop - start);
  reader.SetTree(c);  // resets the reader
  cout << "Start: " << start << " Stop: " << stop << "  diff: " << stop-start << endl;
}

void durationCheckBLT()
{
  TH1D *dt = new TH1D("dt","dt",100,-50,50);
  TH1D *dtPackets = new TH1D("dtPackets","dtPackets",100,-50,50);
  vector<int> runs;
  int rundummy;
  ifstream runList("./runs/p3lqk-complete.txt");
  while(runList >> rundummy) runs.push_back(rundummy);
  runList.close();
  long prevStopUnix=0;
  for (auto run : runs)
  {
    GATDataSet ds;
  	string runPath = ds.GetPathToRun(run,GATDataSet::kBuilt);
    // cout << runPath << endl;
    if (runPath == "") {
      cout << "File doesn't exist.  Continuing ..." << endl;
      continue;
    }
  	TChain *blt = new TChain("MGTree");
    if (!blt->Add(runPath.c_str())){
  		cout << "File exists, but chain is messed up.  Continuing ...\n";
  		continue;
  	}
    if (blt->GetEntries() == 0){
      cout << "No entries for run " << run << ".  Continuing ...\n";
      continue;
    }
    TTreeReader reader(blt);
    // TTreeReaderValue<MJTRun> bRun(reader,"run");
    TTreeReaderValue<double> bTime(reader,"fTime");
    TTreeReaderValue<long> bTimeStart(reader,"fStartTime");
    TTreeReaderValue<long> bTimeStop(reader,"fStopTime");
    reader.SetEntry(0);

    // does the method above of getting these always work?
    long startUnix = (*bTimeStart);
    long stopUnix = (*bTimeStop);

    // ian says these are worthless until the data is reprocessed
    // double startTS = (double)(bRun->GetStartTimeStamp())*1.e-8;
    // double stopTS = (double)(bRun->GetStopTimeStamp())*1.e-8;

    double firstTime = (*bTime)*1.e-9;
    reader.SetEntry(blt->GetEntries()-1);
    double lastTime = (*bTime)*1.e-9;

    double diff = (lastTime-firstTime) - (double)(stopUnix - startUnix);
    long diffPackets = startUnix-prevStopUnix;
    string whichIsAhead = "";
    if (diff > 0) whichIsAhead = "dig";
    else if (diff < 0) whichIsAhead = "cpu";
    dt->Fill(diff);
    dtPackets->Fill(diffPackets);

    printf("%i  %li  %li  %-5li  %-10.2f  %-10.2f  %-10.2f  %-8.1f  %s  %-8li  %.0f\n", run,startUnix,stopUnix,stopUnix-startUnix,firstTime,lastTime,lastTime-firstTime,diff,whichIsAhead.c_str(),diffPackets,diff-(double)diffPackets);
    delete blt;
    prevStopUnix = stopUnix;
  }
  TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",800,600);
  c1->SetLogy();
  dt->GetXaxis()->SetTitle("gretina duration - unix duration");
  dt->Draw();
  c1->Print("gretUnixDuration.pdf");
  dtPackets->GetXaxis()->SetTitle("start packet - prev stop packet");
  dtPackets->Draw();
  c1->Print("diffPackets.pdf");
}

void durationCheckVETO()
{
  vector<int> runs;
  int rundummy;
  ifstream runList("./runs/p3jdy-complete.txt");
  while(runList >> rundummy) runs.push_back(rundummy);
  runList.close();
  for (auto run : runs)
  {
    GATDataSet ds;
  	string runPath = ds.GetPathToRun(run,GATDataSet::kBuilt);
    // cout << runPath << endl;
    if (runPath == "") {
      cout << "File doesn't exist.  Continuing ..." << endl;
      continue;
    }
  	TChain *vetoChain = new TChain("VetoTree");
    if (!vetoChain->Add(runPath.c_str())){
  		cout << "File exists, but chain is messed up.  Continuing ...\n";
  		continue;
  	}
    if (vetoChain->GetEntries() == 0){
      cout << "No entries for run " << run << ".  Continuing ...\n";
      continue;
    }
    TTreeReader reader(vetoChain);
    // TTreeReaderValue<MJTRun> bRun(reader,"run");
    TTreeReaderValue<double> bTime(reader,"fTime");
    TTreeReaderValue<long> bTimeStart(reader,"fStartTime");
    TTreeReaderValue<long> bTimeStop(reader,"fStopTime");
    reader.SetEntry(0);

    // does the method above of getting these always work?
    long startUnix = (*bTimeStart);
    long stopUnix = (*bTimeStop);
    // works for p3lqk!
    // works for p3kjr!
    // works for p3jdy!

    printf("%i  %li  %li  %li\n",run,startUnix,stopUnix,stopUnix-startUnix);
  }
}

void ds3skimCheck()
{
  TChain *skim = new TChain("skimTree");
  skim->Add("/project/projectdirs/majorana/data/mjd/surfmjd/analysis/skim/DS3/gatrev_153453544/*.root");
  vector<int> runList;
  int cts = skim->Draw("run:trapENFDBSGCal:channel","trapENFDBSGCal > 6000 && channel%2==1","GOFF");
  for (int i = 0; i < cts; i++){
    int run = skim->GetV1()[i];
    double ene = skim->GetV2()[i];
    int chn = skim->GetV3()[i];
    printf("%i  %-7.1f  %i\n",run,ene,chn);
    runList.push_back(run);
  }
  for (auto i : runList) cout << i << endl;
  // 17417
  // 17460
  // 17460
  // 17467
  // 17604
  // 17604
  // 17608
  // 17608
  // 17608
  // 17675
  // 17675
  // 17675
  // 17677
  // 16991
  // 16991
  // 16991
  // 16991
  // 17109
}

void clockResetCheck()
{
  // 16797 - 1st DS3 run
  int runNum = 17324;

  GATDataSet ds(17324);

  TChain *vetoChain = new TChain("vetoTree");
  if (!vetoChain->Add(TString::Format("./avout/veto_run%i.root",runNum))){
    cout << "File doesn't exist.  Exiting ...\n";
    return;
  }
  TTreeReader reader(vetoChain);
	TTreeReaderValue<MJVetoEvent> vetoEventIn(reader,"vetoEvent");
  reader.Next();
  MJVetoEvent veto = *vetoEventIn;
  cout << "First Veto time: " << veto.GetTimeSec() << ", packet " << veto.GetScalerIndex() << endl;
  MJVetoEvent first = veto;

  TChain *gat = ds.GetGatifiedChain(false);
  vector<double> *timestamp = 0;
  gat->SetBranchAddress("timestamp",&timestamp);
  gat->GetEntry(0);
  cout << "First Ge time: " << timestamp->at(0)*1.e-8 << endl;

  TChain *builtChain = ds.GetBuiltChain(false);
  MGTEvent *evt=0;
  MGVDigitizerData *dig=0;
  builtChain->SetBranchAddress("event",&evt);
  double bTimeFirst=0, bTimeBefore=0, bTimeAfter=0;
  uint64_t bItr=0,bIndex=0;
  while (true)
  {
    builtChain->GetEntry(bItr);
    // CAUTION: Are you sure the built data isn't flushing the buffer at the beginning?
    // Does it matter if it is flushing?
    if (evt->GetNDigitizerData() > 0)
    {
      dig = evt->GetDigitizerData(0);
      bIndex = dig->GetIndex();
      if (bTimeFirst == 0)
        bTimeFirst = ((double)dig->GetTimeStamp())*1.e-8;
      if ((int)bIndex < first.GetScalerIndex())
        bTimeBefore = ((double)dig->GetTimeStamp())*1.e-8;
      bTimeAfter = ((double)dig->GetTimeStamp())*1.e-8;
      // printf("%li  ind %lu  first %.3f  before %.3f  after %.3f\n", bItr,bIndex,bTimeFirst,bTimeBefore,bTimeAfter);
    }
    if ((int)bIndex > first.GetScalerIndex()) break;
    if ((int)bIndex > builtChain->GetEntries()) break;
    bItr++;
  }
  double bVetoTime = (bTimeAfter + bTimeBefore)/2.;
  double scalerOffset = bVetoTime - bTimeFirst;
  double scalerUnc = (bTimeAfter - bTimeBefore)/2.;

  printf("First veto time from Ge timestamps: %.3f +/- %.3f sec.  Time after start of run: %.1f\n",bVetoTime,scalerUnc,scalerOffset);

  cout << "Time from veto - time from ge: " << veto.GetTimeSec() - bVetoTime << endl;

  // now get the thresholds for all the enabled channels.

  // Get the thresholds for all the enabled channels in this run.

  map<int,int> threshMap;
  MJTChannelMap *chanMap = ds.GetChannelMap();
  MJTChannelSettings *chanSet = ds.GetChannelSettings();
  vector<uint32_t> enabChans = chanSet->GetEnabledIDList();
  for (auto i : enabChans)
  {
    uint32_t crate = chanMap->GetInt(i,"kVME");
    uint32_t card = chanMap->GetInt(i,"kCardSlot");
    uint32_t sigChan = 0;
    if (i%2==0) sigChan = chanMap->GetInt(i,"kChanHi");
    else if (i%2==0) sigChan = chanMap->GetInt(i,"kChanLo");
    uint32_t thresh = chanSet->GetInt("TRAP Threshold",crate,card,sigChan,"ORGretina4MModel");
    threshMap.insert({i,thresh});
    // printf("%i:  Crate %i  Card %i  Sig.Chan %i  Thresh %i\n",i,crate,card,sigChan,thresh);
  }

}

int main()
{
  // gROOT->ProcessLine(".x ~/env/MJDClintPlotStyle.C");
  // durationCheckGAT();
  // durationCheckBLT();
  // durationCheckVETO();
  // ds3skimCheck();
  clockResetCheck();

}