// MJD Data Set Livetime Calculator.
// Future: Add corrections from LN fill, muon veto cut, and others.
//
// Clint Wiseman, USC/Majorana
// 10/6/2016

#include <iostream>
#include <map>
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "GATDataSet.hh"
#include "DataSetInfo.hh"
#include "MJVetoEvent.hh"

using namespace std;

double vetoReduction(GATDataSet &ds, int dsNum);

int main(int argc, char** argv)
{
	if (argc < 2) {
		cout << "Usage: ./ds_livetime [dataset number]\n";
		return 1;
	}
	int dsNum = stoi(argv[1]);

  // this must match the run sequences in DataSetInfo.hh
  map<int,int> dsMap = {{0,76},{1,51},{2,7},{3,24},{4,22},{5,80}};
  GATDataSet ds;
  for (int i = 0; i < dsMap[dsNum]; i++) {
    LoadDataSet(ds, dsNum, i);
  }
  double vetoDeadTime = vetoReduction(ds, dsNum);

  // NOTE: This assumes the veto'd time is within the run boundary.
  // DS-0 3159.18
  // Total dead time: 3159.18 seconds from 1623 muon candidate events.
  // DS-1 3323.95
  // Total dead time: 3323.95 seconds from 3023 muon candidate events.
  // DS-2 572.603
  // Total dead time: 572.603 seconds from 532 muon candidate events.
  // DS-3 1093.16
  // Total dead time: 1093.16 seconds from 1067 muon candidate events.
  // DS-4 9986.62
  // Total dead time: 9986.62 seconds from 1034 muon candidate events.
  // DS-5 FIXME
  // Total dead time: -9.07416e+06 seconds from 2607 muon candidate events.

  // double totalLiveTime = ds.GetRunTime()/1e9; // this can take a while
  // double adjLiveTime = totalLiveTime - vetoDeadTime;
  // cout << Form("DS-%i  Total %.2f sec, Veto Dead %.2f sec.  Adjusted: %.2f sec\n",dsNum,totalLiveTime,vetoDeadTime,adjLiveTime);
}

double vetoReduction(GATDataSet &ds, int dsNum)
{
  TChain *vetoChain = NULL;

  // Make the muon list, exactly the same way as we do in skim_mjd_data
  if(vetoChain==NULL && dsNum!=4) {
    vetoChain = ds.GetVetoChain();
    cout << "Found " << vetoChain->GetEntries() << " veto entries.  Creating muon list ...\n";
  }
  vector<int> muRuns;
  vector<int> muTypes;
  vector<double> muRunTStarts;
  vector<double> muTimes;
  vector<double> muUncert;
  if (dsNum != 4)
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
  else if (dsNum==4) LoadDS4MuonList(muRuns,muRunTStarts,muTimes,muTypes,muUncert);
  size_t nMu = muTimes.size();
  if(nMu == 0) {
    cout << "couldn't load mu data" << endl;
    return 0;
  }
  cout << "Muon list has " << muRuns.size() << " entries.\n";

  double deadTimeTotal = 0;
  for (int i = 0; i < (int)muRuns.size(); i++)
  {
    // this matches what's in skim_mjd_data
    double deadTime = 1. + 2 * fabs(muUncert[i]); // 1 second window, uncertainties on both ends.
    if (dsNum==4) deadTime = 4. + 4. * fabs(muUncert[i]);  // larger ~10s window for ds-4
    deadTimeTotal += deadTime;
    // printf("%i  %i  %i  %.0f  %.3f +/- %.3f  dead %.3f\n",i,muRuns[i],muTypes[i],muRunTStarts[i],muTimes[i],muUncert[i],deadTimeTotal);
    // FIXME: ds-5 muon uncertainties are being reported as negative
  }
  cout << "Total dead time: " << deadTimeTotal << " seconds from " << muRuns.size() << " muon candidate events.\n";

  return deadTimeTotal;
}