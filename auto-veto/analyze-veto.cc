
// You don't need all these #includes

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
#include "TF1.h"

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

void CheckTimingResolution();

int main(int argc, const char** argv)
{
  CheckTimingResolution();
}

void CheckTimingResolution()
{
  int dsNum = 5;
  map<int,int> dsMap = {{0,76},{1,51},{3,24},{4,22},{5,80}};

  TH1D *hDelta = new TH1D("hDelta","hDelta",1000,-0.001,0.001);

  // for (int i = 0; i <= dsMap[dsNum]; i++) {
  for (int i = 0; i <= 0; i++)
  {
    cout << "Loading DS-" << dsNum << " run sequence " << i << endl;
    GATDataSet ds;
    LoadDataSet(ds, dsNum, i);
    TChain *vetoChain = ds.GetVetoChain();
    cout << "Found " << vetoChain->GetEntries() << endl;

    TTreeReader vetoReader(vetoChain);
    TTreeReaderValue<double> scaler(vetoReader,"fTimeSec");
    TTreeReaderValue<double> sbc(vetoReader,"fTimeSBC");
    TTreeReaderValue<bool> badScaler(vetoReader,"fBadScaler");

    double prevScaler=0, prevSBC=0;
    bool prevBadScaler = true;
  	while(vetoReader.Next())
  	{
      // Calculate deltas and histogram them
      if (!(*badScaler) && !prevBadScaler)
      {
        double deltaScaler = *scaler - prevScaler;
        double deltaSBC = *sbc - prevSBC;
        double delta = deltaScaler - deltaSBC;
        hDelta->Fill(delta);
      }

      prevScaler=*scaler;
      prevSBC=*sbc;
      prevBadScaler=*badScaler;
    }
  }

  // gROOT->ProcessLine(".x ~/env/MJDClintPlotStyle.C");
  TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",800,600);
  hDelta->SetTitle("");
  hDelta->Draw();
  hDelta->Fit("gaus","q");
  double mean = hDelta->GetFunction("gaus")->GetParameter(1);
  double sigma = hDelta->GetFunction("gaus")->GetParameter(2);

  cout << "mean: " << mean << "  sigma: " << sigma << endl;
  c1->Print("./output/deltaScalerSBC.pdf");


}