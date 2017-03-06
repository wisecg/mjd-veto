#include <iostream>

#include "TFile.h"
#include "TROOT.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TGraph.h"

using namespace std;


int main()
{
	gROOT->ProcessLine(".x ~/env/MJDClintPlotStyle.C");
	// gROOT->ProcessLine( "gErrorIgnoreLevel = 2001;");

  TCanvas *can = new TCanvas("can","Bob Ross's Canvas",800,600);

  // DS-0
  // Chan 629 looks uncalibrated? Maybe need to exclude it.
  // TFile *f1 = new TFile("./data/skimDS1_HitsOver2630.root");
  // TTree *t1 = (TTree*)f1->Get("skimTree");
  // t1->Draw("dtmu_s","muVeto");
  // can->Print("DS1_MuonHits_dtmu.pdf");
  // TH1D *h1 = new TH1D("h1","h1",100,2500,13000);
  // h1->GetXaxis()->SetTitle("trapENFCal (hit)");
  // t1->Project("h1","trapENFCal","gain==1");
  // h1->Draw();
  // TH1D *h2 = new TH1D("h2","h2",100,2500,13000);
  // t1->Project("h2","trapENFCal","gain==1 && muVeto");
  // h2->SetLineColor(kRed);
  // h2->Draw("same");
  // can->Print("DS1_Over2630_hitSpec.pdf");
  // TAxis *ax1 = h1->GetXaxis();
  // TAxis *ax2 = h2->GetXaxis();
  // printf("DS1: Total hits > 2630: %.0f  Muon hits > 2630: %.0f  Hits > 7000: %.0f  Muon hits > 7000: %.0f\n", h1->Integral(ax1->FindBin(2630),ax2->FindBin(13000)), h2->Integral(ax1->FindBin(2630),ax2->FindBin(13000)), h1->Integral(ax1->FindBin(7000),ax2->FindBin(13000)), h2->Integral(ax1->FindBin(7000),ax2->FindBin(13000)));

  // DS-1
  TFile *f1 = new TFile("./data/skimDS1_HitsOver2630.root");
  TTree *t1 = (TTree*)f1->Get("skimTree");
  t1->Draw("dtmu_s","muVeto");
  can->Print("DS1_MuonHits_dtmu.pdf");
  TH1D *h1 = new TH1D("h1","h1",100,2500,13000);
  h1->GetXaxis()->SetTitle("trapENFCal (hit)");
  t1->Project("h1","trapENFCal","gain==1");
  h1->Draw();
  TH1D *h2 = new TH1D("h2","h2",100,2500,13000);
  t1->Project("h2","trapENFCal","gain==1 && muVeto");
  h2->SetLineColor(kRed);
  h2->Draw("same");
  can->Print("DS1_Over2630_hitSpec.pdf");
  TAxis *ax1 = h1->GetXaxis();
  TAxis *ax2 = h2->GetXaxis();
  printf("DS1: Total hits > 2630: %.0f  Muon hits > 2630: %.0f  Hits > 7000: %.0f  Muon hits > 7000: %.0f\n", h1->Integral(ax1->FindBin(2630),ax2->FindBin(13000)), h2->Integral(ax1->FindBin(2630),ax2->FindBin(13000)), h1->Integral(ax1->FindBin(7000),ax2->FindBin(13000)), h2->Integral(ax1->FindBin(7000),ax2->FindBin(13000)));

  /*
  // DS-2
  TFile *f2 = new TFile("skimDS3_HitsOver2630.root");
  TTree *t2 = (TTree*)f2->Get("skimTree");
  t2->Draw("dtmu_s","muVeto && dtmu_s < 0.002");
  can->Print("DS3_MuonHits_dtmu.pdf");
  TH1D *h3 = new TH1D("h3","h3",100,2500,13000);
  h3->GetXaxis()->SetTitle("trapENFDBSGCal (hit)");
  t2->Project("h3","trapENFDBSGCal","gain==1");
  h3->Draw();
  TH1D *h4 = new TH1D("h4","h4",100,2500,13000);
  t2->Project("h4","trapENFDBSGCal","gain==1 && muVeto");
  h4->SetLineColor(kRed);
  h4->Draw("same");
  can->Print("DS3_Over2630_hitSpec.pdf");
  TAxis *ax3 = h3->GetXaxis();
  TAxis *ax4 = h4->GetXaxis();
  printf("DS3: Total hits > 2630: %.0f  Muon hits > 2630: %.0f  Hits > 7000: %.0f  Muon hits > 7000: %.0f\n", h3->Integral(ax3->FindBin(2630),ax4->FindBin(13000)), h4->Integral(ax3->FindBin(2630),ax4->FindBin(13000)), h3->Integral(ax3->FindBin(7000),ax4->FindBin(13000)), h4->Integral(ax3->FindBin(7000),ax4->FindBin(13000)));


  // DS-3
  TFile *f3 = new TFile("skimDS3_HitsOver2630.root");
  TTree *t3 = (TTree*)f3->Get("skimTree");
  t3->Draw("dtmu_s","muVeto && dtmu_s < 0.002");
  can->Print("DS3_MuonHits_dtmu.pdf");
  TH1D *h3 = new TH1D("h3","h3",100,2500,13000);
  h3->GetXaxis()->SetTitle("trapENFDBSGCal (hit)");
  t3->Project("h3","trapENFDBSGCal","gain==1");
  h3->Draw();
  TH1D *h4 = new TH1D("h4","h4",100,2500,13000);
  t3->Project("h4","trapENFDBSGCal","gain==1 && muVeto");
  h4->SetLineColor(kRed);
  h4->Draw("same");
  can->Print("DS3_Over2630_hitSpec.pdf");
  TAxis *ax3 = h3->GetXaxis();
  TAxis *ax4 = h4->GetXaxis();
  printf("DS3: Total hits > 2630: %.0f  Muon hits > 2630: %.0f  Hits > 7000: %.0f  Muon hits > 7000: %.0f\n", h3->Integral(ax3->FindBin(2630),ax4->FindBin(13000)), h4->Integral(ax3->FindBin(2630),ax4->FindBin(13000)), h3->Integral(ax3->FindBin(7000),ax4->FindBin(13000)), h4->Integral(ax3->FindBin(7000),ax4->FindBin(13000)));

  // DS-4
  TFile *f3 = new TFile("skimDS4_HitsOver2630.root");
  TTree *t3 = (TTree*)f3->Get("skimTree");
  t3->Draw("dtmu_s","muVeto");
  can->Print("DS4_MuonHits_dtmu.pdf");
  TH1D *h5 = new TH1D("h5","h5",100,2500,13000);
  h5->GetXaxis()->SetTitle("trapENFCal (hit)");
  t3->Project("h5","trapENFCal","gain==1");
  h5->Draw();
  TH1D *h6 = new TH1D("h6","h6",100,2500,13000);
  t3->Project("h6","trapENFCal","gain==1 && muVeto");
  h6->SetLineColor(kRed);
  h6->Draw("same");
  can->Print("DS4_Over2630_hitSpec.pdf");
  TAxis *ax5 = h5->GetXaxis();
  TAxis *ax6 = h6->GetXaxis();
  printf("DS4: Total hits > 2630: %.0f  Muon hits > 2630: %.0f  Hits > 7000: %.0f  Muon hits > 7000: %.0f\n", h5->Integral(ax5->FindBin(2630),ax6->FindBin(13000)), h6->Integral(ax5->FindBin(2630),ax6->FindBin(13000)), h5->Integral(ax5->FindBin(7000),ax6->FindBin(13000)), h6->Integral(ax5->FindBin(7000),ax6->FindBin(13000)));

  // DS-5
  // NOTE: it looks like you need to check each file -- 644 is giving a ton of extra hits (it's not calibrated?) and you need to exclude it.  Isn't it a pulser monitor channel?
  TFile *f4 = new TFile("skimDS5_HitsOver2630.root");
  TTree *t4 = (TTree*)f4->Get("skimTree");
  t4->Draw("dtmu_s","muVeto");
  can->Print("DS5_MuonHits_dtmu.pdf");
  TH1D *h7 = new TH1D("h7","h7",100,2500,13000);
  h7->GetXaxis()->SetTitle("trapENFCal (hit)");
  t4->Project("h7","trapENFCal","gain==1");
  h7->Draw();
  TH1D *h8 = new TH1D("h8","h8",100,2500,13000);
  t4->Project("h8","trapENFCal","gain==1 && muVeto");
  h8->SetLineColor(kRed);
  h8->Draw("same");
  can->Print("DS5_Over2630_hitSpec.pdf");
  TAxis *ax7 = h7->GetXaxis();
  TAxis *ax8 = h8->GetXaxis();
  printf("DS5: Total hits > 2630: %.0f  Muon hits > 2630: %.0f  Hits > 7000: %.0f  Muon hits > 7000: %.0f\n", h7->Integral(ax7->FindBin(2630),ax8->FindBin(13000)), h8->Integral(ax7->FindBin(2630),ax8->FindBin(13000)), h7->Integral(ax7->FindBin(7000),ax8->FindBin(13000)), h8->Integral(ax7->FindBin(7000),ax8->FindBin(13000)));

  */
}
