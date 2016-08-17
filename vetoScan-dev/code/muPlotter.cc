#include "vetoScan.hh"
using namespace std;

void muPlotter(string arg)
{
	//This can also handle more than 1 run
	// TChain *t = new TChain("vetoEvent");
	// t->AddFile(arg.c_str());

	TFile *f = new TFile(arg.c_str());
	TTree *t = (TTree*)f->Get("vetoEvent");
	
	// need to get the start/stop times for a big histo

	// TH1D *hqt = new TH1D("hqt","QDCTotal",100,0,55000);
	// TH1D *hqt2 = new TH1D("hqt2","QDCTotal+cut",100,0,55000);
	TH2D *hqm = new TH2D("hqm","qdctotal vs multip.",100,0,55000,32,0,32);

	char cut1[500];
	sprintf(cut1,"CoinType[1]==1");
	TCut tcut1 = cut1;

	// t->Project("hqt","events.totE");
	// t->Project("hqt2","events.totE",tcut1);
	t->Project("hqm","events.multip:events.totE",tcut1);

	TCanvas* c1 = new TCanvas("c1","Bob Ross's Canvas",800,600);
     
	c1->cd(1);
	// c1->SetLogy();
	// hqt2->SetLineColor(kRed);
	// hqt2->SetLineWidth(2.0);
	// hqt2->Draw();
	// hqt->Draw("same");

	hqm->Draw("COL");
	hqm->SetStats(0);

	c1->Print("./output/plotter.pdf");

	delete t;
	f->Close();
     
}