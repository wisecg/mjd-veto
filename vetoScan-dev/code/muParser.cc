#include "vetoScan.hh"
using namespace std;
void muParser(string arg) 
{
	/*
	TChain *v = new TChain("vetoEvent");
	
	// Can handle more than one input file
	// v->AddFile(arg.c_str());
	v->AddFile("./output/KJR_Total.root");
	// v->AddFile("./output/P3KJR_R01.root");
	// v->AddFile("./output/P3KJR_R02.root");

	MJVetoEvent *event = NULL;
	int CoinType[32] = {0};
	int CutType[32] = {0};
	double LEDfreq = 0;
	double LEDrms = 0;
	int highestMultip = 0;
	Long64_t start = 0;
	Long64_t stop = 0;
	v->SetBranchAddress("events",&event);
	v->SetBranchAddress("CoinType[32]",CoinType);
	v->SetBranchAddress("CutType[32]",CutType);
	v->SetBranchAddress("LEDfreq",&LEDfreq);
	v->SetBranchAddress("LEDrms",&LEDrms);
	v->SetBranchAddress("highestMultip",&highestMultip);
	v->SetBranchAddress("start",&start);
	v->SetBranchAddress("stop",&stop);
	int vEntries = v->GetEntries();
	cout << "Found " << vEntries << " entries.\n";

	v->GetEntry(0);
	long firstStart = start;
	v->GetEntry(vEntries-1);
	long lastStop = stop;
	double days = (double)(lastStop - firstStart)/(86400);
	double binsPerDay = 1;  //1,0.5,0.2
	int bins = (int)(binsPerDay*days);
	TH1D *muRate = new TH1D("muRate","",bins,firstStart,lastStop);
	TH1D *muRate_w = new TH1D("muRate_w","",bins,firstStart,lastStop);

	double muonTime = 0;
	int muonCount = 0;
	long liveTime = 0;
	int prevRun = 0;
	struct tm * ptm;
	
	// this should be tied to binsPerDay?
	int day = 0;
	int prevDay = 0;
	int muonsThisDay = 0;
	double liveTimeDay = 0;
	double liveFracDay = 0;

	for (long i = 0; i < vEntries; i++) 
	{
		v->GetEntry(i);

		// should not need a multiplicity cut ...
		// P3JDY really drove the number down too.
		if ((CoinType[1]||CoinType[2]) && event->GetMultip() < 24 && !event->GetBadScaler()) 
		{
			muonCount++;
			muonTime = start + event->GetTimeSec();
			muRate->Fill(muonTime);			

			// find the day
			time_t muonStamp = (time_t)muonTime;
			ptm = gmtime ( &muonStamp );
			day = ptm->tm_mday;


			if (day!=prevDay)
			{
				// if (event->run != prevRun) ;;;;
				
				liveFracDay = liveTimeDay/86400;
				muRate_w->Fill(muonsThisDay,liveFracDay);
				
				muonsThisDay = 0;
				liveFracDay = 0;
			}
			else {
				muonsThisDay++;
			}

			prevDay = day;
		}

		// live time
		// still assuming the veto was completely live during the run
		if (event->GetRun() != prevRun) liveTime += stop-start;

		prevRun = event->GetRun();
	}

	// don't forget to fill the last day of muonsThisDay.


	// output and plotting
	//
	MJDB::MJSlowControlsDoc doc;
	cout << doc.GetGMTString(firstStart) << "  " << doc.GetGMTString(lastStop) << endl;
	string startDate = doc.GetGMTString(firstStart);
	string stopDate = doc.GetGMTString(lastStop);
	double liveDays = (double)liveTime/86400;

	cout << "\nDates covered: " << startDate << " to " << stopDate << " (GMT)" << endl;
	printf("Total days : %.2f,  Total Live Time : %.2f.  Live-time fraction : %.2f%%.\n",
		days,liveDays,100*(liveDays/days));
	printf("Muon rate : Of %i entries, found %i muons in %.2f live days .....\n\t%.2f muons/day , %.2f muons/4hrs\n\n"
		,vEntries,muonCount,liveDays,muonCount/liveDays,muonCount/(liveDays*6));

	TCanvas *c = new TCanvas("c","Bob Ross's Canvas",800,600);
	muRate->GetXaxis()->SetTimeOffset(0,"gmt");
	muRate->GetXaxis()->SetTimeFormat("%m/%d/%y");
  	muRate->GetXaxis()->SetTimeDisplay(1);
  	muRate->GetXaxis()->SetNdivisions(-506);	//-503
  	muRate->GetXaxis()->SetLabelOffset(0.02);
  	muRate->GetXaxis()->SetTitleOffset(1.5);

  	char xTitle[200];
  	sprintf(xTitle,"[%.2f days / bin]",1/binsPerDay);
  	muRate->GetXaxis()->SetTitle(xTitle);

  	char yTitle[200];
  	sprintf(yTitle,"Muons / %.1f days",1/binsPerDay);
  	muRate->GetYaxis()->SetTitle(yTitle);

  	char Title[200];
  	sprintf(Title,"Muons since %s (%.1f days)",doc.GetGMTString(firstStart).c_str(),days);
  	muRate->SetTitle(Title);
  	muRate->SetStats(0);
	muRate->Draw();
	
	c->Print("./output/muRate.pdf");
	c->Print("./output/muRate.C");
	*/	
}