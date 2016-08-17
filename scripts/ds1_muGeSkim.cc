// Calculate veto-germanium coincidence rate
// for DS1 skim files.
// Clint Wiseman, USC
// 5/29/16

void ds1_muGeSkim()
{
	TChain *skim = new TChain("skimTree");
	skim->Add("~/dev/datasets/ds1/*.root");
	long entries = skim->GetEntries();
	printf("Scanning skim files ... \nFound %li entries.\n",entries);

	// skim file branches.
	// May 2016 version.
	unsigned int gatrev=0,EventDC1Bits=0;
	int run=0,iEvent=0,mH=0,mL=0;
	double startTime=0,stopTime=0,sumEH=0,sumEL=0;
	vector<int> *iHit=0,*channel=0,*P=0,*D=0,*gain=0,*mageID=0,*detID=0,*dateMT=0,*muType=0;
	vector<double> *mAct_g=0,*tloc_s=0,*time_s=0,*timeMT=0,*trapECal=0,*trapENFCal=0,*aenorm=0,*kvorrT=0,*toe=0,*dcrSlope95=0,*dcrSlope99=0,*trapETailMin=0,*dtmu_s=0,*trap4usMax=0,*t150=0;
	vector<bool> *isEnr=0,*isNat=0,*isGood=0,*isLNFill=0,*badScaler=0,*muVeto1ms=0;
	// vector<string> *detName;
	if (1) {
		vector<unsigned int> *wfDCBits=0;
		skim->SetBranchAddress("gatrev",&gatrev);
		skim->SetBranchAddress("EventDC1Bits",&EventDC1Bits);
		skim->SetBranchAddress("run",&run);
		skim->SetBranchAddress("iEvent",&iEvent);
		skim->SetBranchAddress("mH",&mH);
		skim->SetBranchAddress("mL",&mL);
		skim->SetBranchAddress("t150",&t150);
		skim->SetBranchAddress("startTime",&startTime);
		skim->SetBranchAddress("stopTime",&stopTime);
		skim->SetBranchAddress("sumEH",&sumEH);
		skim->SetBranchAddress("sumEL",&sumEL);
		skim->SetBranchAddress("iHit",&iHit);
		skim->SetBranchAddress("channel",&channel);
		skim->SetBranchAddress("P",&P);
		skim->SetBranchAddress("D",&D);
		skim->SetBranchAddress("gain",&gain);
		skim->SetBranchAddress("mageID",&mageID);
		skim->SetBranchAddress("detID",&detID);
		skim->SetBranchAddress("dateMT",&dateMT);
		skim->SetBranchAddress("muType",&muType);
		skim->SetBranchAddress("mAct_g",&mAct_g);
		skim->SetBranchAddress("tloc_s",&tloc_s);
		skim->SetBranchAddress("time_s",&time_s);
		skim->SetBranchAddress("timeMT",&timeMT);
		skim->SetBranchAddress("trapECal",&trapECal);
		skim->SetBranchAddress("trapENFCal",&trapENFCal);
		skim->SetBranchAddress("trap4usMax",&trap4usMax);
		skim->SetBranchAddress("aenorm",&aenorm);
		skim->SetBranchAddress("kvorrT",&kvorrT);  // in GAT: trirt100nsft10nsMax / trapECal 
		skim->SetBranchAddress("toe",&toe);
		skim->SetBranchAddress("dcrSlope95",&dcrSlope95);
		skim->SetBranchAddress("dcrSlope99",&dcrSlope99);
		skim->SetBranchAddress("trapETailMin",&trapETailMin);	// events with positive trapETailMin are removed
		skim->SetBranchAddress("dtmu_s",&dtmu_s);
		skim->SetBranchAddress("isEnr",&isEnr);
		skim->SetBranchAddress("isNat",&isNat);
		skim->SetBranchAddress("isGood",&isGood);
		skim->SetBranchAddress("isLNFill",&isLNFill);
		skim->SetBranchAddress("badScaler",&badScaler);
		skim->SetBranchAddress("muVeto1ms",&muVeto1ms);
	}

	TH1D *vetoSpec_hit = new TH1D("h0","",129,100,13000);	// 100 keV/bin
	TH1D *vetoSpec_sum = new TH1D("h0","",39,1000,40000);	// 1000 keV/bin

	int gEventsOver2650 = 0;	// total events over 2650
	int vEventsOver2650 = 0;	// only veto coincidences
	for (long i = 0; i < entries; i++)
	{
		skim->GetEntry(i);
		double dtmu = dtmu_s->at(0);

		if (dtmu > -0.2e-3 && dtmu < 1 && EventDC1Bits==0)
		{
			if (mL !=0) {
				printf("run %i  dtmu %.3f  mL %i  sumEL %.2f\n",run,dtmu,mL,sumEL);
				vetoSpec_sum->Fill(sumEL);
			}

			int chans = (int)channel->size();
			for (int j = 0; j < chans; j++)
			{
				if ((int)channel->at(j)%2==1)	// low gain channels only
				{
					vetoSpec_hit->Fill(trapENFCal->at(j));
				}
			}

			if (sumEL > 2650) vEventsOver2650++;
		}
		else if (!badScaler) cout << " bad scaler ! " << endl;

		if (sumEL > 2650 && EventDC1Bits==0) gEventsOver2650++;
	}

	printf("Done with scan.\n");
	printf("Total events over 2650: %i  Veto-Coin events over 2650: %i\n",gEventsOver2650,vEventsOver2650);

	// ============= make some plots ==============

	TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",800,600);
	c1->SetLogy();
	vetoSpec_hit->GetXaxis()->SetTitle("Hit Energy (keV)");
	vetoSpec_hit->GetXaxis()->SetTitleOffset(1.1);
	vetoSpec_hit->GetYaxis()->SetTitle("Counts (100 keV / bin)");
	vetoSpec_hit->GetYaxis()->SetTitleOffset(0.9);
	vetoSpec_hit->Draw();
	c1->Update();
	c1->Print("vetoSpec_DS1_hit.pdf");

	TCanvas *c2 = new TCanvas("c2","Bob Ross's Canvas",800,600);
	vetoSpec_sum->GetXaxis()->SetTitle("Sum Energy (keV)");
	vetoSpec_sum->GetXaxis()->SetLabelSize(25);
	vetoSpec_sum->GetYaxis()->SetTitle("Counts (1 MeV / bin)");
	vetoSpec_sum->Draw();
	c2->Update();
	c2->Print("vetoSpec_DS1_sum.pdf");


}