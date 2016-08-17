void ds0_muGeSkim()
{
	TChain *skim = new TChain("skimTree");
	skim->Add("~/dev/datasets/ds0/*.root");
	long entries = skim->GetEntries();

	// initialize
	int run=0,event=0,mH=0,mL=0;
	double startTime=0,stopTime=0,sumEH=0,sumEL=0;
	vector<int> *channel=0,*gain=0,*str=0,*row=0,*mageID=0,*detID=0,*muType=0;
	vector<bool> *isEnr=0,*isGood=0,*badScaler=0,*muVeto1ms=0;
	vector<double> *trapECal=0,*aenorm=0,*toe=0,*trapETailMin=0,*tloc_s=0,*time_s=0,*mAct_g=0,*dtmu_s=0;
	if(1){
		skim->SetBranchAddress("run",&run);
		skim->SetBranchAddress("startTime",&startTime);
		skim->SetBranchAddress("stopTime",&stopTime);
		skim->SetBranchAddress("event",&event);
		skim->SetBranchAddress("trapECal",&trapECal);
		skim->SetBranchAddress("aenorm",&aenorm);
		skim->SetBranchAddress("toe",&toe);
		skim->SetBranchAddress("trapETailMin",&trapETailMin);
		skim->SetBranchAddress("sumEH",&sumEH);
		skim->SetBranchAddress("mH",&mH);
		skim->SetBranchAddress("sumEL",&sumEL);
		skim->SetBranchAddress("mL",&mL);
		skim->SetBranchAddress("tloc_s",&tloc_s);
		skim->SetBranchAddress("time_s",&time_s);
		skim->SetBranchAddress("channel",&channel);
		skim->SetBranchAddress("gain",&gain);
		skim->SetBranchAddress("str",&str);
		skim->SetBranchAddress("row",&row);
		skim->SetBranchAddress("mageID",&mageID);
		skim->SetBranchAddress("detID",&detID);
		// skim->SetBranchAddress("detName",&detName);
		skim->SetBranchAddress("isEnr",&isEnr);
		skim->SetBranchAddress("mAct_g",&mAct_g);
		skim->SetBranchAddress("isGood",&isGood);
		skim->SetBranchAddress("dtmu_s",&dtmu_s);
		skim->SetBranchAddress("muType",&muType);
		skim->SetBranchAddress("badScaler",&badScaler);
		skim->SetBranchAddress("muVeto1ms",&muVeto1ms);
	}

	TH1D *mult_all = new TH1D("mult_all","",21,0,21);
	TH1D *mult_veto = new TH1D("mult_veto","",21,0,21);

	TH1D *vetoSpec_hit = new TH1D("hit","",129,100,13000);	// 100 keV/bin
	TH1D *vetoSpec_sum = new TH1D("sum","",39,1000,40000);	// 1000 keV/bin

	int gEventsOver2650 = 0;	// total events over 2650
	int vEventsOver2650 = 0;	// only veto coincidences
	for (long i = 0; i < entries; i++)
	{
		skim->GetEntry(i);

		bool good = 0;
		double dtmu = dtmu_s->at(0);
		bool bads = badScaler->at(0);

		for (int j = 0; j < (int)channel->size(); j++)
		{
			if (isGood->at(j)) good=1;
		}

		if(good && sumEL > 1000) {
			mult_all->Fill(mL);

			if (dtmu > -0.2e-3 && dtmu < 1)
			{
				mult_veto->Fill(mL);
			}

			// what's that weird event?
			if (mH > 15){
				printf("run %i  event %i  dtmu_s %.2f  badScaler %i  mL %i  sumEL %.2f\n",run,event,dtmu,bads,mL,sumEL);
			}
		}

		if (dtmu > -0.2e-3 && dtmu < 1 && good)
		{
			if (mL !=0) {
				// printf("run %i  dtmu %.3f  mL %i  sumEL %.2f\n",run,dtmu,mL,sumEL);
				vetoSpec_sum->Fill(sumEL);
			}

			int chans = (int)channel->size();
			for (int j = 0; j < chans; j++)
			{
				if ((int)channel->at(j)%2==1)	// low gain channels only
				{
					vetoSpec_hit->Fill(trapECal->at(j));
				}
			}

			if (sumEL > 2650) vEventsOver2650++;
		}
		else if (!badScaler) cout << " bad scaler ! " << endl;

		if (sumEL > 2650 && good) gEventsOver2650++;
	}

	printf("Done with scan.\n");
	printf("Total events over 2650: %i  Veto-Coin events over 2650: %i\n",gEventsOver2650,vEventsOver2650);

	// ================== make some plots ====================

	TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",800,600);
	c1->SetLogy();
	mult_all->GetXaxis()->SetTitle("LG Multiplicity"); // note for caption: Sum-E > 1MeV events only.
	mult_all->GetYaxis()->SetTitle("Counts");
	mult_all->Draw();

	mult_veto->SetLineColor(kRed);
	mult_veto->Draw("same");

	TLegend* leg1 = new TLegend(0.5,0.6,0.87,0.85);
	leg1->AddEntry(mult_all,"All Events","l");
	leg1->AddEntry(mult_veto,"Veto-Coin Events","l");
	leg1->Draw();

	c1->Update();
	c1->Print("vetoMult_DS0_hit.pdf");

	TCanvas *c2 = new TCanvas("c2","Bob Ross's Canvas",800,600);
	c2->SetLogy();
	vetoSpec_hit->GetXaxis()->SetTitle("Hit Energy (keV)");
	vetoSpec_hit->GetXaxis()->SetTitleOffset(1.1);
	vetoSpec_hit->GetYaxis()->SetTitle("Counts (100 keV / bin)");
	// vetoSpec_hit->GetYaxis()->SetTitleOffset(1.1);
	vetoSpec_hit->Draw();
	c2->Update();
	c2->Print("vetoSpec_DS0_hit.pdf");

	TCanvas *c3 = new TCanvas("c3","Bob Ross's Canvas",800,600);
	vetoSpec_sum->GetXaxis()->SetTitle("Sum Energy (keV)");
	vetoSpec_sum->GetXaxis()->SetLabelSize(25);
	vetoSpec_sum->GetYaxis()->SetTitle("Counts (1 MeV / bin)");
	vetoSpec_sum->GetYaxis()->SetTitleOffset(0.9);
	vetoSpec_sum->Draw();

	TLegend* leg2 = new TLegend(0.3,0.6,0.87,0.87);
	leg2->AddEntry(vetoSpec_sum,"Veto-Coin Events > 1 MeV","");
	leg2->Draw();


	c3->Update();
	c3->Print("vetoSpec_DS0_sum.pdf");

}