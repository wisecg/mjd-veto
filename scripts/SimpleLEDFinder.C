// ==================================================
// Processing Functions
// ==================================================
long GetStartUnixTime(GATDataSet ds)
{
	TChain *c = ds.GetVetoChain();
	MJTRun *runInfo = new MJTRun();
	c->SetBranchAddress("run",&runInfo);
	c->GetEntry(0);
	return (long)runInfo->GetStartTime();
}

long GetStopUnixTime(GATDataSet ds)
{
	TChain *c = ds.GetVetoChain();
	MJTRun *runInfo = new MJTRun();
	c->SetBranchAddress("run",&runInfo);
	c->GetEntry(0);
	return (long)runInfo->GetStopTime();
}

// ==================================================
// Analysis
// ==================================================
// Check if the veto data is getting high multiplicity events.
void SimpleLEDFinder(string file) 
{
	int run;
	long start;
	long stop;
	int duration;
	ifstream InputList;
  	InputList.open(file.c_str());
	while(!InputList.eof())
	{
		InputList >> run;
		GATDataSet ds(run);
		start = GetStartUnixTime(ds);
		stop = GetStopUnixTime(ds);
		duration = stop - start;
			
		// standard veto initialization block
		TChain *v = ds.GetVetoChain();
		long vEntries = v->GetEntries();
		MJTRun *vRun = new MJTRun();
		MGTBasicEvent *vEvent = new MGTBasicEvent(); 
		unsigned int mVeto = 0;
		uint32_t vBits = 0;
		v->SetBranchAddress("run",&vRun);
		v->SetBranchAddress("mVeto",&mVeto);
		v->SetBranchAddress("vetoEvent",&vEvent);
		v->SetBranchAddress("vetoBits",&vBits);
		v->GetEntry(0);
		MJTVetoData *vData[32];
		for (int i=0; i<32; i++) { vData[i] = dynamic_cast<MJTVetoData*>(vEvent->GetDetectorData()->At(i)); }
		long SBCOffset = vData[0]->GetSBCTSsec() - start;
		//long SBCOffset = -15049;
		//vRun->ListRunBits();	// Dump run bits

		// Do a very rough estimate of the number of LED events 
		// and output the frequency.
		int multipCounter = 0;
		int highestmultip = 0;
		for (int i = 0; i < (int)vEntries; i++)
		{
			v->GetEntry(i);
			if (mVeto > 26) multipCounter++;
			if (mVeto < highestmultip) highestmultip = mVeto;
		}
		//printf("Run: %i  Frequency: %.2f\n",run,(double)multipCounter/duration);
		if (multipCounter == 0) printf("No LED's!  Run: %i  Highest multiplicity found: %i\n",run,highestmultip);
	}

}
