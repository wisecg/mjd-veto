void vetoSimple()
{
	int run = 6960;
	GATDataSet ds(run);

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

	int thresh[32] = {134,126,107,102,163,122,121,114,128,101,120,114,105,113,103,9999,
		67,164,106,135,88,93,142,73,136,117,93,123,159,176,130,96};

	for,(int i = 0; i < vEntries; i++) 
	{
		v->GetEntry(i);

		// Fill MJVetoEvent object
		MJVetoEvent veto;
		veto.SetSWThresh(thresh);	
    	veto.WriteEvent(i,vRun,vEvent,vBits,run,true); // true: force-write event with errors.
	    	
    	printf("i: %i  multiplicity: %i  total energy: %i\n",i,veto.GetMultip(),veto.GetTotE());
	}
}