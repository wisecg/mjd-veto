// Check if the veto data is getting any high multiplicity events.
// C. Wiseman.

#include "vetoScan.hh"

void vetoLEDFinder(string file) 
{
	// From muFinder:
	// LED Cut Parameters (C-f "Display Cut Parameters" below.)
	// double LEDWindow = 0.1;	
	// int LEDMultipThreshold = 10;  // "multipThreshold" = "highestMultip" - "LEDMultipThreshold"
	// int LEDSimpleThreshold = 5;   // used when LED frequency measurement is bad.

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

		// Do a very rough estimate of the number of LED events 
		// and output the frequency.
		int multipCounter = 0;
		int highestmultip = 0;
		for (int i = 0; i < (int)vEntries; i++)
		{
			v->GetEntry(i);
			if (mVeto > 26) multipCounter++;
			if ((int)mVeto > highestmultip) highestmultip = mVeto;
		}
		printf("Run: %i  Approx Freq: %.2f  Max multip: %i\n",run,(double)multipCounter/duration,highestmultip);
		if (multipCounter == 0) printf("No LED's!  Run: %i  Highest multiplicity found: %i\n",run,highestmultip);
	}

}
