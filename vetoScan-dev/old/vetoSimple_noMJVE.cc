#include "vetoScan.hh"

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
	// v->GetEntry(0);
	// MJTVetoData *vData[32];
	// for (int i=0; i<32; i++) { vData[i] = dynamic_cast<MJTVetoData*>(vEvent->GetDetectorData()->At(i)); }

	// loop over entries
	for (int i = 0; i < vEntries; i++) {
		v->GetEntry(i);

		MJTVetoData *vData[32];
		for (int i=0; i<32; i++) { vData[i] = dynamic_cast<MJTVetoData*>(vEvent->GetDetectorData()->At(i)); }

		// get timestamps
		vData[0]->GetTimeStamp()

		vData[]

		// get qdc value

		// get errors
	}
}