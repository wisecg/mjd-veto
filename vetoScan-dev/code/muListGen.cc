// Generate a list of muon-induced events from muFinder ROOT output.
// Used in the DEMONSTRATOR Veto Cut.
//
// Clint Wiseman, USC/Majorana
// 4/22/16

#include "vetoScan.hh"
using namespace std;
void muListGen(string arg) 
{
	TChain *v = new TChain("vetoEvent");
	
	// Input chain: can handle more than one file.
	// v->AddFile(arg.c_str());
	v->AddFile("./output/DS1_01.root");
	v->AddFile("./output/DS1_02.root");
	v->AddFile("./output/DS1_03.root");
	v->AddFile("./output/DS1_04.root");
	v->AddFile("./output/DS1_05.root");
	v->AddFile("./output/DS1_06.root");

	// Output: Text file muon list (used in skim files)
	string Name = "DS1";
	string outName = "./output/MuonList_"+Name+".txt";
	ofstream MuonList(outName.c_str());

	// initialize muFinder ROOT output
	MJVetoEvent *event = NULL;
	double LEDfreq=0,LEDrms=0,LEDWindow=0,xTime=0,x_deltaT=0,x_LEDDeltaT=0;
	int multipThreshold=0,highestMultip=0,LEDMultipThreshold=0,LEDSimpleThreshold=0,PlaneHitCount=0;
	Long64_t start=0,stop=0;
	int CoinType[32] = {0};
	int CutType[32] = {0};
	int PlaneHits[12] = {0};
	int PlaneTrue[12] = {0};
	v->SetBranchAddress("events",&event);
	// v->SetBranchAddress("rEntry",&rEntry);	// not implemented for DS1 yet.
	v->SetBranchAddress("LEDfreq",&LEDfreq);
	v->SetBranchAddress("LEDrms",&LEDrms);
	v->SetBranchAddress("multipThreshold",&multipThreshold);
	v->SetBranchAddress("highestMultip",&highestMultip);
	v->SetBranchAddress("LEDWindow",&LEDWindow);
	v->SetBranchAddress("LEDMultipThreshold",&LEDMultipThreshold);
	v->SetBranchAddress("LEDSimpleThreshold",&LEDSimpleThreshold);
	v->SetBranchAddress("start",&start);
	v->SetBranchAddress("stop",&stop);
	v->SetBranchAddress("xTime",&xTime);
	v->SetBranchAddress("x_deltaT",&x_deltaT);
	v->SetBranchAddress("x_LEDDeltaT",&x_LEDDeltaT);
	v->SetBranchAddress("CoinType[32]",CoinType);
	v->SetBranchAddress("CutType[32]",CutType);
	v->SetBranchAddress("PlaneHits[12]",PlaneHits);
	v->SetBranchAddress("PlaneTrue[12]",PlaneTrue);
	v->SetBranchAddress("PlaneHitCount",&PlaneHitCount);
	int vEntries = v->GetEntries();
	cout << "Found " << vEntries << " entries.\n";
	
	// loop over entries
	int counter = 0;
	for (long i = 0; i < vEntries; i++) 
	{
		v->GetEntry(i);

		// Write a text file
		char buffer[200];
		if (CoinType[1] || CoinType[0]) 
		{
			counter++;
			int type;
			if (CoinType[0]) type = 1;
			if (CoinType[1]) type = 2;
			sprintf(buffer,"%i %lli %.8f %i %i\n",event->GetRun(),start,xTime,type,event->GetBadScaler());
			MuonList << buffer;
			// cout << buffer;
		}
	}
	cout << "Found " << counter << " candidates\n";

	// end of routine.
	MuonList.close();
}