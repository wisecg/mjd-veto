#include "vetoScan.hh"

void muDisplayList(string file)
{
	// Can handle more than one input file
	TChain *v = new TChain("vetoEvent");
	v->AddFile(file.c_str());

	// Set up output files
	string Name = file;
	Name.erase(Name.find_last_of("."),string::npos);
	Name.erase(0,Name.find_last_of("\\/")+1);
	string outFile = "./output/vList_"+Name+".txt";
	cout << "Writing veto hit list: " << outFile << endl;

	// Output a veto hit list
	ofstream hitList(outFile.c_str());

	// Initialize output from muFinder
	MJVetoEvent *event = NULL;
	int CoinType[32] = {0};
	int CutType[32] = {0};
	double LEDfreq = 0;
	double LEDrms = 0;
	int highestMultip = 0;
	Long64_t start = 0;
	Long64_t stop = 0;
	double xTime = 0;
	double x_deltaT = 0;
	v->SetBranchAddress("events",&event);
	v->SetBranchAddress("CoinType[32]",CoinType);
	v->SetBranchAddress("CutType[32]",CutType);
	v->SetBranchAddress("LEDfreq",&LEDfreq);
	v->SetBranchAddress("LEDrms",&LEDrms);
	v->SetBranchAddress("highestMultip",&highestMultip);
	v->SetBranchAddress("start",&start);
	v->SetBranchAddress("stop",&stop);
	v->SetBranchAddress("xTime",&xTime);
	v->SetBranchAddress("x_deltaT",&x_deltaT);
	int vEntries = v->GetEntries();
	cout << "Found " << vEntries << " entries.\n";
	v->GetEntry(0);
	int numPanels = highestMultip;

	for (int i = 0; i < vEntries; i++)
	{
		v->GetEntry(i);

		// hit list format:
		// run entry QEC time qdc1 ... qdc32

		if (CoinType[0]==1) 
		{
			hitList << event->GetRun() << " " << i << " " << event->GetSEC() << " " << xTime << " ";
			for (int j=0;j<numPanels;j++)  
			{
				if (event->GetQDC(j) >= event->GetSWThresh(j))
					hitList << event->GetQDC(j) << " ";
				else 
					hitList << 0 << " ";
			}
			hitList << endl;
		}
	}

}