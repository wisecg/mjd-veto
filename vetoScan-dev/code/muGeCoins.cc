// Find muon-Ge coincidences.
// Clint Wiseman, USC/Majorana

#include "vetoScan.hh"

using namespace std;

void muGeCoins(string Input) 
{
	cout << "This is going to get complicated.\n";
/*
	// Input a list of run numbers
	ifstream InputList(Input.c_str());
	if(!InputList.good()) {
		cout << "Couldn't open " << Input << endl;
		return;
	}

	// Input a ROOT file with veto information
	string Name = Input;
	Name.erase(Name.find_last_of("."),string::npos);
	Name.erase(0,Name.find_last_of("\\/")+1);
	char vFileName[200];
	cout << " Name is: " << Name << endl;
	sprintf(vFileName,"./output/%s.root",Name.c_str());
	TFile *f1 = new TFile(vFileName);       
    TTree *vEvent = (TTree*)f1->Get("vetoEvent");
    MJVetoEvent *veto;
	Long64_t start;
	Long64_t stop;
	int CoinType[32];
	int CutType[32];
	int highestMultip = 0;
	double LEDfreq = 0;
	double LEDrms = 0;
	double xTime = 0;
	double x_deltaT = 0;
    vEvent->SetBranchAddress("events",&veto);
	vEvent->SetBranchAddress("CoinType[32]",&CoinType);
	vEvent->SetBranchAddress("CutType[32]",&CutType);
	vEvent->SetBranchAddress("LEDfreq",&LEDfreq);
	vEvent->SetBranchAddress("LEDrms",&LEDrms);
	vEvent->SetBranchAddress("highestMultip",&highestMultip);
	vEvent->SetBranchAddress("start",&start);
	vEvent->SetBranchAddress("stop",&stop);
	vEvent->SetBranchAddress("xTime",&xTime);
	vEvent->SetBranchAddress("x_deltaT",&x_deltaT);
	Long64_t vEntries = vEvent->GetEntries();
	Long64_t vEntry = 0;
	cout << "Veto File has " << vEntries << " entries.\n";

  	// Create TEntryList with a particular coincidence type, possibly from command line.
  	vEvent->Draw(">>elist", "CoinType[0]==1", "entrylist");
    TEntryList *elist = (TEntryList*)gDirectory->Get("elist");
    elist->Print("");	// "": print the name of the tree and file, "all": print all the entry numbers
	long listEntries = elist->GetN();  
	printf("Created TEntryList from vEvent with %li entries.\n",listEntries);

	// Output a ROOT file.
	Char_t OutputFile[200];
	sprintf(OutputFile,"./output/MGC_%s.root",Name.c_str());
	TFile *RootFile = new TFile(OutputFile, "RECREATE"); 	
	TH1::AddDirectory(kFALSE); // Global flag: "When a (root) file is closed, all histograms in memory associated with this file are automatically deleted."
	bool isCoin = false;	// coin with veto
	int gEntry = 0;
	double gRun = 0;
	double gTimeSec = 0;    
	double deltaT = 0;
	double gESum = 0;
	GATMIJ *mij = new GATMIJ();
	MJVetoEvent out;
	MGTEvent* gEvent = new MGTEvent();    
  	vector<double>* trapECal = 0;
	vector<double>* timestamp = 0;
	vector<double>* channel = 0;
    vector<double>* blrwfFMR50 = 0;
    vector<double>* rawWFMin = 0;
    vector<double>* rawWFMax = 0;
	TTree *geCoins = new TTree("geCoins","Ge Events");
    geCoins->Branch("vetoEvent","MJVetoEvent",&out,32000,1);
    geCoins->Branch("isCoin",&isCoin);	
    geCoins->Branch("deltaT",&deltaT,"deltaT/D");
	geCoins->Branch("geEvent","MGTEvent",&gEvent,32000,0);
	geCoins->Branch("mij",&mij);
    geCoins->Branch("gRun",&gRun);
    geCoins->Branch("gEntry",&gEntry);	
    geCoins->Branch("gTimeSec",&gTimeSec);  
    geCoins->Branch("gESum",&gESum);		
    geCoins->Branch("trapECal",&trapECal);
    geCoins->Branch("timestamp",&timestamp);
    geCoins->Branch("channel",&channel);
	geCoins->Branch("blrwfFMR50",&blrwfFMR50);
	geCoins->Branch("rawWFMin",&rawWFMin);
	geCoins->Branch("rawWFMax",&rawWFMax);

	cout << "9\n";

    // Analysis Parameters
	// double GClock = 100000000;	// gretina clock speed
	// double loWindow = 0.1;  
	// double hiWindow = 3602; 		
	// double e_cut = 2650;	// =50 keV, eliminate low-energy noise. =2650 keV, skim high-E events.
	// double e_floor = 4;     // lowest energy to include (not currently used)

    // Time windows around an event on the TEventList 
  	// Watch out for double counting!
  	// The prototype was 4*sigma = 12 seconds
  	// Default: change to 7 seconds b/c LED period is 7sec.
  	// if adjustWindow is set, will adjust time window depending on scaler being good/bad.
  	// bool adjustWindow = true;	
  	// double hi = 7; // sec
  	// double lo = 7; 

    // Loop over files, control it with the TEntryList
  	// We scan ALL the vEvent entries with the same run number 
  	// while we have the corresponding file loaded.
  	int run = 0;
  	cout << "4\n";
	while(!InputList.eof())
	{
		InputList >> run;
		printf("\n ======= Scanning new file: run %i ======\n",run);

		// Grab all entries in the TEntryList with veto.run == run.
		list<int> VEntThisRun;
		for (int q = 0; q < listEntries; q++) {
			vEntry = elist->GetEntry(q);
			cout << "Got vEntry: " << vEntry << endl;
			
			if (vEntry < vEntries) vEvent->GetEntry(vEntry);
			else cout << " you're a dumbass! \n";
		// 	cout << "Got VEvent " << endl;	
		// 	if(veto->run == run && vEntry != -1) {
		// 		VEntThisRun.push_back((int)vEntry);
		// 		printf("V-List: run:%-5i  vEntry:%lli  xTime:%-5.5f  Bad:%i  CT0:%i  CT1:%i  CT2:%i\n"
		// 			,run,vEntry,xTime,veto->badScaler,CoinType[0],CoinType[1],CoinType[2]);
		// 	}
		// 	if (vEntry == -1) break; // done with elist
		}
		// if (VEntThisRun.size()==0) {
		// 	printf("  Found 0 vEvent entries.  Continuing ...\n");
		// 	// continue; 
		// }
		// else printf("Found %i entries in the vEvent.\n",(int)VEntThisRun.size());
		// // int coins = 0;
	}

	// all done!
	RootFile->Close();
*/
}