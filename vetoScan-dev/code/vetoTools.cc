// ==================================================
// Processing Functions
// ==================================================

#include "vetoScan.hh"
using namespace std;

void Test()
{
	string s1("../somepath/somemorepath/somefile.ext");
	string s2("..\\somepath\\somemorepath\\somefile.ext");
	cout << s1.substr(0, s1.find_last_of("\\/")) << endl;
	cout << s2.substr(0, s2.find_last_of("\\/")) << endl;

	string s3("./runs/P3K93_VetoOnlyRuns.txt");
	cout << s3.substr(s3.find_last_of("\\/")+1,string::npos) << endl;
}

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

int GetNumFiles(string arg)
{
	int run = 0; 

	ifstream InputList;
	InputList.open(arg.c_str());
	if(!InputList.good()) {
    	cout << "Couldn't open " << arg << " !" << endl;
    	return 0;
    }

	int filesToScan = 0;
  	while(!InputList.eof()) { 
  		InputList >> run; 
  		filesToScan++; 
  	}
  	cout << "Scanning " << filesToScan << " files." << endl;	  	
  	InputList.close();
  	
  	return filesToScan;
}

// ROOT color wheel: https://root.cern.ch/root/html/TColor.html
int color(int i)
{
	if (i == 0) return kRed-2;  
	if (i == 1) return kBlue;  
	if (i == 2) return kBlue+6;
	if (i == 3) return kAzure; 
	if (i == 4) return kMagenta;
	if (i == 5) return kGreen;
	if (i == 6) return kSpring;
	return 0;
}

// For tagging plane-based coincidences.
// This uses a zero-indexed map.
int PanelMap(int i){
	// 0: Lower Bottom
	// 1: Upper Bottom
	// 2: Inner Top
	// 3: Outer Top
	// 4: Inner North
	// 5: Outer North
	// 6: Inner South
	// 7: Outer South
	// 8: Inner West
	// 9: Outer West
	// 10: Inner East
	// 11: Outer East 
			
	if 		(i == 0) return 0;  // L-bot 1
	else if (i == 1) return 0;  
	else if (i == 2) return 0;
	else if (i == 3) return 0;
	else if (i == 4) return 0;
	else if (i == 5) return 0;  // L-bot 6

	else if (i == 6) return 1;  // U-bot 1
	else if (i == 7) return 1;
	else if (i == 8) return 1;
	else if (i == 9) return 1;
	else if (i == 10) return 1;
	else if (i == 11) return 1; // U-bot 6

	else if (i == 17) return 3; // Top outer
	else if (i == 18) return 3; // Top outer
	else if (i == 20) return 2; // Top inner
	else if (i == 21) return 2; // Top inner

	else if (i == 15) return 5; // North outer
	else if (i == 16) return 5; // North outer
	else if (i == 19) return 4; // North inner
	else if (i == 23) return 4; // North inner

	else if (i == 24) return 6; // South inner
	else if (i == 25) return 7; // South outer
	else if (i == 26) return 6; // South inner
	else if (i == 27) return 7; // South outer
	
	else if (i == 12) return 8; // West inner
	else if (i == 13) return 8; // West inner
	else if (i == 14) return 9; // West outer
	else if (i == 22) return 9; // West outer
	
	else if (i == 28) return 10; // East inner
	else if (i == 29) return 11; // East outer
	else if (i == 30) return 10; // East inner
	else if (i == 31) return 11; // East outer

	else return -1;
}

// MJVetoEvent "error filter" - analysis codes skip events which fail
// the cuts here.  Returns true if there is a bad error.
bool CheckForBadErrors(MJVetoEvent veto, int entry, int isGood, bool verbose) 
{
	bool badError = false;

	if (isGood != 1) {
		int error[18] = {0};
		veto.UnpackErrorCode(isGood,error);
		
		// search for particular errors and set a flag.
		for (int q=0; q<18; q++) 
		{
			// 4: don't skip bad-scaler events
			// 10: P3K93: don't skip "event count doesn't match ROOT entry" events
			// 7 & 11: don't skip hw count mismatches or scaler != qdc event count.  Due to problems with continuous running mode.

			// if (q!=4 && q!=10 && error[q] == 1) badError = true;
			if (q != 4 && q != 7 && q != 10 && q != 11 && q != 12 && error[q] == 1) badError = true;
		}
		
		if (badError) {
			if (verbose == true) {
				cout << "Skipped Entry: " << entry << endl;
				veto.Print();
				cout << endl;
			}
			return badError;
		}
	}
	else return badError;	// no bad errors found

	return false;
}

// Place threshold 35 qdc above pedestal location.
// For now, "deactivate" does nothing.
int FindQDCThreshold(TH1F *qdcHist, int panel, bool deactivate) 
{
	// new method
	int manualMaxBin = -1;
	int manualMaxVal = -1;
	int locBin = -1;
	int locVal = -1;
	qdcHist->GetXaxis()->SetRangeUser(0,500);
	int nbinsx = qdcHist->GetXaxis()->GetNbins();
	for (int i = 0; i < nbinsx; i++)
	{
		locBin = qdcHist->GetBin(i);
 		locVal = qdcHist->GetBinContent(locBin);
		// if (locVal > 0) printf("panel: %i  bin: %i  value: %i  manualMaxBin: %i\n",panel,locBin,locVal,manualMaxBin);
		if (locVal > manualMaxVal) 
		{
			manualMaxBin = locBin;
			manualMaxVal = locVal;
		}
	}
	// old method
	// int firstNonzeroBin = qdcHist->FindFirstBinAbove(10,1);
	// qdcHist->GetXaxis()->SetRange(firstNonzeroBin-10,firstNonzeroBin+50);
	// int oldMaxBin = qdcHist->GetMaximumBin();
	// int oldMaxVal = qdcHist->GetBinContent(oldMaxBin);
	// printf("Panel: %i  manual - maxBin: %i  maxVal: %i\t old - maxBin: %i  maxVal: %i\n",panel,manualMaxBin,manualMaxVal,oldMaxBin,oldMaxVal);

	// Set QDC Threshold (using new method) 
	double xval = qdcHist->GetXaxis()->GetBinCenter(manualMaxBin);
	return xval+35;
}

int* GetQDCThreshold(string file, int *arr, string name)
{
	string Name = "";

	// Can either use file name or specify a threshold with "name".
	if (name == "") 
	{
		// strip off path and extension
		Name = file;
		Name.erase(Name.find_last_of("."),string::npos);
		Name.erase(0,Name.find_last_of("\\/")+1);
	}
	else {
		Name = name;
	}

	// open the threshold dictionary
	// generated by vetoThreshFinder
	//
	ifstream InputList("vetoSWThresholds.txt");
	if(!InputList.good()) {
    	cout << "Couldn't open vetoSWThresholds.txt " << endl;
    	return NULL;
    }
	int result[32] = {0};
	string buffer;
	char val[200];
	bool foundRange = false;
	while (!InputList.eof() && !foundRange)
	{
		InputList >> buffer;
		if (buffer == Name) {
			cout << "Found SW threshold values for: " << Name << endl;
			foundRange = true;
			for (int i = 0; i < 32; i++) {
				InputList >> val;
				result[i] = atoi(val);
			}
		}
	}
	if (!foundRange) {
		cout << "Didn't find SW threshold values for this range. \n Using defaults..." << endl;
		for (int j=0;j<32;j++) result[j]=500;
	}
	memcpy(arr,result,sizeof(result));

	return arr;
}

double InterpTime(int entry, vector<double> times, vector<double> entries, vector<bool> badScaler)
{
	if ((times.size() != entries.size()) || times.size() != badScaler.size()) 
	{
		cout << "Vectors are different sizes!\n";  
		if (entry >= (int)times.size()) 
			cout << "Entry is larger than number of entries in vector!\n";
		return -1; 
	}
	
	double iTime = 0;

	double lower = 0;
	double upper = 0;
	if (!badScaler[entry]) iTime = times[entry];
	else 
	{
		for (int i = entry; i < (int)entries.size(); i++) 
		{
			if (badScaler[i] == 0) { upper = times[i]; break; }
		}
		for (int j = entry; j > 0; j--)
		{
			if (badScaler[j] == 0) { lower = times[j]; break; }
		}

		iTime = (upper + lower)/2.0;
	}

	return iTime;
}