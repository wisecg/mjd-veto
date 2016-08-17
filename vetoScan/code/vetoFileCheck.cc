/*
	vetoFileCheck.C
	Clint Wiseman, USC/Majorana
	August 2015.

	Macro to check a list of run numbers and determine if the built files exist.
	Will also check if files have been blinded and aren't readable.
*/

#include "vetoScan.hh"

using namespace std;

void vetoFileCheck(string Input, string partNum, bool checkBuilt, bool checkGat, bool checkGDS)
{
	// Input a list of run numbers
	ifstream InputList(Input.c_str());
	if(!InputList.good()) {
    	cout << "Couldn't open " << Input << endl;
    	return;
    }

	int run;
	char GATFile[200];
	char BuiltFile[200];
	char path[200];
	double durationTotal = 0;
	while(!InputList.eof())
	{
		InputList >> run;

		if (partNum != "") {
			sprintf(path,"/global/project/projectdirs/majorana/data/mjd/surfmjd/data");
			sprintf(BuiltFile,"%s/built/%s/OR_run%u.root",path,partNum.c_str(),run); 
			sprintf(GATFile,"%s/gatified/%s/mjd_run%u.root",path,partNum.c_str(),run); 
		}
		else cout << "Warning!  Empty part number!" << endl;

		if (checkBuilt && partNum != "") {
			TFile *f1 = new TFile(BuiltFile);
    		f1->Close();
    		delete f1;
    	}
    	if (checkGat && partNum != "") {
    		TFile *f2 = new TFile(GATFile);
    		f2->Close();
    		delete f2;
    	}
    	if (checkGDS){
    		GATDataSet *ds = new GATDataSet(run);
    		cout << ds->GetRunTime() << endl;
    		// TChain *b = ds->GetBuiltChain();
    		// cout << "Built file: " << b->GetEntries() << " entries\n";
    		// TChain *g = ds->GetGatifiedChain();
    		// cout << "Gatified file: " << g->GetEntries()<< " entries\n";
    		delete ds;
    	}

		// Check also that the duration is not corrupted!
		if (checkBuilt) {
			Float_t duration = 0;
			TChain *MGTree = new TChain("MGTree");
			MGTree->AddFile(BuiltFile);
			MJTRun *MyRun = new MJTRun();
			MGTree->SetBranchAddress("run",&MyRun);
			MGTree->GetEntry(0);
			duration = MyRun->GetStopTime() - MyRun->GetStartTime();
			//if (duration >= 3595 && duration <= 3605) cout << run << endl;
			if (duration <= 0 || duration > 4000 ) {
				printf("\nRun %i has duration %.0f, skipping file!\n\n",run,duration);
				continue;
			}
			durationTotal+=duration;
			delete MGTree;
			delete MyRun;
		}
    }
	cout << "Total duration: " << durationTotal << " seconds." << endl; 
}
