/*
	CheckFiles.C
	Clint Wiseman, USC/Majorana
	August 2015.

	Macro to check a list of run numbers and determine if the built files
	exist.
	Will also check if files have been blinded and aren't readable.

	Usage: 
	root[0] .X CheckFiles.C 
	root[0] .X CheckFiles.C ("InputFile_WithNoExtension")
*/

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>

#include <TChain.h>
#include <TFile.h>

using namespace std;

void CheckFiles(string Input = ""){

	int mode = 1; // switch: 0 for local files, 1 for pdsf files

	// Input a list of run numbers
	if (Input == "") Char_t InputName[200] = "builtVeto_DebugList";
	else Char_t InputName[200] = Input.c_str();
	Char_t InputFile[200];
	sprintf(InputFile,"%s.txt",InputName);
	ifstream InputList;
	InputList.open(InputFile);
	Char_t GATFile[200];
	Char_t BuiltFile[200];

	Int_t run;
	while(!InputList.eof()){

		// initialize 
		InputList >> run;
		if (mode==0) sprintf(BuiltFile,"~/dev/datasets/builtVeto/OR_run%i.root",run);
		else if (mode==1) sprintf(BuiltFile,"/global/project/projectdirs/majorana/data/mjd/surfmjd/data/built/P3JDY/OR_run%u.root",run); 

		if (mode==0) sprintf(GATFile,"~/dev/datasets/builtVeto/mjd_run%i.root",run);
		else if (mode==1) sprintf(GATFile,"/global/project/projectdirs/majorana/data/mjd/surfmjd/data/gatified/P3JDY/mjd_run%u.root",run); 


		//cout << "Checking for built & gatified runs, run " << run << endl;

		// if file doesn't exist, ROOT will fail to open it.
		TFile *f1 = new TFile(BuiltFile);
    	f1->Close();

		// Check also that the duration is not corrupted!
		Float_t duration = 0;
		TChain *MGTree = new TChain("MGTree");
		MGTree->AddFile(BuiltFile);
		MJTRun *MyRun = new MJTRun();
		MGTree->SetBranchAddress("run",&MyRun);
        MGTree->GetEntry(0);
        duration = MyRun->GetStopTime() - MyRun->GetStartTime();
        if (duration <= 0 || duration > 4000 ) {
        	printf("\nRun %i has duration %.0f, skipping file!\n\n",run,duration);
        	continue;
        }

    	TFile *f2 = new TFile(GATFile);
    	f2->Close();
    }

}
