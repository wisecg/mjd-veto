#include "vetoScan.hh"

// adapted from Rene Brun's copytree.C
void GrabVetoTree(string file) 
{
	int run = 0;

	ifstream InputList;
  	InputList.open(file.c_str());
	while(!InputList.eof())
	{
		InputList >> run;
	
		char File[200];
		sprintf(File,"/global/project/projectdirs/majorana/data/mjd/surfprot/data/built/P3END/OR_run%u.root",run); 

		cout << "Scanning run " << run << endl;

		TFile *oldfile = new TFile(File);
		TTree *oldtree = (TTree*)oldfile->Get("VetoTree");
		
		MGTBasicEvent *vEvent = new MGTBasicEvent(); 
		
		if (oldtree->GetListOfBranches()->FindObject("vetoEvent")) 
			oldtree->SetBranchAddress("vetoEvent",&vEvent);
		else { cout << "Couldn't find branch: \"vetoEvent\"" << endl; break; }

		char Out[200];
		sprintf(Out,"./output/OR_run%u.root",run);
		TFile *newfile = new TFile(Out,"recreate");
		TTree *newtree = oldtree->CloneTree();

		newtree->Print();
		newfile->Write();
		delete oldfile;
		delete newfile;
	}
}