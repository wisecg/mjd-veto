/*
	builtVeto.C
	Clint Wiseman, USC/Majorana
	June 2015.
	
 => This code can be run on PDSF or locally.  It takes a .txt file of run numbers
	as an input argument, and uses the name of the text file to generate output
	in the folder ./output

 => Recommended: When scanning a new input file of run numbers on PDSF, run CheckFiles.C 
 	to make sure files exist and have not been blinded.  
 	This code has a tendency to quit unexpectedly when it encounters a "bad" file.
 
 => builtVeto uses the same sorting method developed in builtVetoSimple.  
 	It is then mainly used for plotting scaler corruption and multiplicity of events.
	Generates ROOT files of histograms, allowing one to look at run-by-run
	scaler corruption in time, and run-by-run multiplicity to look for 
	changes in the system.

 => builtVetoCal.C is a bit more advanced, and (among other things) uses
	custom threshold values for each veto panel.

	Usage:
	CINT: root[0] .X builtVeto.C ("Filename_list_of_run_numbers")  <--- NO .TXT extension.
	bash: root -b -q -l builtVeto.C ("The_filename_without_extension")
*/

const int numPanels = 32;

// Structure to hold a complete matched veto event 
struct VEvent {
    Int_t run;				// run number
    Int_t vEnt; 			// qdc eventCount (i.e. entry)
    Int_t sEnt; 			// scaler entry
    Int_t card[3];			// 0: scaler ID, 1: qdc card 1 (panels 1-16), 2: qdc card 2 (panels 17-32)
    Float_t sTime;			// scaler time
    Bool_t sTimeBad;		// bad scaler entry
    Float_t lTime;			// LED time
    Float_t eTime;			// entry time
    Int_t QDC[numPanels];	// store qdc values for both cards in the same array 
    Int_t IsUnderThreshold[numPanels];
    Int_t IsOverflow[numPanels];
};

// Structure to hold a qdc entry
struct QEvent{
    Int_t runNumber;
    UInt_t crate;
    UInt_t card;
    UInt_t EventCount;
	Int_t QDC[16];
    Int_t IsUnderThreshold[16];
    Int_t IsOverflow[16];
};


void builtVeto(string Input = ""){

	int mode = 1;				// switch: 0 for local files, 1 for pdsf files
	int card1 = 13;
	int card2 = 18;

	// Input a list of run numbers
	if (Input == "") Char_t InputName[200] = "builtVeto_DebugList";
	else Char_t InputName[200] = Input.c_str();
	Char_t InputFile[200];
	sprintf(InputFile,"%s.txt",InputName);
	ifstream InputList;
	InputList.open(InputFile);
	Char_t TheFile[200];

	// Set up output file(s)
	Char_t OutputFile[200];
	sprintf(OutputFile,"%s.root",InputName);
	TFile *RootFile = new TFile(OutputFile, "RECREATE");	
  	TH1::AddDirectory(kFALSE); // Global flag: "When a (root) file is closed, all histograms in memory associated with this file are automatically deleted."

	//=== Global counters / variables / plots ===

		Int_t run = 0;
		Float_t duration = 0;
		
		// get number of files in dataset for the TGraph
		Int_t filesToScan = 0;
		Int_t filesScanned = 0;
	  	while(!InputList.eof()) { InputList >> run; filesToScan++; }
	  	cout << "Scanning " << filesToScan << " files." << endl;	  	
	  	InputList.close();
	  	InputList.open(InputFile);
	  	run=0;
	 	TGraph *SCorruption = new TGraph(filesToScan);
	 	Int_t BadTSTotal = 0;

	 	TH1D *TotalCorruptionInTime = new TH1D("TotalCorruptionInTime","corrupted entries during run (entry method)",(Int_t)3600/5,0,3600);
	 	TotalCorruptionInTime->GetXaxis()->SetTitle("time (5 sec / bin)");
	 	Bool_t PlotCorruptedEntries = true; // flag for plotting corrupted entries in time for EACH RUN

	 	TH1D *TotalMultiplicity = new TH1D("TotalMultiplicity","Events over threshold",32,0,32);
	 	TotalMultiplicity->GetXaxis()->SetTitle("number of panels hit");
	 	Bool_t PlotMultiplicity = true;	// flag to plot multiplicity for EACH RUN

		
	//=== End ===



	// Loop over files
	while(!InputList.eof()){

		// initialize 
		InputList >> run;
		if (mode==0) sprintf(TheFile,"~/dev/datasets/builtVeto/OR_run%i.root",run);
		else if (mode==1) sprintf(TheFile,"/global/project/projectdirs/majorana/data/mjd/surfmjd/data/built/P3JDY/OR_run%u.root",run); 
		TChain *VetoTree = new TChain("VetoTree");
		VetoTree->AddFile(TheFile);
		TChain *MGTree = new TChain("MGTree");
		MGTree->AddFile(TheFile);
		MJTRun *MyRun = new MJTRun();
		MGTree->SetBranchAddress("run",&MyRun);
		Long64_t nentries = VetoTree->GetEntries();
		MGTBasicEvent *b = 0;
		//MGTBasicEvent b;
		VetoTree->SetBranchAddress("vetoEvent",&b);
		const int vd_size = 16;	
		MJTVetoData *vd[vd_size];
		VetoTree->GetEntry(0);
		for (int i=0; i<vd_size; i++) { vd[i] = dynamic_cast<MJTVetoData*>(b->GetDetectorData()->At(i)); }
        MGTree->GetEntry(0);
        duration = MyRun->GetStopTime() - MyRun->GetStartTime();

    	//=== Single-file counters / variables / plots

			QEvent qdc = {0};
			QEvent prevqdc = {0};
				prevqdc.EventCount=-1;
			VEvent veto = {0};
			Bool_t EventMatch = false;

			Int_t BadTSInFile = 0;
			Float_t corruption = 0;
			if (PlotCorruptedEntries) {
				TH1D *CorruptionInTime = new TH1D("CorruptionInTime","corrupted entries during run (entry method)",(Int_t)duration/5,0,(Int_t)duration);
				CorruptionInTime->GetXaxis()->SetTitle("time (5 sec / bin)");
			}
			Int_t numPanelsHit = 0;
			if (PlotMultiplicity) {
				TH1D *OneRunMultiplicity = new TH1D("multiplicity","multiplicity of veto entries",32,0,32);
  				OneRunMultiplicity->GetXaxis()->SetTitle("number of panels hit");
  			}
			
		//=== End ===

		// Loop over VetoTree entries
		printf("Now scanning run %i: %lli entries, %.2f sec.  \n",run,nentries,duration);
		for (int i = 0; i < nentries; i++) {
			VetoTree->GetEntry(i);



			// move data into QEvent structure
			qdc.runNumber=run;
			qdc.crate=vd[0]->GetCrate();
			qdc.card=vd[0]->GetCard();
			qdc.EventCount=vd[0]->GetEventCount();
			int k = 0;
			for (int j = 0; j<16; j++)	{
				k = vd[j]->GetChannel();
				qdc.QDC[k]=vd[j]->GetAmplitude();
				qdc.IsUnderThreshold[k]=vd[j]->IsUnderThreshold();
				qdc.IsOverflow[k]=vd[j]->IsOverflow();
			}
			
			// check qdc data after moving into QEvent structure
			//printf("QDC -- run: %i  Entry: %i  crate: %i  card: %i  EventCount: %i \n",qdc.runNumber,i,qdc.crate,qdc.card,qdc.EventCount);
			//cout << "QDC: "; for (int k = 0; k<16; ++k) { cout << qdc.QDC[k] << "  "; } cout << endl;
			//cout << "UTh: "; for (int k = 0; k<16; ++k) { cout << qdc.IsUnderThreshold[k] << "  "; } cout << endl;
			//cout << "Ovr: "; for (int k = 0; k<16; ++k) { cout << qdc.IsOverflow[k] << "  "; } cout << endl;	
			
			// check entry numbers
			//printf("Entry: %i  card: %i  EventCount: %i",i,qdc.card,qdc.EventCount);
			//printf("  ||  ScalerCount: %i  TimeStamp: %.5f  IsBadTs: %i \n",
			//	vd[0]->GetScalerCount(),vd[0]->GetTimeStamp()/1E8,vd[0]->IsBadTS());
			
			// set flag if current qdc entry has same EventCount as previous one.
			EventMatch = false;
			if (qdc.EventCount == prevqdc.EventCount) { 
				EventMatch = true; 
				//printf("EventMatch true.  qdc.EventCount:%i  prevqdc.EventCount:%i \n",qdc.EventCount,prevqdc.EventCount);
			}
			else if (abs(qdc.EventCount - prevqdc.EventCount) > 1 && i > 2) {
				printf(" EventCount mismatch! Run:%i current:%i previous:%i card:%i prev.card:%i  Breaking at %.0f%% through file.\n",run,i,qdc.card,prevqdc.card,((Float_t)i/nentries)*100);
				break;
			}
			
			if (EventMatch) {
				
				// move matching events into VEvent structure and incorporate scaler info.
				veto.run = qdc.runNumber;
				veto.vEnt = qdc.EventCount;
				veto.sEnt = vd[0]->GetScalerCount();
				veto.sTime = vd[0]->GetTimeStamp()/1E8;
				veto.sTimeBad = vd[0]->IsBadTS();
				veto.eTime = ((Float_t)i/nentries)*duration;	
				// case 1
				if (prevqdc.card==card1 && qdc.card==card2) {
					veto.card[0]=vd[0]->GetScalerID();
					veto.card[1]=prevqdc.card;
					veto.card[2]=qdc.card;
					for (int k = 0; k<16; k++) {
						veto.QDC[k]=prevqdc.QDC[k];
						veto.QDC[16+k]=qdc.QDC[k];
						veto.IsUnderThreshold[k]=prevqdc.IsUnderThreshold[k];
						veto.IsUnderThreshold[16+k]=qdc.IsUnderThreshold[k];
						veto.IsOverflow[k]=prevqdc.IsOverflow[k];
						veto.IsOverflow[16+k]=qdc.IsOverflow[k];
					}
				}
				// case 2
				else if (prevqdc.card==card2 && qdc.card==card1) {
					veto.card[0]=vd[0]->GetScalerID();
					veto.card[1]=qdc.card;
					veto.card[2]=prevqdc.card;
					for (int k = 0; k<16; k++) {
						veto.QDC[k]=qdc.QDC[k];
						veto.QDC[16+k]=prevqdc.QDC[k];
						veto.IsUnderThreshold[k]=qdc.IsUnderThreshold[k];
						veto.IsUnderThreshold[16+k]=prevqdc.IsUnderThreshold[k];
						veto.IsOverflow[k]=qdc.IsOverflow[k];
						veto.IsOverflow[16+k]=prevqdc.IsOverflow[k];
					}
				}
				else if (prevqdc.card==-1) { cout << "Previous Card was 0, EventMatch: " << EventMatch << endl;  continue; }
				else { printf("Failed to match!  Run: %i  VetoTree entry: %lli  Card:%i  Prev.Card:%i  Breaking at %.0f%% through file.\n",run,i,qdc.card,prevqdc.card,((Float_t)i/nentries)*100); break; }

				// check VEvent data
				//printf("run:%i  vEnt:%i  sEnt:%i  card0:%i  card1:%i  card2:%i  sTime:%.5f  sTimeBad:%i\n",
				//	veto.run,veto.vEnt,veto.sEnt,veto.card[0],veto.card[1],veto.card[2],veto.sTime,veto.sTimeBad);
				//cout << "QDC: "; for (int k = 0; k<numPanels; ++k) { cout << veto.QDC[k] << "  "; } cout << endl;
				//cout << "UTh: "; for (int k = 0; k<numPanels; ++k) { cout << veto.IsUnderThreshold[k] << "  "; } cout << endl;
				//cout << "Ovr: "; for (int k = 0; k<numPanels; ++k) { cout << veto.IsOverflow[k] << "  "; } cout << endl;

				//=====================BEGIN ACTUAL GODDAMMED ANALYSIS=================

				if (veto.sTimeBad) { 
					BadTSInFile++;
					TotalCorruptionInTime->Fill(veto.eTime);
					if (PlotCorruptedEntries) CorruptionInTime->Fill(veto.eTime);
				}

				// multiplicity of panels above threshold
		        for (int k=0; k<numPanels; k++) { if (!veto.IsUnderThreshold[k]) numPanelsHit++; }
		        TotalMultiplicity->Fill(numPanelsHit);			
		        if (PlotMultiplicity) OneRunMultiplicity->Fill(numPanelsHit);



				//=====================END ACTUAL GODDAMMED ANALYSIS===================

			} // end EventMatch condition

			// Save qdc into prevqdc before getting next VetoTree entry.
			prevqdc = qdc;
			EventMatch = false;
			numPanelsHit=0;

		}	// End loop over VetoTree entries.

		// === END OF FILE Output & Plotting ===
		
		corruption = ((Float_t)BadTSInFile/nentries)*100;
		printf(" Corrupted scaler entries: %i of %lli, %.3f %%.\n",BadTSInFile,nentries,corruption);
		if(run>45000000) SCorruption->SetPoint(filesScanned,run-45000000,corruption);
		else SCorruption->SetPoint(filesScanned,run,corruption);

		if (PlotCorruptedEntries) {
			char outfile1[200];	
			sprintf(outfile1,"CorruptionInTime_Run%i",run);
			CorruptionInTime->Write(outfile1,TObject::kOverwrite);
		}
		if (PlotMultiplicity) {
			char outfile2[200];
			sprintf(outfile2,"Multiplicity_Run%i",run);
			OneRunMultiplicity->Write(outfile2,TObject::kOverwrite);
		}

		// ==========================

		delete VetoTree;
		delete MGTree;
		filesScanned++;
	} // End loop over files.

	

	// === END OF SCAN Output & Plotting ===
	printf("Finished loop over files.\n");

	TCanvas *c1 = new TCanvas("c1", "Bob Ross's Canvas",600,600);
	c1->SetGrid();
	SCorruption->SetMarkerColor(4);
	SCorruption->SetMarkerStyle(21);
	SCorruption->SetMarkerSize(0.5);
	SCorruption->SetTitle("Corruption in scaler card");
	SCorruption->GetXaxis()->SetTitle("Run");
	SCorruption->GetYaxis()->SetTitle("% corrupted events");
	SCorruption->Draw("ALP");	
	SCorruption->Write("ScalerCorruption",TObject::kOverwrite);

	TotalCorruptionInTime->Write("TotalCorruptionInTime",TObject::kOverwrite);

	TotalMultiplicity->Write("TotalMultiplicity",TObject::kOverwrite);
	// ==========================

	RootFile->Close();
	cout << "Wrote ROOT file." << endl;
}
