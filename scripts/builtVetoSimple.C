/*
	builtVetoSimple.C
	Clint Wiseman, USC/Majorana
	June 2015.
	
 => This code is run on local data files.
 => The user specifies the input text file of runs inside the code.
 => Makes a TGraph of corruption in the scaler as an example.
 
 => The main purpose of this code was to develop the "sorting" of separate 
    QDC cards into one consistent structure, VEvent.  This helps to separate processing
    from analysis in the code.

 => The veto data should always come out from ORCA in groups of 3: Scaler, QDC1, QDC2.
 	Once it goes through the MJOR event builder, the scaler data is merged into the QDC entries.
 	The built data (usually/always) contains TWO QDC entries associated with the same veto event.
    In this code, the sorting method checks to see if the "previous" QDC event is the same
    as the current one, using the variable "EventCount".
    
 => There have been difficulties with the eventCount variable.  Usually
    it goes like 0,0,1,1,2,2,3,3 etc.  The QEvent structure has in the past been used 
    to *force* matching when eventCount goes like 0,1,0,2,1,3... or similar, by saving a 
    std::list of un-matched QEvents in memory, and using a "fifo" queue for sorting.

	Usage:
	CINT: root[0] .X builtVetoSimple.C
	bash: root -b -q -l builtVetoSimple.C
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


void builtVeto(){

	// Input a list of run numbers
	Char_t InputFile[200] = "builtVeto_DebugList.txt";
	ifstream InputList;
	InputList.open(InputFile);
	Int_t run;
	Char_t TheFile[200];

	//=== Global counters / variables / plots ===
		
		// get number of files in dataset for the TGraph
		Int_t filesToScan = 0;
		Int_t filesScanned = 0;
	  	while(!InputList.eof()) { InputList >> run; filesToScan++; }
	  	cout << "Scanning " << filesToScan << " files." << endl;	  	
	  	InputList.close();
	  	InputList.open(InputFile);
	  	run=0;

	 	TGraph *SCorruption = new TGraph(filesToScan);
		
	//=== End ===


	// Loop over files
	while(!InputList.eof()){

		//=== Single-file counters / variables / plots
			QEvent qdc = {0};
			QEvent prevqdc = {0};
				prevqdc.EventCount=-1;
			VEvent veto = {0};
			Bool_t EventMatch = false;
		//=== End ===

		InputList >> run;
		sprintf(TheFile,"~/dev/datasets/builtVeto/OR_run%i.root",run);

		TFile *f = new TFile(TheFile);
		TTree *VetoTree = (TTree*)f->Get("VetoTree");

		MGTBasicEvent *b = 0;
		VetoTree->SetBranchAddress("vetoEvent",&b);
		Long64_t nentries = VetoTree->GetEntries();
		cout << "Found " << nentries << " entries." << endl;

		const int vd_size = 16;	
		MJTVetoData *vd[vd_size];
		VetoTree->GetEntry(0);
		for (int i=0; i<vd_size; i++) { vd[i] = dynamic_cast<MJTVetoData*>(b->GetDetectorData()->At(i)); }
		cout << "Initialized MJTVetoData" << endl;

		// Loop over VetoTree entries
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
				veto.eTime = (Float_t)i/nentries;	// needs the duration: pull in fStopTime-fStartTime from MGTEventTree.
				// case 1
				if (prevqdc.card==11 && qdc.card==18) {
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
				else if (prevqdc.card==18 && qdc.card==11) {
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
				else { printf("Failed to match!  VetoTree entry: %lli  Card:%i  Prev.Card:%i  \n",i,qdc.card,prevqdc.card); break; }

				// check VEvent data
				printf("run:%i  vEnt:%i  sEnt:%i  card0:%i  card1:%i  card2:%i  sTime:%.5f  sTimeBad:%i\n",
					veto.run,veto.vEnt,veto.sEnt,veto.card[0],veto.card[1],veto.card[2],veto.sTime,veto.sTimeBad);
				//cout << "QDC: "; for (int k = 0; k<numPanels; ++k) { cout << veto.QDC[k] << "  "; } cout << endl;
				//cout << "UTh: "; for (int k = 0; k<numPanels; ++k) { cout << veto.IsUnderThreshold[k] << "  "; } cout << endl;
				//cout << "Ovr: "; for (int k = 0; k<numPanels; ++k) { cout << veto.IsOverflow[k] << "  "; } cout << endl;

				//=====================BEGIN ANALYSIS=================



				//=====================END ANALYSIS===================

			} // end EventMatch condition

			// Save qdc into prevqdc before getting next VetoTree entry.
			prevqdc = qdc;
			EventMatch = false;

		}	// End loop over VetoTree entries.

	filesScanned++;
	} // End loop over files.

	// === END OF SCAN Output & Plotting ===
	
	// ==========================

	// Close input/output files.
}