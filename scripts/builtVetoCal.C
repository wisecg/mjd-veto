/*
	builtVetoCal.C
	Clint Wiseman, USC/Majorana
	July 2015.
	
 => This code can be run on PDSF or locally.  It takes a .txt file of run numbers
	as an input argument, and uses the name of the text file to generate output
	in the folder ./output

 => Recommended: When scanning a new input file of run numbers on PDSF, run CheckFiles.C 
 	to make sure files exist and have not been blinded.  
 	This code has a tendency to quit unexpectedly when it encounters a "bad" file.
 
 => builtVetoCal is an extension of builtVeto, with the primary goal of finding
 	the peaks in the QDC spectrum from LED flashers embedded in the veto panels.


	Usage:
	CINT: root[0] .X builtVeto.C ("Filename_without_extension")  <--- NO .TXT extension.
	bash: root -b -q -l builtVeto.C ("Filename_without_extension")
	Compiled mode: 
	.X builtVeto.C++ ("Filename_without_extension")
	(must also comment in the #includes and the int main function)
*/

/*
// Compiled mode fails locally, but seems to work on PDSF.
// Don't forget to comment in/out the #ifndef bit at the end of this file.

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TPad.h>
#include <TVirtualPad.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>

// MUST load MDGOMJ classes in ROOT before calling these headers
#include "MJTRun.hh"
#include "MJTVetoData.hh"
#include "MGTBasicEvent.hh"

using namespace std;
*/
const int numPanels = 32;

// Structure to hold a complete matched veto event 
struct VEvent {
    Int_t run;				// run number
    Int_t vEnt; 			// qdc eventCount (i.e. entry)
    Int_t sEnt; 			// scaler entry
    Int_t card[3];		// 0: scaler ID, 1: qdc card 1 (panels 1-16), 2: qdc card 2 (panels 17-32)
    Float_t sTime;		// scaler time
    Bool_t sTimeBad;		// bad scaler entry
    Float_t lTime;		// LED time
    Float_t eTime;		// entry time
    Int_t QDC[numPanels];	// store qdc values for both cards in the same array 
    Int_t IsUnderThreshold[numPanels];
    Int_t IsOverflow[numPanels];

    // constructor
    VEvent(){
    	run=0; vEnt=0; sEnt=0;
    	card[0]=0;card[1]=0;card[2]=0;
    	sTime=0; sTimeBad=0; lTime=0; eTime=0;
    	for(int i=0;i<numPanels;i++) {
    		QDC[i]=0; IsUnderThreshold[i]=0; IsOverflow[i]=0;
    	}
    }
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

    // constructor
    QEvent(){
    	runNumber=0; crate=0; card=0; EventCount=0;
    	for(int i=0;i<16;i++){
    		QDC[i]=0; IsUnderThreshold[i]=0; IsOverflow[i]=0;	
    	}
    }
};

// global pointers for qdc histograms.
TH1F *hRawQDC[numPanels];  
TH1F *hCutQDC[numPanels];
TH1F *hThreshQDC[numPanels];
TH1F *hLEDCutQDC[numPanels];

void builtVetoCal(string Input = ""){

	int mode = 1;		// switch: 0 for local files, 1 for pdsf files
	UInt_t card1 = 13;	// 11: prototype, 13: module 1
	UInt_t card2 = 18;
	Bool_t useThresh = true; // if true, also enables fitting LED peaks

	// "low" qdc threshold values
	//Int_t thresh[numPanels] = {123,115,95,93,152,115,105,103,119,91,108,103,94,107,95,167,
	//	53,150,89,120,65,85,132,62,130,101,80,108,145,164,119,82};

	// "high" qdc threshold values
	Int_t thresh[numPanels] = {124,117,96,93,155,115,112,105,120,91,109,108,95,112,96,168,
		63,157,100,127,72,100,140,65,145,125,82,112,151,168,122,94};

	// Input a list of run numbers
	Char_t InputName[200];
	string def = "builtVeto_DebugList"; // default
	if (Input == "") strcpy(InputName,def.c_str());
	else strcpy(InputName,Input.c_str());
	Char_t InputFile[200];
	sprintf(InputFile,"%s.txt",InputName);
	ifstream InputList;
	InputList.open(InputFile);
	Char_t TheFile[200];

	// Set up output file(s)
	Char_t OutputFile[200];
	sprintf(OutputFile,"./output/%s.root",InputName);
	TFile *RootFile = new TFile(OutputFile, "RECREATE");	
  	TH1::AddDirectory(kFALSE); // Global flag: "When a (root) file is closed, all histograms in memory associated with this file are automatically deleted."

	ofstream calib;
	Char_t CalibFile[200];
	sprintf(CalibFile,"./output/%s_CalibrationTable.txt",InputName);
	calib.open(CalibFile);


	//=== Global counters / variables / plots ===

		Int_t run = 0;
		Float_t duration = 0;
		Float_t durationTotal = 0;
		
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
	 	Bool_t PlotCorruptedEntries = false; // flag for plotting corrupted entries in time for EACH RUN

	 	TH1D *TotalMultiplicity = new TH1D("TotalMultiplicity","Events over threshold",32,0,32);
	 	TotalMultiplicity->GetXaxis()->SetTitle("number of panels hit");
	 	Bool_t PlotMultiplicity = false;	// flag to plot multiplicity for EACH RUN

	 	const Int_t nqdc_bins=1400;  // this gives 3 qdc / bin
		const Float_t ll_qdc=0.;
		const Float_t ul_qdc=4200.;
		Char_t hname[50];
		for (Int_t i=0; i<numPanels; i++){
			sprintf(hname,"hRawQDC%d",i);
			hRawQDC[i] = new TH1F(hname,hname,nqdc_bins,ll_qdc,ul_qdc);
			sprintf(hname,"hCutQDC%d",i);
			hCutQDC[i] = new TH1F(hname,hname,nqdc_bins,ll_qdc,ul_qdc);
			sprintf(hname,"hThreshQDC%d",i);
			hThreshQDC[i] = new TH1F(hname,hname,500,ll_qdc,500);
			sprintf(hname,"hLEDCutQDC%d",i);
			hLEDCutQDC[i] = new TH1F(hname,hname,500,ll_qdc,500);
		}

		Long64_t CountsBelowThresh[numPanels] = {0};
		Long64_t TotalCounts[numPanels] = {0};		

		
	//=== End ===



	// Loop over files
	while(!InputList.eof()){

		// initialize 
		InputList >> run;
		if (mode==0) sprintf(TheFile,"~/dev/datasets/muFinder/OR_run%i.root",run);
		else if (mode==1) sprintf(TheFile,"/global/project/projectdirs/majorana/data/mjd/surfmjd/data/built/P3JDY/OR_run%u.root",run); 
		TChain *VetoTree = new TChain("VetoTree");
		VetoTree->AddFile(TheFile);
		TChain *MGTree = new TChain("MGTree");
		MGTree->AddFile(TheFile);
		MJTRun *MyRun = new MJTRun();
		MGTree->SetBranchAddress("run",&MyRun);
		Long64_t nentries = VetoTree->GetEntries();
		MGTBasicEvent *b = new MGTBasicEvent; 
		//MGTBasicEvent b;
		VetoTree->SetBranchAddress("vetoEvent",&b);
		const int vd_size = 16;	
		MJTVetoData *vd[vd_size];
		VetoTree->GetEntry(0);
		for (int i=0; i<vd_size; i++) { vd[i] = dynamic_cast<MJTVetoData*>(b->GetDetectorData()->At(i)); }
        MGTree->GetEntry(0);
        duration = MyRun->GetStopTime() - MyRun->GetStartTime();
        if (duration < 0) {
        	printf("\nRun %i has duration %.0f, skipping file!\n\n",run,duration);
        	continue;
        }
        durationTotal += duration;

    	//=== Single-file counters / variables / plots

			QEvent qdc;
			QEvent prevqdc;
				prevqdc.EventCount=-1;
			VEvent veto;
			Bool_t EventMatch = false;

			Int_t BadTSInFile = 0;
			Float_t corruption = 0;
			Int_t numPanelsHit = 0;

			// only write these if their bools are set = true
			TH1D *CorruptionInTime = new TH1D("CorruptionInTime","corrupted entries during run (entry method)",(Int_t)duration/5,0,(Int_t)duration);
			CorruptionInTime->GetXaxis()->SetTitle("time (5 sec / bin)");
			TH1D *OneRunMultiplicity = new TH1D("multiplicity","multiplicity of veto entries",32,0,32);
			OneRunMultiplicity->GetXaxis()->SetTitle("number of panels hit");
		
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
				else if ((int)prevqdc.card==-1) { cout << "Previous Card was 0, EventMatch: " << EventMatch << endl;  continue; }
				else { printf("Failed to match!  Run: %i  VetoTree entry: %i  Card:%i  Prev.Card:%i  Breaking at %.0f%% through file.\n",run,i,qdc.card,prevqdc.card,((Float_t)i/nentries)*100); break; }

				// check VEvent data
				//printf("run:%i  vEnt:%i  sEnt:%i  card0:%i  card1:%i  card2:%i  sTime:%.5f  sTimeBad:%i\n",
				//	veto.run,veto.vEnt,veto.sEnt,veto.card[0],veto.card[1],veto.card[2],veto.sTime,veto.sTimeBad);
				//cout << "QDC: "; for (int k = 0; k<numPanels; ++k) { cout << veto.QDC[k] << "  "; } cout << endl;
				//cout << "UTh: "; for (int k = 0; k<numPanels; ++k) { cout << veto.IsUnderThreshold[k] << "  "; } cout << endl;
				//cout << "Ovr: "; for (int k = 0; k<numPanels; ++k) { cout << veto.IsOverflow[k] << "  "; } cout << endl;

				//=====================BEGIN ACTUAL GODDAMMED ANALYSIS=================

				// scaler corruption
				if (veto.sTimeBad) { 
					BadTSInFile++;
					TotalCorruptionInTime->Fill(veto.eTime);
					if (PlotCorruptedEntries) CorruptionInTime->Fill(veto.eTime);
				}

				// loop over panels
				for (int k=0; k<numPanels; k++) {
					
					// test lowered panel-by-panel thresholds
					if (veto.QDC[k]<thresh[k]) CountsBelowThresh[k]++;
					TotalCounts[k]++;

					// plot qdc entries above threshold
					hRawQDC[k]->Fill(veto.QDC[k]);
					if (veto.QDC[k]<500) hThreshQDC[k]->Fill(veto.QDC[k]);
					if (useThresh && veto.QDC[k]>thresh[k]) hCutQDC[k]->Fill(veto.QDC[k]);

					// count multiplicity
					if (useThresh) { if (veto.QDC[k]>thresh[k]) numPanelsHit++; }
					else { if (!veto.IsUnderThreshold[k]) numPanelsHit++; }
				}
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

	// estimate reduction in data
	if (useThresh) {
		Long64_t total=0, below=0;
		double reduction;
		for (int i=0; i<numPanels; i++) {
			below += CountsBelowThresh[i];
			total += TotalCounts[i];
		}
		reduction = ((double)below/total)*100;
		double rateBefore = (double)total/durationTotal; 
		double rateAfter = (double)(total-below)/durationTotal;
		printf("Counts below thresh make up %.2f%% of total entries scanned, over %.2f seconds.\n",reduction,durationTotal);
		printf("This is %lli of %lli events.  Rate reduction would be from %.2f Hz to %.2f Hz.\n",below,total,rateBefore,rateAfter);
	}

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

	// QDC Plots & Calibration Table:
	// gaus: A gaussian with 3 parameters: f(x) = p0*exp(-0.5*((x-p1)/p2)^2)).
	TF1 *fits[numPanels];
	TCanvas *vcan0 = new TCanvas("vcan0","cut & fitted veto QDC, panels 1-32",0,0,800,600);
	vcan0->Divide(8,4,0,0);
	TCanvas *vcan1 = new TCanvas("vcan1","veto QDC thresholds, panels 1-32",0,0,800,600);
	vcan1->Divide(8,4,0,0);	
	Char_t buffer[2000];
	printf("\n Calibration Table:\n  Panel / Mean,error / Sigma,error / Chi-square/NDF (~1?) / LED Peak Pos.\n");
	for (Int_t i=0; i<numPanels; i++){
		
		vcan0->cd(i+1);
		TVirtualPad *vpad0 = vcan0->cd(i+1); vpad0->SetLogy();
		hThreshQDC[i]->Write(); // write the low-QDC part of the spectrum separately
		if (useThresh) {
			hCutQDC[i]->Write();    // write the cut QDC
			hCutQDC[i]->Fit("gaus","q");
			hCutQDC[i]->Draw();
			fits[i] = hCutQDC[i]->GetFunction("gaus");
			if (hCutQDC[i]->GetEntries() > 0) {
				Float_t NDF = (Float_t)fits[i]->GetNDF();
				if (fits[i]->GetNDF() == 0) NDF = 0.0000001;
				sprintf(buffer,"%i  %.1f  %.1f  %.1f  %.1f  %.1f",i,
					fits[i]->GetParameter(1),fits[i]->GetParError(1),		// mean (LED center)
					fits[i]->GetParameter(2),fits[i]->GetParError(2),		// sigma
					fits[i]->GetChisquare()/NDF); // X^2/NDF (closer to 1, better fit.)
				cout << "  " << buffer << endl;
				calib << buffer << endl;
			}
		}
		else {
			hRawQDC[i]->Write();		// write the raw QDC without fitting
		}	
    	
   	// plot low-QDC range separately
  		vcan1->cd(i+1);
  		TVirtualPad *vpad1 = vcan1->cd(i+1); vpad1->SetLogy();
  		hThreshQDC[i]->Draw();
	}

	// Output canvasses of interest.
	Char_t OutputName[200];
	
	sprintf(OutputName,"./output/%s_VetoQDC.C",InputName);
	vcan0->Print(OutputName);

	sprintf(OutputName,"./output/%s_ThreshVetoQDC.C",InputName);
	vcan1->Print(OutputName);
	
	// ==========================

	calib.close();
	RootFile->Close();
	cout << "Wrote ROOT file." << endl;
}

/*
// this should be commented in if the code is run in compiled mode
// should be able to take multiple inputs, but I can't figure out the ROOT syntax:
// .x builtVetoCal.C++ ("M1BG_debug" "M1BG_Range5")  FAILS.
#ifndef __CINT__
int main(int argc, const char* argv[]){
	size_t size;
	if (argc>0) {
		for(int i = 0; i < argc; i++) {
			printf( "Input arg %d: %s\n", i, argv[i] );
			size = sizeof argv[i]/sizeof(size_t);
			string Input(argv[i],size);
			builtVetoCal(Input);
		}
	}
	else {
		string Input = "";
		builtVetoCal(Input);
	}
}
#endif
*/