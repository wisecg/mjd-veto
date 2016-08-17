/*
	runhealth.C
	Andrew Lopez, UTK/Majorana
	January 2016

 => This code can be run on PDSF.  It takes a .txt file as an input argument,
	and uses the name of the text file to generate output
	
=> The input .txt file MUST be called runhealth_list.txt
	The input file can have run numbers from various directories, an example of the format is below
	
		@DirectoryName1
		10000
		10001
		10002
		@DirectoryName2
		20000
		20001
		20002
		@DirectoryName3
		30000
		30001
		30002
		!DONOTDELETE
	
	Usage:
	CINT: root[0] .X runhealth.C ()
	bash: root -b -q -l runhealth.C ()
*/

#ifndef __CINT__
#include <iostream>
#include <fstream>
#include <sstream>
#include <fstream>
#include <string>
#include <stdlib.h>

#include "TFile.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include "TF1.h"

#include "GATDataSet.hh"
#include "MJVetoEvent.hh"
#include "MJTRun.hh"
#include "MJTVetoData.hh"
#include "MGTBasicEvent.hh"
#endif

using namespace std;
const int numPanels = 32;

// ==================================================
// Processing Functions
// ==================================================
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
	string InputLine;
	InputList.open(arg.c_str());
	if(!InputList.good()) {
    	cout << "Couldn't open " << arg << " !" << endl;
    	return 0;
    }

	int filesToScan = 0;
  	while (true){
		int run = 0;
    	InputList >> InputLine;
		if (InputList.eof()) break;
  		if (InputLine[0] == '0' || InputLine[0] == '1' || InputLine[0] == '2' || InputLine[0] == '3' || InputLine[0] == '4' || InputLine[0] == '5' || InputLine[0] == '6' || InputLine[0] == '7' || InputLine[0] == '8' || InputLine[0] == '9'){
			run = atoi(InputLine.c_str());
			filesToScan++;
		}
  	}
  	cout << "Scanning " << filesToScan << " files." << endl;	  	
  	InputList.close();
  	return filesToScan;
}

//Begin program
void runhealth(){
	
	// Set up output file
	string file = "runhealth_list.txt";
	string filename = file;
	filename.erase(filename.find_last_of("."), string::npos); 
	char OutputFile[200];
	sprintf(OutputFile,"%s.root",filename.c_str());
	TFile *RootFile = new TFile(OutputFile, "RECREATE");	
  	TH1::AddDirectory(kFALSE); // Global flag: "When a (root) file is closed, all histograms in memory associated with this file are automatically deleted."
	int filesToScan = GetNumFiles(file);
	

	//define global variables/arrays/plots
  	int run = 0;			// run number
	const int ledcut = 20; //define led cut as > 20 panels hit
	int filesScanned = 0;	// counter	
	Int_t ledcount = 0;
	Int_t threshold[numPanels] = {0}; //recalculated led threshold for runs investigated
	Float_t totalduration = 0;
	Int_t totalnentries = 0;
	Int_t totalledcount = 0;
	Int_t totalledcountPanel[numPanels] = {0}; 
	Float_t ledqdcsquaredsum[numPanels] = {0}; 
	Float_t ledqdcsum[numPanels] = {0}; 	 

	//mjvbits event counters
	Int_t kmcevcount = 0; //kMissingChannels counter
	Int_t kecevcount = 0; //kExtraChannels counter
	Int_t ksoevcount = 0; //kScalerOnly counter
	Int_t kbtsevcount = 0; //kBadTimeStamp counter
	Int_t koosevcount = 0; //kQDCOutOfSequence counter
	Int_t kdcevcount = 0; //kDuplicateChannels counter
	Int_t khwcmevcount = 0; //kHWCountMismatch counter
	Int_t kLEDsOffevcount = 0;	//counts how many events are NOT LEDs
	Int_t ksoerrorevcount = 0; //counts how many events with missing QDC fail to trigger kScalerOnly
	
	//mjvbits file counters
	Int_t kmcfilecount = 0; //kMissingChannels file counter
	Int_t kecfilecount = 0; //kExtraChannels file counter
	Int_t ksofilecount = 0; //kScalerOnly file counter
	Int_t kbtsfilecount = 0; //kBadTimeStamp file counter
	Int_t koosfilecount = 0; //kQDCOutOfSequence file counter
	Int_t kdcfilecount = 0; //kDuplicateChannels file counter
	Int_t khwcmfilecount = 0; //kHWCountMismatch file counter
	Int_t kLEDsOffFilecount = 0; //counts how many files have ZERO LEDs
	Int_t ksoerrorfilecount = 0; //counts how many files with missing QDC fail to trigger kScalerOnly 
	
	//led (low) qdc threshold values from findThresh.C
	Int_t ledthresh[numPanels] = {545, 442, 642, 478, 818, 525, 508, 959, 729, 674, 626, 553, 458, 438, 670, 784, 348, 412, 486, 443, 613, 600, 382, 480, 632, 388, 562, 444, 480, 596, 661, 515};
	
	
	
	//prepare histograms and graphs
	TH1D *TotalMultiplicity = new TH1D("TotalMultiplicity","Events over threshold",numPanels,0,numPanels);
	TotalMultiplicity->GetXaxis()->SetTitle("number of panels hit");
	TH1F *hRawQDC[numPanels];
	const int nqdc_bins=4200;
	const int ll_qdc = 0;
	const int ul_qdc=4200;
	Char_t hname[50];
	for (Int_t i=0; i<numPanels; i++){
		sprintf(hname,"hRawQDC%d",i);
		hRawQDC[i] = new TH1F(hname,hname,nqdc_bins,ll_qdc,ul_qdc);
		hRawQDC[i]->GetXaxis()->SetTitle("QDC Value");
	}
	
	TGraph *gRunvFreq;
	gRunvFreq = new TGraph(filesToScan);
	TH1F *hDTFile;
	hDTFile = new TH1F[filesToScan];
	TGraph *gMultvTimeFile;
	gMultvTimeFile = new TGraph[filesToScan];
	
	
	//start file loop
	// Loop over files in dataset
	ifstream InputList;
  	InputList.open(file.c_str());
	string InputLine;
	Char_t TheFile[200];
	char run_dir[] = "empty";
	if (InputList.is_open()){
		while(true){	
			int run = 0;

			
			//------read and organize the input from the input file (directory, run #, nonsense)
			{
			InputList >> InputLine;
			if (InputList.eof()) break;

			if (InputLine.at(0) == '@'){	
				InputLine.erase(0,1);
				InputLine.copy(run_dir,6,0);
			}
			if (InputLine[0] == '0' || InputLine[0] == '1' || InputLine[0] == '2' || InputLine[0] == '3' || InputLine[0] == '4' || InputLine[0] == '5' || InputLine[0] == '6' || InputLine[0] == '7' || InputLine[0] == '8' || InputLine[0] == '9'){
				run = atoi(InputLine.c_str());
			}
			if (InputLine.at(0) == '!') ; //do nothing
			

			cout << "directory = " << run_dir << endl;
			cout << "run = " << run << endl;
			}
			//---------finished reading out input
			
			
			
			
			//begin analysis if run != 0, which means the input wasn't a directory or the eof line
			if (run != 0) {
				sprintf(TheFile,"/global/project/projectdirs/majorana/data/mjd/surfmjd/data/built/%s/OR_run%u.root",run_dir,run); 
				
				//initialize veto data
				TChain *VetoTree = new TChain("VetoTree");
				VetoTree->AddFile(TheFile);
				Long64_t nentries = VetoTree->GetEntries();
				MJTRun *VetoRun = new MJTRun();
				MGTBasicEvent *vetoEvent = new MGTBasicEvent();
				UInt_t mVeto = 0;
		
				//set branch addresses
				VetoTree->SetBranchAddress("run",&VetoRun);
				VetoTree->SetBranchAddress("mVeto",&mVeto);
				VetoTree->SetBranchAddress("vetoEvent",&vetoEvent);
				uint32_t vBits = 0;
				VetoTree->SetBranchAddress("vetoBits",&vBits);
				
							
				// Unsigned int from MGTypes.hh -- kData=0, kTest=1, kCalibration=2, kMC=3, kUndefined=4																	  
				printf("Run Type: %u\n",VetoRun->GetRunType()); 
				
				//single file counters/ variables/ plots
				sprintf(hname,"hDTFile%d", run);
				hDTFile[filesScanned] = new TH1F(hname,hname,100,0,10);
				sprintf(hname,"gMultvTimeFile%d", run);
				gMultvTimeFile[filesScanned] = new TGraph(nentries);
				gMultvTimeFile[filesScanned].SetName(hname);

				bool IsEmpty = false;
				long duration = 0;
				long start = 0;
				long stop = 0;
				Int_t ledcount = 0;
				Float_t Tevent = 0;
				Float_t Tevent_prev = 0;
				Float_t deltaT = 0;

				
				//veto error bit file bits/counters
				bool kmcfilebit = false;
				bool kecfilebit = false;
				bool ksofilebit = false;
				bool kbtsfilebit = false;
				bool koosfilebit = false;
				bool kdcfilebit = false;
				bool khwcmfilebit = false;
				bool ksoerrorfilebit = false;
				
				GATDataSet ds(run);
				start = GetStartUnixTime(ds);
				stop = GetStopUnixTime(ds);
				duration = stop - start;
				totalduration += duration;
				totalnentries += nentries;

			
				// Loop over VetoTree entries
				printf("Now scanning run %i: %lli entries, %.2f sec.  \n",run,nentries,duration);
				if (nentries == 0) IsEmpty = true;
				for (int z = 0; z < nentries; z++) {
					VetoTree->GetEntry(z);
					
					//single entry variables
					bool isLED = false;
					Int_t lednumPanelsHit = 0;
					Float_t time = vetoEvent->GetTime()/1E9; //event time in seconds
					bool kmcbit = false;
					bool kecbit = false;
					bool ksobit = false;
					bool kbtsbit = false;
					bool koosbit = false;
					bool kdcbit = false;
					bool khwcmbit = false;
					bool kLEDsOffbit = false;
					bool ksoerrorbit = false;
					
					
//[A. Begin] access mjtVetoData and sort into arrays
				{
					// Access the MJTVetoData objects "vd"
					MJTVetoData *vd[numPanels];
					for (int i=0; i<numPanels; i++) { vd[i] = dynamic_cast<MJTVetoData*>(vetoEvent->GetDetectorData()->At(i)); }
					
					//********************************************************************
					//sort data into arrays
					// Sort raw data into arrays and then display.
					// This may not be totally necessary, but makes hit pattern analysis easier
					//   to match to the physical veto panel locations, and Yuri's wiring diagrams.
					// Most things are cast to int's.  
					// Original types can be found in MJTVetoData.hh and MGDetectorData.hh if necessary.
					const int card1 = 13;
					const int card2 = 18;
					int Card[numPanels] = {0};
					int QDC[numPanels] = {0};
					int IsUnderThreshold[numPanels] = {0};
					int IsOverflow[numPanels] = {0};
					int ID[numPanels] = {0};
					long Index[numPanels] = {0};
					int k = 0;
					for (int j = 0; j<numPanels; j++)	{
						if (vd[j]){
							k = vd[j]->GetChannel();	// goes from 0 to 15
							if (vd[j]->GetCard() == card1) {
								Card[k] = vd[j]->GetCard();
								QDC[k] = (int)vd[j]->GetAmplitude();
								IsUnderThreshold[k] = (int)vd[j]->IsUnderThreshold();
								IsOverflow[k] = (int)vd[j]->IsOverflow();
								ID[k] = vd[j]->GetID();
								Index[k] = (Long_t)vd[j]->GetIndex();
							}
							else if (vd[j]->GetCard() == card2) {
								Card[16+k] = vd[j]->GetCard();
								QDC[16+k] = (int)vd[j]->GetAmplitude();
								IsUnderThreshold[16+k] = (int)vd[j]->IsUnderThreshold();
								IsOverflow[16+k] = (int)vd[j]->IsOverflow();	
								ID[16+k] = vd[j]->GetID();
								Index[16+k] = (Long_t)vd[j]->GetIndex();
							}
						}
					}
				}
//[A. End]
					
						

					//being led analysis
					//fill event DT histograms
					Tevent = time;
					deltaT = Tevent - Tevent_prev;
					Tevent_prev = Tevent;
					hDTFile[filesScanned].Fill(deltaT);
				


				
					//Identify LEDs
					for (int k = 0; k<numPanels; k++){
	
						hRawQDC[k]->Fill(QDC[k]);
							if (QDC[k] > ledthresh[k]) lednumPanelsHit++;
					}
					//cout << "lednumPanelsHit = " << lednumPanelsHit << endl;
					
					//fill event multiplicity vs time 
					if (time < 5000){
						gMultvTimeFile[filesScanned].SetPoint(z,time,lednumPanelsHit);
					}	
					
					TotalMultiplicity->Fill(lednumPanelsHit);
					if (lednumPanelsHit > ledcut) {
						isLED = true;
						ledcount++;
						totalledcount++;
					}
					//calc rms 
					if (isLED){
						
						for (int k = 0; k<numPanels; k++){			
							totalledcountPanel[k]++;					
							ledqdcsum[k] += QDC[k];					
							ledqdcsquaredsum[k] += QDC[k]*QDC[k];	
						}											
					
					}
					
					//vbits
					if (MJBits::GetBit(vBits, MJVetoBits::kMissingChannels)) {kmcevcount++; kmcbit = true; kmcfilebit = true;}
					if (MJBits::GetBit(vBits, MJVetoBits::kExtraChannels)) {kecevcount++; kecbit = true; kecfilebit = true;}
					if (MJBits::GetBit(vBits, MJVetoBits::kScalerOnly)) {ksoevcount++; ksobit = true; ksofilebit = true;}
					if (MJBits::GetBit(vBits, MJVetoBits::kBadTimeStamp)) {kbtsevcount++; kbtsbit = true; kbtsfilebit = true;}
					if (MJBits::GetBit(vBits, MJVetoBits::kQDCOutOfSequence)) {koosevcount++; koosbit = true; koosfilebit = true;}
					if (MJBits::GetBit(vBits, MJVetoBits::kDuplicateChannel)) {kdcevcount++; kdcbit = true; kdcfilebit = true;}
					if (MJBits::GetBit(vBits, MJVetoBits::kHWCountMismatch)) {khwcmevcount++; khwcmbit = true; khwcmfilebit = true;}
					if (!isLED) {kLEDsOffevcount++; kLEDsOffbit = true;}
					
					for (Int_t i=0; i<numPanels; i++){
						if(QDC[i] == 0 && !ksobit) {ksoerrorbit = true; ksoerrorfilebit = true;}
					}	
					if (ksoerrorbit) ksoerrorevcount++;
//[B. Begin] Print out results for each bad entry
{
/*
					//simple printout 
					if (kmcbit){ 
					printf(" New Error: Bit kMissingChannels  Count: %i  Run: %i  Entry: %li  Time: %.5f  %.2f%%.\n"
					,kmcevcount,run,z,time,((double)z/nentries)*100);} 
						
					if (kecbit) {
						printf(" New Error: Bit kExtraChannels  Count: %i  Run: %i  Entry: %li  Time: %.5f  %.2f%%.\n"
					,kecevcount,run,z,time,((double)z/nentries)*100); }		
				
					if (ksobit) {
						printf(" New Error: Bit kScalerOnly  Count: %i  Run: %i  Entry: %li  Time: %.5f  %.2f%%.\n"
					,ksoevcount,run,z,time,((double)z/nentries)*100); }
						
					if (kbtsbit) {
						printf(" New Error: Bit kBadTimeStamp  Count: %i  Run: %i  Entry: %li  Time: %.5f  %.2f%%.\n"
					,kbtsevcount,run,z,time,((double)z/nentries)*100); }

					if (koosbit) {
						printf(" New Error: Bit kQDCOutOfSequence  Count: %i  Run: %i  Entry: %li  Time: %.5f  %.2f%%.\n"
					,koosevcount,run,z,time,((double)z/nentries)*100); }

					if (kdcbit) {
						printf(" New Error: Bit kDuplicateChannels  Count: %i  Run: %i  Entry: %li  Time: %.5f  %.2f%%.\n"
					,kdcevcount,run,z,time,((double)z/nentries)*100); }

					if (khwcmbit) {
						printf(" New Error: Bit kHWCountMismatch  Count: %i  Run: %i  Entry: %li  Time: %.5f  %.2f%%.\n"
					,khwcmevcount,run,z,time,((double)z/nentries)*100); } 
*/		
/*				
					// detailed printout
					if (kmcbit){
						printf(" New Error: Bit kMissingChannels  Count: %i  Run: %i  Entry: %li  Time: %.5f  %.2f%%.\n"
						,kmcevcount,run,z,time,((double)z/nentries)*100); 
					
						// Info for fIndex and EventCount mismatches.
						printf("   Channels: %lu  scaler_Ind: %li  qdc1_Ind: %li  qdc2_Ind: %li  s-QDC1: %li  s-qdc2: %li  SEC: %li  QEC: %li\n\n"
						,vetoEvent->GetNDetectorData(),vd[0]->GetScalerIndex(),Index[0],Index[31]
						,vd[0]->GetScalerIndex() - Index[0], vd[0]->GetScalerIndex() - Index[31]
						,vd[0]->GetScalerCount(),vd[0]->GetEventCount());
					}

					if (kecbit){
						printf(" New Error: Bit kExtraChannels  Count: %i  Run: %i  Entry: %li  Time: %.5f  %.2f%%.\n"
						,kecevcount,run,z,time,((double)z/nentries)*100); 
					
						// Info for fIndex and EventCount mismatches.
						printf("   Channels: %lu  scaler_Ind: %li  qdc1_Ind: %li  qdc2_Ind: %li  s-QDC1: %li  s-qdc2: %li  SEC: %li  QEC: %li\n\n"
						,vetoEvent->GetNDetectorData(),vd[0]->GetScalerIndex(),Index[0],Index[31]
						,vd[0]->GetScalerIndex() - Index[0], vd[0]->GetScalerIndex() - Index[31]
						,vd[0]->GetScalerCount(),vd[0]->GetEventCount());
					}
					
					if (ksobit){
						printf(" New Error: Bit kScalerOnly  Count: %i  Run: %i  Entry: %li  Time: %.5f  %.2f%%.\n"
						,ksoevcount,run,z,time,((double)z/nentries)*100); 
					
						// Info for fIndex and EventCount mismatches.
						printf("   Channels: %lu  scaler_Ind: %li  qdc1_Ind: %li  qdc2_Ind: %li  s-QDC1: %li  s-qdc2: %li  SEC: %li  QEC: %li\n\n"
						,vetoEvent->GetNDetectorData(),vd[0]->GetScalerIndex(),Index[0],Index[31]
						,vd[0]->GetScalerIndex() - Index[0], vd[0]->GetScalerIndex() - Index[31]
						,vd[0]->GetScalerCount(),vd[0]->GetEventCount());
					}
					
					if (kbtsbit){
						printf(" New Error: Bit kBadTimeStamp  Count: %i  Run: %i  Entry: %li  Time: %.5f  %.2f%%.\n"
						,kbtsevcount,run,z,time,((double)z/nentries)*100); 
					
						// Info for fIndex and EventCount mismatches.
						printf("   Channels: %lu  scaler_Ind: %li  qdc1_Ind: %li  qdc2_Ind: %li  s-QDC1: %li  s-qdc2: %li  SEC: %li  QEC: %li\n\n"
						,vetoEvent->GetNDetectorData(),vd[0]->GetScalerIndex(),Index[0],Index[31]
						,vd[0]->GetScalerIndex() - Index[0], vd[0]->GetScalerIndex() - Index[31]
						,vd[0]->GetScalerCount(),vd[0]->GetEventCount());
					}
					
					if (koosbit){
						printf(" New Error: Bit kQDCOutOfSequence  Count: %i  Run: %i  Entry: %li  Time: %.5f  %.2f%%.\n"
						,koosevcount,run,z,time,((double)z/nentries)*100); 
					
						// Info for fIndex and EventCount mismatches.
						printf("   Channels: %lu  scaler_Ind: %li  qdc1_Ind: %li  qdc2_Ind: %li  s-QDC1: %li  s-qdc2: %li  SEC: %li  QEC: %li\n\n"
						,vetoEvent->GetNDetectorData(),vd[0]->GetScalerIndex(),Index[0],Index[31]
						,vd[0]->GetScalerIndex() - Index[0], vd[0]->GetScalerIndex() - Index[31]
						,vd[0]->GetScalerCount(),vd[0]->GetEventCount());
					}
					
					if (kdcbit){
						printf(" New Error: Bit kDuplicateChannels  Count: %i  Run: %i  Entry: %li  Time: %.5f  %.2f%%.\n"
						,kdcevcount,run,z,time,((double)z/nentries)*100); 
					
						// Info for fIndex and EventCount mismatches.
						printf("   Channels: %lu  scaler_Ind: %li  qdc1_Ind: %li  qdc2_Ind: %li  s-QDC1: %li  s-qdc2: %li  SEC: %li  QEC: %li\n\n"
						,vetoEvent->GetNDetectorData(),vd[0]->GetScalerIndex(),Index[0],Index[31]
						,vd[0]->GetScalerIndex() - Index[0], vd[0]->GetScalerIndex() - Index[31]
						,vd[0]->GetScalerCount(),vd[0]->GetEventCount());
					}
					
					if (khwcmbit){
						printf(" New Error: Bit kHWCountMismatch  Count: %i  Run: %i  Entry: %li  Time: %.5f  %.2f%%.\n"
						,khwcmevcount,run,z,time,((double)z/nentries)*100); 
					
						// Info for fIndex and EventCount mismatches.
						printf("   Channels: %lu  scaler_Ind: %li  qdc1_Ind: %li  qdc2_Ind: %li  s-QDC1: %li  s-qdc2: %li  SEC: %li  QEC: %li\n\n"
						,vetoEvent->GetNDetectorData(),vd[0]->GetScalerIndex(),Index[0],Index[31]
						,vd[0]->GetScalerIndex() - Index[0], vd[0]->GetScalerIndex() - Index[31]
						,vd[0]->GetScalerCount(),vd[0]->GetEventCount());
					}

					if (ksoerrorbit){
						printf(" New Error: Bit No QDC entries but bit kScalerOnly not thrown  Count: %i  Run: %i  Entry: %li  Time: %.5f  %.2f%%.\n"
						,ksoerrorevcount,run,z,time,((double)z/nentries)*100); 
					
						// Info for fIndex and EventCount mismatches.
						printf("  kScalerOnly Bit: %i  QDC[0]: %f  QDC[1]: %f  QDC[2]: %f  QDC[3]: %f  QDC[4]: %f  QDC[5]: %f  QDC[6]: %f  QDC[7]: %f  QDC[8]: %f  QDC[9]: %f  QDC[10]: %f  QDC[11]: %f  QDC[12]: %f  QDC[13]: %f  QDC[14]: %f  QDC[15]: %f  QDC[16]: %f  QDC[17]: %f  QDC[18]: %f  QDC[19]: %f  QDC[20]: %f  QDC[21]: %f  QDC[22]: %f  QDC[23]: %f  QDC[24]: %f  QDC[25]: %f  QDC[26]: %f  QDC[27]: %f  QDC[28]: %f  QDC[29]: %f  QDC[30]: %f  QDC[31]: %f\n\n"
						,ksobit, QDC[0], QDC[1], QDC[2], QDC[3], QDC[4], QDC[5], QDC[6], QDC[7], QDC[8], QDC[9], QDC[10], QDC[11], QDC[12], QDC[13], QDC[14], QDC[15], QDC[16], QDC[17], QDC[18], QDC[19], QDC[20], QDC[21], QDC[22], QDC[23], QDC[24], QDC[25], QDC[26], QDC[27], QDC[28], QDC[29], QDC[30], QDC[31]);
					}	
*/					
				}
//[B. End]		
				
					
				} //End loop over entries
				
				//ledfreq stats
				gRunvFreq->SetPoint(filesScanned,run,float(ledcount)/float(duration));	
				
				if (kmcfilebit) kmcfilecount++;
				if (kecfilebit) kecfilecount++;
				if (ksofilebit) ksofilecount++;
				if (kbtsfilebit) kbtsfilecount++;
				if (koosfilebit) koosfilecount++;
				if (kdcfilebit) kdcfilecount++;
				if (khwcmfilebit) khwcmfilecount++;
				if (ledcount == 0) kLEDsOffFilecount++;
				if (ksoerrorfilebit) ksoerrorfilecount++;


				printf(" Run: %i   # of LEDs: %i \n",run,ledcount);
				filesScanned++;
				//return run to zero before reading next line (so directory won't be read as a run)
				run = 0;
			
			
			} //end loop over runs
		
		} //end loop over InputList (exits this loop when end of input list file is reached)
		InputList.close();

	} //end of InputList if statement
	
	//get histogram statistics
	cout << " i | pedestal | mean | sigma | threshold " << endl;
	for (Int_t i=0; i<numPanels; i++){
		Float_t pedestal = hRawQDC[i]->GetMaximumBin();	
		
		Float_t mean = ledqdcsum[i]/float(totalledcountPanel[i]); 				
		Float_t rmssquared = ledqdcsquaredsum[i]/float(totalledcountPanel[i]);	
		Float_t sigma = sqrt(rmssquared - mean*mean);				
		
		threshold[i] = pedestal + 2.0*sigma;
		
		cout << i << " | " << pedestal << " | " << mean << " | " << sigma << " | " << threshold[i] << endl;
		
	}	
	
	//write out graphs and histograms
	TotalMultiplicity->Write("TotalMultiplicity",TObject::kOverwrite);
	
	gRunvFreq->SetTitle("run number vs ledcount/duration");
	gRunvFreq->GetXaxis()->SetTitle("Run Number");
	gRunvFreq->GetYaxis()->SetTitle("Measured LED Frequency");
	gRunvFreq->SetMarkerColor(4);
	gRunvFreq->SetMarkerStyle(21);
	gRunvFreq->SetMarkerSize(0.5);
	gRunvFreq->SetLineColorAlpha(kWhite,0);
	gRunvFreq->Write("gRunvFreq",TObject::kOverwrite);
	
	
	TDirectory *rawqdc = RootFile->mkdir("RawQDC");
	TDirectory *FileDT = RootFile->mkdir("FileDT");
	TDirectory *FileMultvTime = RootFile->mkdir("FileMultvTime");
	
	for (Int_t i = 0; i<numPanels; i++){
		RootFile->cd("RawQDC");
		hRawQDC[i]->Write();
	}	
	
	for (Int_t i = 0; i<filesToScan; i++){
		RootFile->cd("FileDT");
		hDTFile[i].SetTitle("Delta T of sequential entries");	//c:569 check
		hDTFile[i].GetXaxis()->SetTitle("Delta T (seconds)");
		hDTFile[i].GetYaxis()->SetTitle("Multiplicity");
		hDTFile[i].Write();
	}	
	
	for (Int_t i = 0; i<filesToScan; i++){
		RootFile->cd("FileMultvTime");
		gMultvTimeFile[i].SetTitle("low threshold number of panels hit vs event #");
		gMultvTimeFile[i].GetXaxis()->SetTitle("Event Number");
		gMultvTimeFile[i].GetYaxis()->SetTitle("Multiplicity");
		gMultvTimeFile[i].SetMarkerColor(4);
		gMultvTimeFile[i].SetMarkerStyle(21);
		gMultvTimeFile[i].SetMarkerSize(0.5);
		gMultvTimeFile[i].SetLineColorAlpha(kWhite,0);
		gMultvTimeFile[i].Write();
	}
	
	RootFile->cd();
	
//[C. Begin] print out summary
	{
	//print out error summary
	printf ("\n \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ END OF SCAN //////////////////////\n");

	printf("Error Summary: \n");
	printf("Runs scanned: %i\n",filesToScan);
	printf(" 0: Missing Channels (< 32 veto datas in event) : %i of %i runs (%.2f %%) and %i of %i events (%.2f %%) \n"
		,kmcfilecount,filesToScan,100*(double)kmcfilecount/filesToScan,kmcevcount,totalnentries,100*(double)kmcevcount/totalnentries);
	printf(" 1: Extra Channels (> 32 veto datas in event) : %i of %i runs (%.2f %%) and %i of %i events (%.2f %%) \n"
		,kecfilecount,filesToScan,100*(double)kecfilecount/filesToScan,kecevcount,totalnentries,100*(double)kecevcount/totalnentries);
	printf(" 2: Duplicate Channels (any channel shows up multiple times) : %i of %i runs (%.2f %%) and %i of %i events (%.2f %%)\n"
		,kdcfilecount,filesToScan,100*(double)kdcfilecount/filesToScan,kdcevcount,totalnentries,100*(double)kdcevcount/totalnentries);
	printf(" 3: Bad Timestamps: %i of %i runs (%.2f %%) and %i of %i events (%.2f %%)\n"
		,kbtsfilecount,filesToScan,100*(double)kbtsfilecount/filesToScan,kbtsevcount,totalnentries,100*(double)kbtsevcount/totalnentries);
	printf(" 4: Scaler found w/ no QDC data: %i of %i runs (%.2f %%) and %i of %i events (%.2f %%)\n"
		,ksofilecount,filesToScan,100*(double)ksofilecount/filesToScan,ksoevcount,totalnentries,100*(double)ksoevcount/totalnentries);
	printf(" 5: fIndex Errors: (QDCIndex - ScalerIndex != 1 or 2) : %i of %i runs (%.2f %%) and %i of %i events (%.2f %%)\n"
		,koosfilecount,filesToScan,100*(double)koosfilecount/filesToScan,koosevcount,totalnentries,100*(double)koosevcount/totalnentries);
	printf(" 6: EventCount Errors: (SEC - QEC != 1 or 2) : %i of %i runs (%.2f %%) and %i of %i events (%.2f %%)\n"
		,khwcmfilecount,filesToScan,100*(double)khwcmfilecount/filesToScan,khwcmevcount,totalnentries,100*(double)khwcmevcount/totalnentries);
		
	printf(" 7: LED's not activated this run : %i of %i runs (%.2f %%) and %i of %i events (%.2f %%)\n"
		,kLEDsOffFilecount,filesToScan,100*(double)kLEDsOffFilecount/filesToScan,kLEDsOffevcount,totalnentries,100*(double)kLEDsOffevcount/totalnentries);

	printf(" 8: No QDC entries but bit kScalerOnly not thrown : %i of %i runs (%.2f %%) and %i of %i events (%.2f %%)\n"
		,ksoerrorfilecount,filesToScan,100*(double)ksoerrorfilecount/filesToScan,ksoerrorevcount,totalnentries,100*(double)ksoerrorevcount/totalnentries);
/*// not sure how to check these bits
		
			
	
	printf(" 9: Unknown card in crate : %i of %i runs (%.2f %%) and %i of %i events (%.2f %%)\n"
		,errorCount9,filesToScan,100*(double)errorCount9/filesToScan,eventError9,nentries,100*(double)eventError9/nentries);
	printf(" 10: Incorrect run number in MJTRun vRun : %i of %i runs (%.2f %%) and %i of %i events (%.2f %%)\n"
		,keccount0,filesToScan,100*(double)keccount0/filesToScan,eventError10,nentries,100*(double)eventError10/nentries);
*/	

	}
//[C. End]

RootFile->Close();
delete[] hDTFile;
delete[] gMultvTimeFile;
delete gRunvFreq;

printf("Total LED Count = %i\n\n Total # of entries %i\n\n Total Duration %f\n\n",totalledcount, totalnentries, totalduration);
cout << "Wrote ROOT file." << endl;

}	//end of program