// For analysis of the new MJD built data format.
// Goal is to give an example of accessing every data member.
// Run in compiled mode with .x vetoReader.C++
// 
// Clint Wiseman, University of South Carolina
// 11/16/2015

/*
	File components:
	1. headerXML
	2. ProcessID
	3. builderID0
	4. MGTree
	5. MGGarbageTree
	6. ChannelMap
	7. ChannelSettings
	8. VetoTree
*/

#ifndef __CINT__
//#include <vector>
#include <iostream>
//#include <fstream>
//#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
//#include "TChain.h"
//#include "TEntryList.h"
//#include "TBranch.h"
//#include "TH1.h"
//#include "CLHEP/Units/SystemOfUnits.h"
//#include "GATMultiplicityProcessor.hh"
#include "GATDataSet.hh"
//#include "MGTWaveform.hh"
#include "MJTRun.hh"
#include "MJTVetoData.hh"
#include "MGTBasicEvent.hh"
using namespace std;
//using namespace CLHEP;
#endif

void vetoReader(){

	// Initialize with standard ROOT methods	
	/*
	TFile *f = new TFile("~/dev/datasets/newModuleOne/built/OR_run6947.root");
	TTree *v = (TTree*)f->Get("VetoTree");
	TTree *b = (TTree*)f->Get("MGTree");
	Long64_t nentries = v->GetEntries();
	cout << "Found " << nentries << " entries." << endl;
	//v->Show(nentries-1);	// for comparison
	//b->Show(0);	// for comparison
	*/

	// Initialize by GATDataSet
	GATDataSet ds(6960);
	TChain *b = ds.GetBuiltChain();
	TChain *v = ds.GetVetoChain();	// Clint's super-fancy addition
	Long64_t nentries = v->GetEntries();
	cout << "Found " << nentries << " veto entries." << endl;
	//v->Show(nentries-1);	// for comparison
	//b->Show(0); // for comparison

	// Must initialize with GATDataSet to use these functions
	// Does not seem to contain veto data.
	// It would be nice to see the HV map information ...
	MJTChannelMap *map = ds.GetChannelMap();
	//map->DumpChannelMap();
	MJTChannelSettings *set = ds.GetChannelSettings();
	//set->DumpSettings();
	//set->DumpEnabledIDs();
	vector<uint32_t> en = set->GetEnabledIDList(); 	// save a vector with the enabled Ge detector list

	cout << "=================================================================" << endl;

	MJTRun *VetoRun = new MJTRun();
	MGTBasicEvent *vetoEvent = new MGTBasicEvent(); 
	UInt_t mVeto = 0;
	
	v->SetBranchAddress("run",&VetoRun);
	v->SetBranchAddress("mVeto",&mVeto);
	v->SetBranchAddress("vetoEvent",&vetoEvent);

	v->GetEntry(0);
	//v->GetEntry(nentries-1);

	// =========================================================================
	// 1. "run"
	// Both MGTree and VetoTree now contain a MJTRun object. (MJTRun.hh)
	
	// This is "packed", returns e.g. 2576980377.
	// Shouldn't need to unpack it, since we can do the test below.
	uint32_t runbits = VetoRun->GetRunBits(); 
	
	// Prints a list of enabled run bits with names to cout.
	cout << "Enabled run bits: " << endl;
	VetoRun->ListRunBits();	

	// Other bits are listed in $MGDODIR/Majorana/MJTypes.hh 
	// Returns a bool (0 or 1)
	cout << "Test for a run bit - Shop Activity: " << VetoRun->GetRunBit(MJRunBits::kMachineShop) << endl; 

	// Unsigned int from MGTypes.hh -- kData=0, kTest=1, kCalibration=2, kMC=3, kUndefined=4																	  
	printf("Run Type: %u\n",VetoRun->GetRunType());  
	printf("Run Number: %i\n",VetoRun->GetRunNumber());

	// Set in /MJOR/MOVetoDataLoader.cc:    fRun->SetRunDescription("Orca run");
	// Returns a string.
	cout << "Run Description: " << VetoRun->GetRunDescription() << endl;  
	
	// Returns unix timestamps (GMT)
	printf("Start time:%li  Stop time: %li\n",VetoRun->GetStartTime(),VetoRun->GetStopTime()); 

	// Both of these return strings.
	// "undefined" is the default in MGRun.cc
	cout << "Parent DAQ label: " << VetoRun->GetParentDAQLabel() << endl;
	cout << "MGDO Conversion Version: " << VetoRun->GetMGDOConversionVersion() << endl;

	// =========================================================================

	// 2. "mVeto"
	// Built-in multiplicity is calculated in MOVetoDataLoader.cc:
	// if(!((MJTVetoData*)(fVetoEvent->GetDetectorData(i)))->IsUnderThreshold())
	// So it's dependent on the IsUnderThreshold software tag set in ORCA.
	//   which is QDC = 500.
	printf("mVeto: %u\n",mVeto);

	// =========================================================================

	// 3. "vetoEvent"
	// Veto data is stored in an MGTBasicEvent filled with MJTVetoData objects.
	
	// From adding all QDC amplitudes in MOVetoDataLoader.cc
	// It looks like it would add the QDC pedestal values as well, which are unphysical
	printf("Total energy? (QDC+pedestal): %.1f \n",vetoEvent->GetETotal());

	// Returns a double: time since the run start (in ns?)
	printf("Time (ns?): %.1f \n",vetoEvent->GetTime());

	// Set in MGTypes.hh -- kReal=0, kPulser=1, kMC=2, kUndefined=3
	printf("Veto event type: %u\n",vetoEvent->GetEventType());

	// Make sure NDetectorData == 32.
	// Returns a size_t.
	printf("NDetectorData: %lu\n",vetoEvent->GetNDetectorData());
	if (vetoEvent->GetNDetectorData() != 32) 
		cout << "Warning! Detector Data is not 32, it's " << vetoEvent->GetNDetectorData() << endl;
	
	// Access the MJTVetoData objects "vd"
	MJTVetoData *vd[32];
	for (int i=0; i<32; i++) { vd[i] = dynamic_cast<MJTVetoData*>(vetoEvent->GetDetectorData()->At(i)); }

	// These values should be the same for every MJTVetoData object in the array
	// Check: If EventCount doesn't match ScalerCount, something has gone wrong.
	printf("Crate: %i  EventCount:%i  ScalerCount: %i\n",
		vd[0]->GetCrate(),vd[0]->GetEventCount(),vd[0]->GetScalerCount());
	printf("ScalerID: %i  ScalerIndex:%llu  TimeStamp: %llu  IsBadTs: %i \n",
		vd[0]->GetScalerID(),vd[0]->GetScalerIndex(),vd[0]->GetTimeStamp(),vd[0]->IsBadTS());
	
	if (vd[0]->GetEventCount() != vd[0]->GetScalerCount()) 
		printf("Warning!  EventCount and ScalerCount don't match!\n");
	
	// The scaler clock is 100MHz = 1E8 counts / sec
	printf("Scaler Time (sec): %.8f\n",vd[0]->GetTimeStamp()/1E8);

	// Display raw data (out of order)
	//   Member functions from MJTVetoData.hh and MGDetectorData.hh (MJTVetoData's base class)
	//   fIndex (ORCA packet number, useful for consistency checks)
	//   fID ("packed" channel number, just like Ge detectors)
	printf("Unsorted MJTVetoData:\n");
	for (int j = 0; j<32; j++)	{
		printf("Ca: %-3i  Ch: %-3i  Amp: %-5.0f  UTh: %-3i  OF: %-3i  ID: %-5i  Index: %-5llu\n",
			vd[j]->GetCard(),vd[j]->GetChannel(),vd[j]->GetAmplitude(),vd[j]->IsUnderThreshold(),
			vd[j]->IsOverflow(),vd[j]->GetID(),vd[j]->GetIndex());
	}
	cout << endl;

	// Sort raw data into arrays and then display.
	// This may not be totally necessary, but makes hit pattern analysis easier
	//   to match to the physical veto panel locations, and Yuri's wiring diagrams.
	// Most things are cast to int's.  
	// Original types can be found in MJTVetoData.hh and MGDetectorData.hh if necessary.
	const int card1 = 13;
	const int card2 = 18;
	int Card[32] = {0};
	int QDC[32] = {0};
	int IsUnderThreshold[32] = {0};
	int IsOverflow[32] = {0};
	int ID[32] = {0};
	long Index[32] = {0};
	int k = 0;
	for (int j = 0; j<32; j++)	{
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
	printf("Sorted MJTVetoData:\n");
	for (int j = 0; j < 32; j++) {
		printf("Ca: %-3i  Panel: %-3i  Amp: %-5i  UTh: %-3i  OF: %-3i  ID:%-5i  Index:%-5li\n",
			Card[j],j,QDC[j],IsUnderThreshold[j],IsOverflow[j],ID[j],Index[j]);
	}

	// Display ORCA packet numbers: Scaler, QDC Card 1, QDC Card 2
	printf("ORCA Packet Numbers - Scaler: %llu  QDC 1 (Ch.%i): %li  QDC 2 (Ch.%i): %li\n",
		vd[0]->GetScalerIndex(),Card[0],Index[0],Card[31],Index[31]);
	
	// Check: QDC index should be no more than 2 greater than scaler index.
	// The event builder has supposedly already done this check.
	if ((Long_t) vd[0]->GetScalerIndex() > Index[0] || (Long_t)vd[0]->GetScalerIndex() > Index[31])
		printf("Warning! Scaler index is larger than QDC index!\n");
	

	// =========================================================================
	

}