#include "vetoScan.hh"

void vetoPerformance(string file)
{	
	// Set up output file
	string filename = file;
	filename.erase(filename.find_last_of("."), string::npos); //remove extension
	filename = filename.substr(filename.find_last_of("\\/")+1,string::npos); // remove path to run

	char OutputFile[200];
	sprintf(OutputFile,"./output/VPC_%s.root",filename.c_str());
	TFile *RootFile = new TFile(OutputFile, "RECREATE");	
  	TH1::AddDirectory(kFALSE); // Global flag: "When a (root) file is closed, all histograms in memory associated with this file are automatically deleted."
	int filesToScan = GetNumFiles(file);

 	bool PlotMultiplicity = true;	// flag to plot multiplicity for EACH RUN
	TH1D *TotalMultiplicity = new TH1D("TotalMultiplicity","Events over threshold",32,0,32);
 	TotalMultiplicity->GetXaxis()->SetTitle("number of panels hit");
	TH1F *hRawQDC[32];  
	TH1F *hCutQDC[32];
	TH1F *hThreshQDC[32];
	TH1F *hLEDCutQDC[32];
 	const int nqdc_bins=4200;
	const int ll_qdc=0;
	const int ul_qdc=4200;
	Char_t hname[50];
	for (Int_t i=0; i<32; i++){
		sprintf(hname,"hRawQDC%d",i);
		hRawQDC[i] = new TH1F(hname,hname,nqdc_bins,ll_qdc,ul_qdc);
		sprintf(hname,"hCutQDC%d",i);
		hCutQDC[i] = new TH1F(hname,hname,nqdc_bins,ll_qdc,ul_qdc);
		sprintf(hname,"hThreshQDC%d",i);
		hThreshQDC[i] = new TH1F(hname,hname,500,ll_qdc,500);
		sprintf(hname,"hLEDCutQDC%d",i);
		hLEDCutQDC[i] = new TH1F(hname,hname,500,ll_qdc,500);
	}

  	// Loop over files in dataset
  	int run = 0;			// run number
  	long start = 0;			// unix timestamps
  	long stop = 0;
	int filesScanned = 0;	// counter
	vector<int> xaxis[7];	// save run numbers for TGraph
	vector<int> yaxis[7];	// save veto bit error counts 

	//int vErrorCountTotal[11] = {0};
	long TotalEntries = 0;
	//long totalDuration = 0;

	// run error statistics
	int errorCount[11] = {0};

	// event error statistics
	int eventError[11] = {0};

	// Keep track of runs with bad errors
	ofstream badFiles;
	badFiles.open("./output/P3KJR_BadRuns.txt");
	
  	ifstream InputList;
  	InputList.open(file.c_str());
	while(!InputList.eof())
	{
		InputList >> run;
		for(int i = 0; i < 7; i++) xaxis[i].push_back(run);
		filesScanned++;
		GATDataSet ds(run);
		start = GetStartUnixTime(ds);
		stop = GetStopUnixTime(ds);
			
		// standard veto initialization block
		// with some error checking
		TChain *v = ds.GetVetoChain();
		long vEntries = v->GetEntries();
		MJTRun *vRun = new MJTRun();
		MGTBasicEvent *vEvent = new MGTBasicEvent(); 
		MJTVetoData *vData[32];
		unsigned int mVeto = 0;
		uint32_t vBits = 0;
		if (v->GetListOfBranches()->FindObject("run")) 
			v->SetBranchAddress("run",&vRun);
		else { cout << "Couldn't find branch: \"run\"" << endl; break; }
		if (v->GetListOfBranches()->FindObject("mVeto")) 
			v->SetBranchAddress("mVeto",&mVeto);
		else { cout << "Couldn't find branch: \"mVeto\"" << endl; break; }
		if (v->GetListOfBranches()->FindObject("vetoEvent")) 
			v->SetBranchAddress("vetoEvent",&vEvent);
		else { cout << "Couldn't find branch: \"vetoEvent\"" << endl; break; }
		if (v->GetListOfBranches()->FindObject("vetoBits")) 
			v->SetBranchAddress("vetoBits",&vBits);
		else { cout << "Couldn't find branch: \"vetoBits\"" << endl; break; }
		v->GetEntry(0);
		TotalEntries += vEntries;

		printf("\n ========= Scanning Run %i: %li entries. =========\n",run,vEntries);

		// only write these if their bools are set = true
		//TH1D *CorruptionInTime = new TH1D("CorruptionInTime","corrupted entries during run",(int)duration/5,0,(int)duration);
		//CorruptionInTime->GetXaxis()->SetTitle("time (5 sec / bin)");
		TH1D *OneRunMultiplicity = new TH1D("multiplicity","multiplicity of veto entries",32,0,32);
		OneRunMultiplicity->GetXaxis()->SetTitle("number of panels hit");		
		
		// Run-level checks
		
		// veto error variables
		bool vBitTripped = false;
		int eventsSkipped = 0;
		//int vPrevErrorCount[11] = {0};
		int vErrorCount[11] = {0};
		bool LEDsOff = true;

		bool error[11] = {false};

		MJVetoEvent prev;

		// Custom SW Threshold (obtained from vetoTimeFinder)
		int thresh[32] = {213,241,311,223,223,279,293,189,301,331,185,231,291,219,273,169,55,141,83,115,89,93,109,73,99,79,69,101,89,109,115,105};

		// Loop over entries
		//for (long i = 0; i < 10; i++) 
		//for (long i = vEntries-10; i<vEntries; i++)
		for (int i = 0; i < vEntries; i++)
		{
			v->GetEntry(i);

			// Fill MJVetoEvent object
			MJVetoEvent veto;
			veto.SetQDCThreshold(thresh);	// kept separate to allow the user to specify their own.
    		int isGood = veto.WriteEvent(vRun,vEvent,vData,vBits,run);
			
			// Check veto bits for errors (defined in MGDO/Majorana/MJTypes)
			if (MJBits::GetBit(vBits, MJVetoBits::kMissingChannels)) 	{vErrorCount[0]++; 	error[0]=true;}
			if (MJBits::GetBit(vBits, MJVetoBits::kExtraChannels)) 		{vErrorCount[1]++; 	error[1]=true;}
			if (MJBits::GetBit(vBits, MJVetoBits::kDuplicateChannel)) 	{vErrorCount[2]++; 	error[2]=true;}
			if (MJBits::GetBit(vBits, MJVetoBits::kBadTimeStamp)) 		{vErrorCount[3]++; 	error[3]=true;}
			if (MJBits::GetBit(vBits, MJVetoBits::kQDCOutOfSequence)) 	{vErrorCount[4]++; 	error[4]=true;}
			if (MJBits::GetBit(vBits, MJVetoBits::kScalerOnly)) 		{vErrorCount[5]++; 	error[5]=true;}
			if (MJBits::GetBit(vBits, MJVetoBits::kHWCountMismatch)) 	{vErrorCount[6]++; 	error[6]=true;}
			if (isGood == -1) {vErrorCount[8]++;   	error[8]=true;}	// QDC entries don't exist but kScalerOnly was not tripped
    		if (isGood == -2) {vErrorCount[9]++;   	error[9]=true;} // Unknown card in crate
    		if (isGood == -3) {vErrorCount[10]++;  	error[10]=true;}// vEvent run number doesn't match the current run

			if (isGood != 1) {
				cout << "isGood:" << isGood << endl;
				veto.Print();
			}
			/*
			printf("Count: %i  Run: %i  Entry: %i  Time: %.5f  %.2f%%.\n"
				,veto.run,i,veto.timeSec,((double)i/vEntries)*100); 
					
			// Info for fIndex and EventCount mismatches.
			printf("   Channels: %lu  scaler_Ind: %li  qdc1_Ind: %li  qdc2_Ind: %li  s-QDC1: %li  s-qdc2: %li  SEC: %li  QEC: %li\n\n"
				,veto.QDCchans,veto.ScalerIndex,veto.QDC1Index,veto.QDC2Index
				,veto.ScalerIndex - veto.QDC1Index, veto.ScalerIndex - veto.QDC2Index
				,veto.SEC,veto.QEC);
			*/


    		/*
			for (int j = 0; j < 10; j++) 
			{ 
				// Show first instance of a bit being thrown.
				//if (vErrorCount[j] == 1 && vPrevErrorCount[j] == 0) 
				
				// Show all instances of bits being thrown.
				if (vErrorCount[j] > 0) 

				// Show the first and second instance of bits being thrown.
				//if ((vErrorCount[j] == 1 && vPrevErrorCount[j] == 0) || (vErrorCount[j] == 2 && vPrevErrorCount[j] == 1))
				{
					vBitTripped = true;
					printf(" New Error: Bit %i  Count: %i  Run: %i  Entry: %i  Time: %.5f  %.2f%%.\n"
						,j,vErrorCount[j],veto.run,i,veto.timeSec,((double)i/vEntries)*100); 
					
					// Info for fIndex and EventCount mismatches.
					printf("   Channels: %lu  scaler_Ind: %li  qdc1_Ind: %li  qdc2_Ind: %li  s-QDC1: %li  s-qdc2: %li  SEC: %li  QEC: %li\n\n"
						,veto.QDCchans,veto.ScalerIndex,veto.QDC1Index,veto.QDC2Index
						,veto.ScalerIndex - veto.QDC1Index, veto.ScalerIndex - veto.QDC2Index
						,veto.SEC,veto.QEC);

					// Show previous event.
					if (i>0) {
						printf("   Last event: Run: %i  Entry: %i  Time: %.5f  %.2f%%.\n"
							,prev.run,i-1,prev.timeSec,((double)(i-1)/vEntries)*100);
						
						printf("   Channels: %lu  scaler_Ind: %li  qdc1_Ind: %li  qdc2_Ind: %li  s-QDC1: %li  s-qdc2: %li  SEC: %li  QEC: %li\n\n"
							,prev.QDCchans,prev.ScalerIndex,prev.QDC1Index,prev.QDC2Index
							,prev.ScalerIndex - prev.QDC1Index, prev.ScalerIndex - prev.QDC2Index
							,prev.SEC,prev.QEC);
					}
				}
			}
			memcpy(vPrevErrorCount, vErrorCount, sizeof(vPrevErrorCount));
			*/

			if (i > 0) prev = veto;
			else prev.Clear();

			if (isGood == 0 || isGood == -1 || isGood == -2 || isGood == -3) {
				eventsSkipped++;
				veto.Clear();
				prev = veto;
				continue;
			}
 
			// Check for LED entries in a simple way
			if (veto.multip > 25) LEDsOff = false;

			// Compare multiplicities
			//printf("mVeto: %i  multip: %i\n",mVeto,veto.multip);

			// Multiplicity and QDC plots
			TotalMultiplicity->Fill(veto.multip);			
	        if (PlotMultiplicity) OneRunMultiplicity->Fill(veto.multip);
			for (int k = 0; k < 32; k++) {
				hRawQDC[k]->Fill(veto.QDC[k]);
				if (veto.QDC[k] < 500) hThreshQDC[k]->Fill(veto.QDC[k]);
				if (veto.QDC[k] > veto.SWThresh[k]) hCutQDC[k]->Fill(veto.QDC[k]);
			}
	
		} // End loop over entries.
		
		printf(" ================= End of Run %i. =================\n",run);

		printf(" Duration: %li seconds.\n",stop-start);

		if (error[0]) { errorCount[0]++; badFiles << run << endl; }
		if (!LEDsOff) { vErrorCount[7] = 1; errorCount[7]++;}
		for (int r = 1; r < 11; r++) {
			if (error[r] && r != 7) errorCount[r]++;
		}
		
		if (PlotMultiplicity) {
			char outfile2[200];
			sprintf(outfile2,"Multiplicity_Run%i",run);
			OneRunMultiplicity->Write(outfile2,TObject::kOverwrite);
		}

		if (vBitTripped) {
			double scalerCorruption = (((double)vErrorCount[3])/vEntries)*100;
			printf(" Skipped %i events. (%.2f%% of total). \n Veto Errors: \n",eventsSkipped,((double)eventsSkipped/vEntries)*100);
			if (vErrorCount[0] > 0) printf(" 0: Missing Channels (< 32 veto datas in event) : %i\n",vErrorCount[0]);
			if (vErrorCount[1] > 0) printf(" 1: Extra Channels (> 32 veto datas in event) : %i\n",vErrorCount[1]);
			if (vErrorCount[2] > 0) printf(" 2: Duplicate Channels (any channel shows up multiple times) : %i\n",vErrorCount[2]);
			if (vErrorCount[3] > 0) printf(" 3: Bad Timestamps: %i of %li entries: %.2f%%.\n",vErrorCount[3],vEntries,scalerCorruption);
			if (vErrorCount[4] > 0) printf(" 4: Scaler found w/ no QDC data: %i\n",vErrorCount[4]);
			if (vErrorCount[5] > 0) printf(" 5: fIndex Errors: (QDCIndex - ScalerIndex != 1 or 2) : %i\n",vErrorCount[5]);
			if (vErrorCount[6] > 0) printf(" 6: EventCount Errors: (SEC - QEC != 1 or 2) : %i\n",vErrorCount[6]);
			if (vErrorCount[7] > 0) printf(" 7: LED's not activated this run : %i\n",vErrorCount[7]);
			if (vErrorCount[8] > 0) printf(" 8: No QDC entries but bit 4 not thrown : %i\n",vErrorCount[8]);
			if (vErrorCount[9] > 0) printf(" 9: Unknown card in crate : %i \n",vErrorCount[9]);
			if (vErrorCount[10]> 0) printf(" 10: Incorrect run number in MJTRun vRun : %i \n",vErrorCount[10]);
			cout << " ====================================================\n";
		}

		for (int x = 0; x < 11; x++) {
			eventError[x] += vErrorCount[x];
		}

		// Push veto bits into yaxis for TGraph
		for (int q = 0; q < 7; q++) {
			yaxis[q].push_back(vErrorCount[q]);
		}

	} // End loop over file list

	printf ("\n \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ END OF SCAN //////////////////////\n");

	printf("Error Summary: \n");
	printf("Runs scanned: %i\n",filesToScan);
	printf(" 0: Missing Channels (< 32 veto datas in event) : %i of %i runs (%.2f %%) and %i of %li events (%.2f %%) \n"
		,errorCount[0],filesToScan,100*(double)errorCount[0]/filesToScan,eventError[0],TotalEntries,100*(double)eventError[0]/TotalEntries);
	
	printf(" 1: Extra Channels (> 32 veto datas in event) : %i of %i runs (%.2f %%) and %i of %li events (%.2f %%) \n"
		,errorCount[1],filesToScan,100*(double)errorCount[1]/filesToScan,eventError[1],TotalEntries,100*(double)eventError[1]/TotalEntries);
	
	printf(" 2: Duplicate Channels (any channel shows up multiple times) : %i of %i runs (%.2f %%) and %i of %li events (%.2f %%)\n"
		,errorCount[2],filesToScan,100*(double)errorCount[2]/filesToScan,eventError[2],TotalEntries,100*(double)eventError[2]/TotalEntries);
	
	printf(" 3: Bad Timestamps: %i of %i runs (%.2f %%) and %i of %li events (%.2f %%)\n"
		,errorCount[3],filesToScan,100*(double)errorCount[3]/filesToScan,eventError[3],TotalEntries,100*(double)eventError[3]/TotalEntries);
	
	printf(" 4: Scaler found w/ no QDC data: %i of %i runs (%.2f %%) and %i of %li events (%.2f %%)\n"
		,errorCount[4],filesToScan,100*(double)errorCount[4]/filesToScan,eventError[4],TotalEntries,100*(double)eventError[4]/TotalEntries);
	
	printf(" 5: fIndex Errors: (QDCIndex - ScalerIndex != 1 or 2) : %i of %i runs (%.2f %%) and %i of %li events (%.2f %%)\n"
		,errorCount[5],filesToScan,100*(double)errorCount[5]/filesToScan,eventError[5],TotalEntries,100*(double)eventError[5]/TotalEntries);
	
	printf(" 6: EventCount Errors: (SEC - QEC != 1 or 2) : %i of %i runs (%.2f %%) and %i of %li events (%.2f %%)\n"
		,errorCount[6],filesToScan,100*(double)errorCount[6]/filesToScan,eventError[6],TotalEntries,100*(double)eventError[6]/TotalEntries);
	
	printf(" 7: LED's not activated this run : %i of %i runs (%.2f %%) and %i of %li events (%.2f %%)\n"
		,errorCount[7],filesToScan,100*(double)errorCount[7]/filesToScan,eventError[7],TotalEntries,100*(double)eventError[7]/TotalEntries);
	
	printf(" 8: No QDC entries but bit 4 not thrown : %i of %i runs (%.2f %%) and %i of %li events (%.2f %%)\n"
		,errorCount[8],filesToScan,100*(double)errorCount[8]/filesToScan,eventError[8],TotalEntries,100*(double)eventError[8]/TotalEntries);
	
	printf(" 9: Unknown card in crate : %i of %i runs (%.2f %%) and %i of %li events (%.2f %%)\n"
		,errorCount[9],filesToScan,100*(double)errorCount[9]/filesToScan,eventError[9],TotalEntries,100*(double)eventError[9]/TotalEntries);
	
	printf(" 10: Incorrect run number in MJTRun vRun : %i of %i runs (%.2f %%) and %i of %li events (%.2f %%)\n"
		,errorCount[10],filesToScan,100*(double)errorCount[10]/filesToScan,eventError[10],TotalEntries,100*(double)eventError[10]/TotalEntries);


	badFiles.close();


	/*
	//TotalCorruptionInTime->Write("TotalCorruptionInTime",TObject::kOverwrite);
	TotalMultiplicity->Write("TotalMultiplicity",TObject::kOverwrite);

	// Plot errors vs. run number
	// Find graph parameters
	int lastxaxis = 0;
	for (int i = 0; i < xaxis[0].size(); i ++ ) 
		lastxaxis = xaxis[0][i];
	int lastyaxis = 0;
	for (int i = 0; i < 7; i++) {
		for (int j = 0; j < yaxis[i].size(); j++) {
			if (yaxis[i][j] > lastyaxis) lastyaxis = yaxis[i][j];
		}
	}
	TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas");
	c1->SetGrid();
	TGraph *g[7];
	TH1F *hr = c1->DrawFrame(xaxis[0][0],0,lastxaxis,lastyaxis);
	hr->SetXTitle("Run");
	hr->SetYTitle("Errors");
	hr->SetTitle("Veto Errors");
	TLegend *legend = new TLegend(0.8,0.58,0.9,0.9);
	for (int i = 0; i < 7; i++){ 
		g[i] = new TGraph(filesScanned,&(xaxis[i][0]),&(yaxis[i][0]));
		g[i]->SetLineColor(color(i));
		g[i]->SetMarkerColor(color(i));
		g[i]->SetMarkerStyle(20);
		g[i]->SetLineWidth(2);	// want to make error bars thicker
		g[i]->Draw("LP");
	}
	legend->AddEntry(g[0],"Missing Channels","lep");
	legend->AddEntry(g[1],"Extra Channels","lep");
	legend->AddEntry(g[2],"Duplicate Channels","lep");
	legend->AddEntry(g[3],"Bad Timestamps","lep");
	legend->AddEntry(g[4],"Scaler No QDC","lep");
	legend->AddEntry(g[5],"fIndex Mismatch","lep");
	legend->AddEntry(g[6],"EventCount Mismatch","lep");
	legend->Draw();
	c1->Write("ErrorPlot",TObject::kOverwrite);
	delete c1;

	// QDC plots & calibration table:
	// gaus: A gaussian with 3 parameters: f(x) = p0*exp(-0.5*((x-p1)/p2)^2)).
	TF1 *fits[32];
	TCanvas *vcan0 = new TCanvas("FittedCutQDC","cut & fitted QDC, panels 1-32",0,0,800,600);
	vcan0->Divide(8,4,0,0);
	TCanvas *vcan1 = new TCanvas("QDCThresholds","veto QDC thresholds, panels 1-32",0,0,800,600);
	vcan1->Divide(8,4,0,0);	
	Char_t buffer[2000];
	printf("\n  LED Peak Information:\n  Panel / Mean,error / Sigma,error / Chi-square/NDF (~1?)\n");
	for (Int_t i=0; i<32; i++)
	{	
		vcan0->cd(i+1);
		TVirtualPad *vpad0 = vcan0->cd(i+1); vpad0->SetLogy();
		hThreshQDC[i]->Write(); // write the low-QDC part of the spectrum separately
		hCutQDC[i]->Write();    // write the cut QDC

		// if we have entries above the cut threshold, fit them
		if (hCutQDC[i]->GetEntries() > 0) {
			hRawQDC[i]->Write();
			hCutQDC[i]->Fit("gaus","q");
			hCutQDC[i]->Draw();
			fits[i] = hCutQDC[i]->GetFunction("gaus");
			float NDF = (float)fits[i]->GetNDF();
			if (fits[i]->GetNDF() == 0) NDF = 0.0000001;
			sprintf(buffer,"%i  %.1f  %.1f  %.1f  %.1f  %.1f",i,
				fits[i]->GetParameter(1),fits[i]->GetParError(1),		// mean (LED center)
				fits[i]->GetParameter(2),fits[i]->GetParError(2),		// sigma
				fits[i]->GetChisquare()/NDF); // X^2/NDF (closer to 1, better fit.)
			cout << "  " << buffer << endl;
		}
		else {
			hRawQDC[i]->Write();	// just write the raw QDC without fitting
		}	    	
   		// plot low-QDC range separately
  		vcan1->cd(i+1);
  		TVirtualPad *vpad1 = vcan1->cd(i+1); vpad1->SetLogy();
  		hThreshQDC[i]->Draw();
	}
	vcan0->Write((filename+"_QDC").c_str());
	vcan1->Write((filename+"_ThreshQDC").c_str());
	*/

	// All done!
	RootFile->Close();

}