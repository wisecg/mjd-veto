#include "vetoScan.hh"
using namespace std;
/*
	Methods of calculating time of veto events:
	1. Scaler
	2. SBC Time (implemented after run 8557, not retroactive)
	3. LED Time (only accurate to last LED)
	4. Entry Time (not very desirable - consistently under or over predicts the time)
	5. Gretina Interpolation time (haven't tested this yet ... probably relies on ORCA packet indexes)
*/

void vetoTimeFinder(string file) 
{
	// Input a list of run numbers
	ifstream InputList(file.c_str());
	if(!InputList.good()) {
    	cout << "Couldn't open " << file << endl;
    	return;
    }
	
	int run = 0;
	while(!InputList.eof())
	{
		InputList >> run;

		// standard initialization
		GATDataSet *ds = new GATDataSet(run);
		TChain *v = ds->GetVetoChain();
		long vEntries = v->GetEntries();
		MJTRun *vRun = new MJTRun();
		MGTBasicEvent *vEvent = new MGTBasicEvent(); 
		unsigned int mVeto = 0;
		uint32_t vBits = 0;
		v->SetBranchAddress("run",&vRun);
		v->SetBranchAddress("mVeto",&mVeto);
		v->SetBranchAddress("vetoEvent",&vEvent);
		v->SetBranchAddress("vetoBits",&vBits);

		printf("\n=========== Scanning Run %i: %li entries. ===========\n",run,vEntries);

		// time stuff.
		v->GetEntry(0);
		long start = (long)vRun->GetStartTime();
		long stop = (long)vRun->GetStopTime();
		double duration = (double)(stop - start);

		// ===================== FIRST LOOP OVER ENTRIES =========================
		// Find the first non-bad scaler entry and offset the SBC time
		// so that the two times are equal.
		//
		// Also calculate the LED frequency so that we can do LEDTime 
		// in the next loop.
		//
		MJVetoEvent first;
		bool foundFirst = false;
		int firstGoodEntry = 0;
		MJVetoEvent prev;
		char hname[200];
		sprintf(hname,"LEDDeltaT_run%i",run);
		TH1F *LEDDeltaT = new TH1F(hname,hname,100000,0,100); // 0.001 sec/bin
		int isGood = 0;
		int highestMultip = 0;	// try to predict how many panels there are for this run.
		int pureLEDcount = 0;

		for (int i = 0; i < vEntries; i++)
		{
			v->GetEntry(i);
			MJVetoEvent veto;
			veto.SetSWThresh();	
	    	isGood = veto.WriteEvent(i,vRun,vEvent,vBits,run);

			// Save the first good entry number for the SBC offset
			if (isGood == 1 && !foundFirst) {
				first = veto;
				foundFirst = true;
				firstGoodEntry = i;
			}
			// Super simple error check
	    	if (isGood != 1) {
	    		continue;
	    	}
	    	// Super simple LED tag
	    	if (veto.GetMultip() > highestMultip && veto.GetMultip() < 33) highestMultip = veto.GetMultip();
			if (veto.GetMultip() >= 20) { 				
				LEDDeltaT->Fill(veto.GetTimeSec()-prev.GetTimeSec());
				pureLEDcount++;
			}
			prev = veto;
		}

		// Find the SBC offset		
		double SBCOffset = first.GetTimeSBC() - first.GetTimeSec();
		printf("First good entry: %i  SBCOffset: %.2f\n",firstGoodEntry,SBCOffset);

		// Find the LED frequency
		printf("\"Pure\" LED count: %i.  Approx rate: %.3f\n",pureLEDcount,pureLEDcount/duration);
		double LEDrms = 0;
		double LEDfreq = 0;
		int dtEntries = LEDDeltaT->GetEntries();
		if (dtEntries > 0) {
			int maxbin = LEDDeltaT->GetMaximumBin();
			LEDDeltaT->GetXaxis()->SetRange(maxbin-100,maxbin+100); // looks at +/- 0.1 seconds of max bin.
			LEDrms = LEDDeltaT->GetRMS();
			LEDfreq = 1/LEDDeltaT->GetMean();
		}
		else {
			printf("Warning! No multiplicity > 20 events!!\n");
			LEDrms = 9999;
			LEDfreq = 9999;
		}
		double LEDperiod = 1/LEDfreq;
		printf("LED_f: %.8f LED_t: %.8f RMS: %8f\n",LEDfreq,LEDperiod,LEDrms);
		delete LEDDeltaT;


		// ===================== SECOND LOOP OVER ENTRIES =========================
		// 
		vector<double> diffs1;
		vector<double> diffs2;
		vector<double> diffs3;
		vector<double> diffs4;
		vector<double> times;

		MJVetoEvent prevLED;
		// int LEDCount = 0;
		// double ledTime = 0;
		// double firstLEDTime = 0;
		bool prevBadScaler = false;
		// double prev_eTime = 0;
		double SBCError = 0;
		double prevSBCError = 0;
		double pureLEDtime = 0;
		int pureLEDcount2 = 0;

		for (int i = 0; i < vEntries; i++)
		// for (int i = 0; i < 10; i++)
		{
			v->GetEntry(i);
			MJVetoEvent veto;
			veto.SetSWThresh();	
	    	isGood = veto.WriteEvent(i,vRun,vEvent,vBits,run);

	    	if (isGood != 1) 
	    	{
				int error[18] = {0};
	    		veto.UnpackErrorCode(isGood,error);
	    		bool worseProblems = false;
	    		for (int q=0; q<18; q++) 
	    		{
	    			// don't skip bad-scaler events
	    			if (q!=4 && error[q] == 1) worseProblems = true;
	    		}
	    		if (worseProblems) {
	    			cout << "Skipped Entry: " << i  << endl;
	    			// if (i < 500) veto.Print();
	    			continue;
	    		}
	    	}

			// ---------------------------------------------
			// 0. ENTRY TIME
			// This is a "last resort" method.
			// It underestimates before halfway through the run, then starts
			// to overestimate.
			//
			double eTime = ((double)i / vEntries) * duration;
			
			if (!veto.GetBadScaler()) diffs3.push_back(veto.GetTimeSec() - eTime);

	    	// ---------------------------------------------
	    	// 1. SBC TIME
			//
			veto.SetTimeSBC(veto.GetTimeSBC() - SBCOffset);
			SBCError = veto.GetTimeSec() - veto.GetTimeSBC();	// keep track of offset

			if (!veto.GetBadScaler() && !prevBadScaler)
			{
				diffs1.push_back(SBCError);
				diffs2.push_back(SBCError-prevSBCError);
				times.push_back(veto.GetTimeSec());
			}
			prevBadScaler = veto.GetBadScaler();
			prevSBCError = SBCError;

	    	// ---------------------------------------------
	    	// 2. "PURE" LED TIME (don't use scaler information)
	    	// Use the number of LED's from the first scan using 
	    	// only multiplicity information.
			//
			if (veto.GetMultip() > highestMultip-5) {
				pureLEDcount2++;
			}
			pureLEDtime = ((double)pureLEDcount2/pureLEDcount) * duration;

			printf("i:%i  scaler: %.2f  ledTime: %.2f   M: %i  count:%i \t diff:%.2f\n",
				i,veto.GetTimeSec(),pureLEDtime,veto.GetMultip(),pureLEDcount2,veto.GetTimeSec()-pureLEDtime);
			
			if (!veto.GetBadScaler() && !prevBadScaler)
			{
				diffs4.push_back(veto.GetTimeSec()-pureLEDtime);
			}
			
			// ---------------------------------------------
			// 2. LED + SCALER METHOD
			// Tag LED's based mainly on the known LED frequency, which requires the scaler.
			// Don't resort to the SBC time.
			//
			/*
			bool TimeCut = false;
			bool IsLED = false;
			bool force = false;

			// Need to know the time since the last LED.
			bool ApproxTime = false;
			double xTime = 0;
			double x_deltaT = 0;
			double xError = 0;

			// 3 MJVetoEvents: veto, prev, and prevLED.
			// Should implement a count for how often these cases happen.
			if (!veto.GetBadScaler() && !prevLED.GetBadScaler()) {
				xTime = veto.GetTimeSec();
				xError = 0;
				x_deltaT = veto.GetTimeSec()-prevLED.timeSec;
			}
			else if (veto.GetBadScaler() && !prev.GetBadScaler() && !prevLED.GetBadScaler()) {
				xTime = prev.timeSec;
				x_deltaT = prev.timeSec-prevLED.timeSec;
				xError = LEDPeriod;
				ApproxTime = true;
			}
			else if (veto.GetBadScaler() && prev.GetBadScaler() && !prevLED.GetBadScaler()) {
				xTime = prevLED.timeSec;
				x_deltaT = eTime - prevLED.timeSec;
				xError = LEDPeriod;
				ApproxTime = true;
			}
			else if (veto.GetBadScaler() && prev.GetBadScaler() && prevLED.GetBadScaler()) {
				xTime = eTime;
				x_deltaT = eTime - prev_eTime;
				xError = LEDPeriod;
				ApproxTime = true;
			}



			// Get the first LED by multiplicity
			if (veto.GetMultip() > 20 && LEDCount == 0 && !veto.GetBadScaler()) {
				firstLEDTime = veto.GetTimeSec();
				cout << "First LED" << endl;
				LEDCount++;
			}
			// Find LED time 
			if (LEDCount == 1) ledTime = firstLEDTime;
			if (LEDCount > 1) ledTime = firstLEDTime + LEDperiod * (LEDCount-1);

			printf(" M:%-3i   LEDs:%-5i   LEDTime:%-10.4f   timeSec:%-10.4f  diff:%-5.4f\n"
				,veto.GetMultip(),LEDCount,ledTime,veto.GetTimeSec(),ledTime-veto.GetTimeSec());


			// Apply time cut
			if (fabs(LEDperiod - x_deltaT) > 5 * LEDrms) 
				TimeCut = true;

			// Force-pass low-multiplicity events with delta-t close to the LED period.
			if (!TimeCut && veto.GetMultip() <= highestMultip - 5) {
				TimeCut = true;
				force = true;
			}

			// Force-fail VERY high-multiplicity events which pass time cut.
			// This blocks the veto system from seeing any event with multiplicity
			// this high as physics.
			else if (TimeCut && veto.GetMultip() >= highestMultip - 2) {		
				TimeCut = false;
				force = true;
			}

			if (!TimeCut) {
				if (force) cout << "F: " << force << endl;
				IsLED = true;
				cout << "LEDCount++" << endl;
				LEDCount++;
				prevLED = veto;
			}
			prev_eTime = eTime;
			*/


		}
		// printf("Found %i LED's.\n",LEDCount);
		
		// Make some plots to check accuracy
		//
		TGraph *g1 = new TGraph(times.size(),&(times[0]),&(diffs1[0]));
		g1->SetTitle("scaler-sbc vs. scaler");
		g1->GetXaxis()->SetTitle("event time (seconds)");
		g1->GetYaxis()->SetTitle("time difference (seconds)");
		g1->GetYaxis()->SetTitleOffset(1.3);

		TGraph *g2 = new TGraph(times.size(),&(times[0]),&(diffs2[0]));
		g2->SetTitle("(scaler-sbc) - previous(scaler-sbc) vs. scaler");
		g2->GetXaxis()->SetTitle("sec");

		TGraph *g3 = new TGraph(times.size(),&(times[0]),&(diffs3[0]));
		g3->SetTitle("scaler-eTime vs. scaler");
		g3->GetXaxis()->SetTitle("sec");
		// g3->SetLineColor(kRed);

		TGraph *g4 = new TGraph(times.size(),&(times[0]),&(diffs4[0]));
		g4->SetTitle("scaler-pureLEDTime vs. scaler");
		g4->GetXaxis()->SetTitle("sec");


		TCanvas *c2 = new TCanvas("c2","Bob Ross's Canvas",2500,1500);
		c2->Divide(2,2);
		c2->cd(1);
		g1->Draw();
		c2->cd(2);
		g2->Draw();
		c2->cd(3);
		g3->Draw();
		c2->cd(4);
		g4->Draw();

		TCanvas *c3 = new TCanvas("c3","Bob Ross's Canvas",2500,1500);
		g1->Draw();

		char name2[200];
		sprintf(name2,"./output/times_%i.pdf",run);
		c2->Print(name2);

		sprintf(name2,"./output/timeSBC_%i.eps",run);
		c3->Print(name2);
		delete c2;
		delete c3;
		delete g1;
		delete g2;

	} // end loop over files

}