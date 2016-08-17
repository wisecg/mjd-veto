// muonDeadTime.cc
// Clint Wiseman, USC/Majorana
// 2/22/16
//
// Calculate the Ge dead time from a muon list.

#include "vetoScan.hh"

void durationChecker(string file)
{
	int run;
	long utc;
	double hitTime;
	int type;  // 1: "over500"		2: "vertical muon"		3:"run gap"
	bool badScaler;
	double duration = 0;
	double durationTotal = 0;

	// Input a muon list file (made by muFinder)
	ifstream InputList(file.c_str());
	if(!InputList.good()) {
    	cout << "Couldn't open " << file << endl;
    	return;
    }
	cout << "Scanning list ..." << endl;
	while(!InputList.eof())
	{
		// InputList >> run >> utc >> hitTime >> type >> badScaler;
		InputList >> run;
		GATDataSet ds(run);
		duration = (double)ds.GetRunTime()/1E9;
		durationTotal += duration;
		printf("%i  %.3f \n",run,duration);
	}
	cout << "List covers " << durationTotal << " seconds of Ge data.\n";
}

void muonDeadTime(string file)
{
	// Input a muon list file (made by muFinder)
	ifstream InputList(file.c_str());
	if(!InputList.good()) {
    	cout << "Couldn't open " << file << endl;
    	return;
    }

	int run;
	long utc;
	double hitTime;
	int type;  // 1: "over500"		2: "vertical muon"		3:"run gap"
	bool badScaler;

	// USED IN DS0: 1 second window
	// double timeBefore = .0002;	// 0.2 ms timeBefore
	// double timeAfter = 1;		// 1 sec after
	// double badScalerWindow = 8; // +/- 8 sec

	// USED IN DS1: SAME!
	double timeBefore = .0002;	// 0.2 ms timeBefore
	double timeAfter = 1;		// 1 sec after
	double badScalerWindow = 8; // +/- 8 sec

	double deadTime = 0;
	double windowOverlaps = 0;
	int woCount = 0;

	int numBadScalers = 0;
	double deadBadScalerTime = 0;

	int numGoodScalers = 0;
	double deadGoodScalerTime = 0;

	// double duration = 0;

	cout << "Scanning list ..." << endl;
	while(!InputList.eof())
	{
		InputList >> run >> utc >> hitTime >> type >> badScaler;

		if (type == 1 || type == 2) {
			if (badScaler) {
				deadTime += 2*badScalerWindow;
				deadBadScalerTime += 2*badScalerWindow;
			}
			else {
				deadTime += timeBefore + timeAfter;
				deadGoodScalerTime += timeBefore + timeAfter;
			}
		}
		if (type == 3) {
			if (badScaler) {
				deadTime += badScalerWindow;
				deadBadScalerTime += badScalerWindow;
			}
			else {
				deadTime += timeAfter;
				deadGoodScalerTime += timeAfter;
			}
		}

		if (badScaler) numBadScalers++;
		else numGoodScalers++;

		// Attempt to find the possible overlapping events.
		// These are usually at the end of a run.
		// But since "duration" is only accurate to an integer,
		// they may not be true overlaps.
		// We don't know exactly when runs end, which makes
		// finding overlaps more difficult than what is done here.
		//
		/*
		GATDataSet ds(run);
		duration = ds.GetRunTime()/CLHEP::second;
		if (!badScaler) {
			if (hitTime+timeAfter > duration) {
				windowOverlaps += hitTime+timeAfter-duration;
				woCount++;
				printf("%i  dur: %.3f  hit: %.3f  dead: %.3f  ovrlap: %.3f  ht+ta-d: %.3f\n"
					,run,duration,hitTime,deadTime,windowOverlaps,hitTime+timeAfter-duration);
			}
			else if (type != 3 && hitTime-timeBefore < 0) {
				windowOverlaps += timeBefore-hitTime;
				printf("%i  dur: %.3f  hit: %.3f  dead: %.3f  ovrlap: %.3f tb-ht: %.3f \n"
					,run,duration,hitTime,deadTime,windowOverlaps,timeBefore-hitTime);
				woCount++;
			}
		}
		else {
			if (hitTime+badScalerWindow > duration) {
				windowOverlaps += hitTime+badScalerWindow - duration;
				woCount++;
				printf("BS %i  dur: %.3f  hit: %.3f  dead: %.3f  ovrlap: %.3f  ht+bsw-d: %.3f\n"
					,run,duration,hitTime,deadTime,windowOverlaps,hitTime+badScalerWindow-duration);
			}
			if (type != 3 && hitTime-badScalerWindow < 0) {
				windowOverlaps += badScalerWindow-hitTime;
				woCount++;
				printf("BS %i  dur: %.3f  hit: %.3f  dead: %.3f  ovrlap: %.3f bsw-ht: %.3f\n"
					,run,duration,hitTime,deadTime,windowOverlaps,badScalerWindow-hitTime);
			}
		}
		*/

	}
	cout << "Dead time due to veto: " << deadTime << " seconds." << endl;

	printf("Bad scalers: %i, %.2f of %.2f sec (%.2f%%)\n",numBadScalers,deadBadScalerTime,deadTime,((double)deadBadScalerTime/deadTime)*100);

	printf("Good scalers: %i, %.2f of %.2f sec (%.2f%%)\n",numGoodScalers,deadGoodScalerTime,deadTime,((double)deadGoodScalerTime/deadTime)*100);

	cout << "Adjustment from " << woCount << " window overlaps: " << windowOverlaps << " seconds" << endl;

	double newDeadTime=deadTime-windowOverlaps;
	cout << "Actual dead time due to veto: " << newDeadTime << " seconds" << endl;
}

