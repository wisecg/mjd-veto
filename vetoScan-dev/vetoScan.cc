// MJD veto analysis suite.
// Uses the November '15 MJD built data format.
// 
// Clint Wiseman, USC/Majorana
// Andrew Lopez, UTK/Majorana

#include "vetoScan.hh"

using namespace std;

static const char Usage[] =
"\n"
"Usage: vetoScan [options].\n"
"     REQUIRED for most routines : \n"
"     -F (--file) ./path/to/your/runList.txt\n\n"
"Additional options:\n"
"     -h (--help) : Print usage info\n"
"     -S (--serial) : Set the part number (P3JDY, etc.)  REQUIRED to use checkFiles.\n"
"     -f (--checkFiles) : Check that files exist.\n"
"                       : Options: `checkBuilt`, `checkGAT`, `checkBoth`, `checkGDS`, `checkAll`\n" 
"                       : (checkBuilt and checkGAT both require PDSF)\n"
"     -H (--findThresh) : Find QDC software thresholds for a set of runs.\n"
"                       : Options: `runs` or `totals`\n"
"     -T (--swThresh) : Set QDC software threshold using `vetoSWThresholds.txt`\n"
"     -m (--muFinder) : Scan runs for muons.\n"
"                     : If -T is specified, user picks which SW thresholds to use.\n"
"                     : Output options: `root`,`list`,`both`\n"
"     -p (--perfCheck) : Veto performance check (data quality).\n"
"                      : Option: `runs`, `totals`\n"
"                      : If -T is specified, user picks which SW thresholds to use.\n"
"     -t (--timeCheck) : Veto time variables check.\n"
"     -u (--duration) : Find duration (in seconds) of file list.\n"
"     -l (--findLED) : Find veto LED events.\n"
"     -d (--dead) : Calculate Ge dead time from a muon list.\n"
"     -o (--plot) : Run muPlotter\n"
"     -r (--parse) : Run muParser\n"
"     -G (--geCoins) : run muGeCoins\n"
"     -D (--dispList) : Create veto hit list for vetoDisplay code\n"
"     -L (--vetoList) : Create veto hit list for DEMONSTRATOR Veto Cut\n"
"     -s (--muSimple) : Run a simplified version of muFinder\n"
"\n";

int main(int argc, char** argv) 
{	
	if (argc == 1) {
		cout << Usage;
		return 0;
	}

	// Parse command line arguments with getopt_long:
	// http://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Option-Example.html
	//
	string file = "", partNum = "", threshName = "";
	bool findMuons=0, perfCheck=0, fileCheck=0, findTime=0,findLED=0,findThresh=0,deadTime=0,durationCheck=0;
	bool muPlot=0, muParse=0,checkBuilt=0,checkGAT=0,checkGDS=0,root=0,list=0;
	bool runBreakdowns=0,geCoins=0,muList=0,vetoCutList=0;
	bool muSimp=0;
	//
	int c;
	int option_index = 0;
	while (1)
	{
		static struct option long_options[] = {
			{"help", no_argument, 0, 'h'},
			{"file", no_argument, 0, 'F'},
			{"serial",no_argument, 0, 'S'},
			{"checkFiles", required_argument, 0, 'f'},
			{"findThresh", no_argument, 0, 'H'},
			{"swThresh", required_argument, 0, 'T'},
			{"muFinder", required_argument, 0, 'm'},
			{"perfCheck", required_argument, 0, 'p'},
			{"timeCheck", no_argument, 0, 't'},
			{"duration", no_argument, 0, 'u'},
			{"findLED", no_argument, 0, 'l'},
			{"dead", no_argument, 0, 'd'},
			{"plot", no_argument, 0, 'o'},
			{"parse", no_argument, 0, 'r'},
			{"geCoins", required_argument, 0, 'G'},
			{"dispList", no_argument,0,'D'},
			{"vetoList", no_argument, 0, 'L'},
			{"muSimple", no_argument, 0, 's'}
		};

		// don't forget to add a new option here too!
		c = getopt_long (argc, argv, "hF:S:f:H:T:m:p:tldorGDLsu",long_options,&option_index);
		if (c == -1) break;

		switch (c)
		{
		case 'h':
			cout << Usage;
			break;
		case 'F':
			file = string(optarg);
			cout << "Using input file: " << file << endl;
			break;
		case 'S':
			partNum = string(optarg);
			cout << "Using part number: " << partNum << endl;
			break;
		case 'f':
			fileCheck=1;
			if (string(optarg) == "checkBuilt") 		checkBuilt=1;
			else if (string(optarg) == "checkGAT")		checkGAT=1;
			else if (string(optarg) == "checkGDS")		checkGDS=1;
			else if (string(optarg) == "checkBoth") {	checkBuilt=1; checkGAT=1; }
			else if (string(optarg) == "checkAll")  {	checkBuilt=1; checkGAT=1; checkGDS=1; }
			printf("Checking for ... built? %i  gatified? %i  GATDataSet? %i\n",checkBuilt,checkGAT,checkGDS);
			break;
		case 'H':
			findThresh=1;
			if (string(optarg) == "runs") runBreakdowns=1;
	    	else if (string(optarg) == "totals") runBreakdowns=0;
			break;
		case 'T':
			threshName = string(optarg);
			cout << "Using thresholds from " << threshName << endl;
			break;
		case 'm': 
			findMuons=1; 
			if (string(optarg) == "root") root=1;
			else if (string(optarg) == "list") list=1;
			else if (string(optarg) == "both") { root=1; list=1; }
			break;
	    case 'p': 
	    	perfCheck=1; 
	    	if (string(optarg) == "runs") runBreakdowns=1;
	    	else if (string(optarg) == "totals") runBreakdowns=0;
	    	break;
		case 't': findTime=1; break;
		case 'l': findLED=1; break;
		case 'd': deadTime=1; break;
		case 'u': durationCheck=1; break;
		case 'o': muPlot=1; break;
		case 'r': muParse=1; break;
		case 'G': geCoins=1; break;
		case 'D': muList=1; break;
		case 'L': vetoCutList=1; break;
		case 's': muSimp=1; break;
		case '?':
		    if (isprint (optopt))  fprintf (stderr, "Unknown option `-%c'.\n", optopt);
		    else fprintf (stderr,"Unknown option character `\\x%x'.\n",optopt);
		    cout << "Call vetoScan -h for usage options.\n";
		    return 1;
		default:
		    abort ();
		}
	}

	// =======================================================
	// Run selected routines
	//
	int thresh[32] = {0};

	if (fileCheck) 	vetoFileCheck(file,partNum,checkBuilt,checkGAT,checkGDS);
	if (findTime)	vetoTimeFinder(file);
	if (findThresh)	vetoThreshFinder(file,runBreakdowns);
	if (perfCheck)	
	{	
		if (threshName != "") GetQDCThreshold(file,thresh,threshName);
		else GetQDCThreshold(file,thresh);
		vetoPerformance(file,thresh,runBreakdowns);
	}
	if (findMuons) 
	{  	
		if (threshName != "") GetQDCThreshold(file,thresh,threshName);
		else GetQDCThreshold(file,thresh);
		muFinder(file,thresh,root,list);
	}
	if (muSimp) 
	{  	
		if (threshName != "") GetQDCThreshold(file,thresh,threshName);
		else GetQDCThreshold(file,thresh);
		muSimple(file,thresh);
	}
	if (findLED) 	vetoLEDFinder(file);
	if (deadTime) 	muonDeadTime(file);
	if (durationCheck) durationChecker(file);
	if (muPlot)		muPlotter(file);
	if (muParse)	muParser(file);
	if (geCoins)	muGeCoins(file);
	if (muList)		muDisplayList(file);
	if (vetoCutList) muListGen(file);

	// =======================================================

	cout << "\nCletus codes good." << endl;
}