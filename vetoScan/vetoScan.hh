// MJD veto analysis suite.
// Uses the November '15 MJD built data format.
// 
// Clint Wiseman, USC/Majorana
// Andrew Lopez, UTK/Majorana

#ifndef VETOSCAN_H_GUARD
#define VETOSCAN_H_GUARD

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <string>
#include "getopt.h"

#include "TFile.h"
#include "TChain.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include "TF1.h"
#include "TLine.h"

#include "MJVetoEvent.hh"
#include "GATDataSet.hh"

using namespace std;

// Processing (defined in vetoTools.cc)
void Test();
long GetStartUnixTime(GATDataSet ds);
long GetStopUnixTime(GATDataSet ds);
int GetNumFiles(string arg);
int color(int i);
int PanelMap(int i);
int* GetQDCThreshold(string file, int *arr, string name = "");
bool CheckForBadErrors(MJVetoEvent veto, int entry, int isGood, bool deactivate);
int FindQDCThreshold(TH1F *qdcHist, int panel, bool verbose);
double InterpTime(int entry, vector<double> times, vector<double> entries, vector<bool> badScaler);

// Analysis
void vetoFileCheck(string file = "", string partNum = "", bool checkBuilt = true, bool checkGat = true, bool checkGDS = false);
void vetoPerformance(string file, int *thresh = NULL, bool runBreakdowns = false);
void vetoThreshFinder(string arg, bool runHistos = false);
void muFinder(string file, int *thresh = NULL, bool root = false, bool list = false);
void muDisplayList(string file);

#endif