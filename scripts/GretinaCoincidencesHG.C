#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TF1.h>
#include <TGraph.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>
#include <stdio.h>
using namespace std;

/*
  TO DO:
  1. Add in coincident event (2sec) tree for list display
  2. Change pretty macros to go to 10 MeV and change binning so they still look pretty.
*/

int GretinaCoincidencesHG()
{  
  int pdsf = 1;				// switch: 0 to analyze local files, 1 to run on pdsf files
  
  // Cut parameters
  double clock = 100000000;
  double loWindow = 360*clock;		// pulser duration: 300 (1st 5 mins) until runs in Jan 2015
  double hiWindow = 3601*clock; 	// Approx. max time (taken from run 45001058): 360096289127
  double e_cut = 50; 			// keV, eliminate low-energy noise

  // Input file (run numbers and muon times)
  // Output file (event information: overflow counters, coincident event info.)
  FILE *input;
  char params[200];
  input = fopen("MuonHitsFinal.txt","r");
  if(input == NULL) {
    perror("Error opening file");
    return(-1);
  }
  
  char outputName[100]="MuonHitsFinalHG";
  char outputRoot[100], outputNoCuts[100], outputCut2[100], outputCut3[100];
  sprintf(outputRoot,"%s.root",outputName);  
  sprintf(outputNoCuts,"%sNoCuts.C",outputName);  
  sprintf(outputCut2,"%sCut2.C",outputName);
  sprintf(outputCut3,"%sCut3.C",outputName);
  TFile *RootFile = new TFile(outputRoot, "RECREATE");	
  TH1::AddDirectory(kFALSE); // Global flag: "When a (root) file is closed, all histograms in memory associated with this file are automatically deleted."
  TTree *highE = new TTree("highE","High-energy events");
  int runNumber = 0;      // watch for double-counting!
  int prevrunNumber = 0;  
  float lookHere = 0;
  float prevlookHere = 0; 
  double e_cal = 0;
  int runNumber = 0;
  int timeSec = 0;
  int chEasy = 0;
  double diffTime = 0;
  highE->Branch("runNumber",&runNumber,"runNumber/I");
  highE->Branch("timeSec",&timeSec,"timeSec/I");   
  highE->Branch("e_cal",&e_cal,"e_cal/D");
  highE->Branch("chEasy",&chEasy,"chEasy/I");
  highE->Branch("diffTime",&diffTime,"diffTime/D");
  TTree *coinEvents4 = new TTree("coinEvents4","Events within 2sec of muon hits");
  coinEvents4->Branch("runNumber",&runNumber,"runNumber/I");
  coinEvents4->Branch("timeSec",&timeSec,"timeSec/I");   
  coinEvents4->Branch("e_cal",&e_cal,"e_cal/D");
  coinEvents4->Branch("chEasy",&chEasy,"chEasy/I");
  coinEvents4->Branch("diffTime",&diffTime,"diffTime/D");

  // Energy histograms:
  // "Full spectrum" only cuts out events > 50keV with good timestamps. No muon timing cuts.
  TH1D *fullSpectrum = new TH1D("fullSpectrum","Total Energy Spectrum (Unique Events)",2000,0,10000);//"bins,lower,upper"
  fullSpectrum->GetXaxis()->SetTitle("Energy [KeV]");
  // 3K spectrum - for 2D histo
  TH1D *Spectrum3k = new TH1D("Spectrum3k","Total Energy Spectrum",3000,0,3000);//"bins,lower,upper"
  Spectrum3k->GetXaxis()->SetTitle("Offline Energy [KeV]");
  // events from -10 sec before muon event to +60 sec after (asymmetrical)
  TH1D *energyCoins70 = new TH1D("energyCoins70","Energy of Ge events near muon hits (-10s to 60s)",2000,0,10000); 
  energyCoins70->GetXaxis()->SetTitle("Offline Energy [KeV]");
  // events inside +/- 100 sec of muon hits (symmetrical)
  TH1D *energyCoins200 = new TH1D("energyCoins200","Energy of Ge events near muon hits (+/- 100sec)",2000,0,10000); 
  energyCoins200->GetXaxis()->SetTitle("Offline Energy [KeV]");
  // events inside +/- 2 sec of muon hits (symmetrical)
  TH1D *energyCoins4 = new TH1D("energyCoins4","Energy of Ge events near muon hits (+/- 2sec)",2000,0,10000);
  energyCoins4->GetXaxis()->SetTitle("Offline Energy [KeV]");
  // Events in each channel, per run (no timing cuts)
  char ename[100]; 
  TH1D *energyByChannel[7]; // call from 0 to 6, maps to chEasy.
  for (int i=0;i<7;i++) {
    sprintf(ename,"ch%dEnergy",i);
    energyByChannel[i]= new TH1D(ename,"",10000,0,10000);
    energyByChannel[i]->GetXaxis()->SetTitle("45000000 + Run Number");
  }
    
  // Time histogram:
  // events inside +/- 100 sec of muon hits (symmetrical)
  TH1D *timeCoins200 = new TH1D("timeCoins200","Time diff between muon and Ge events (+/- 100sec)", 200,-100, 100);
  timeCoins200->GetXaxis()->SetTitle("Time of Ge event relative to muon hit [sec]");  
  // events from -10 sec before muon event to +60 sec after (asymmetrical)  
  TH1D *timeCoins70 = new TH1D("timeCoins70","Time diff between muon and Ge events (-10s to 60s)", 70,-10, 60);
  timeCoins70->GetXaxis()->SetTitle("Time of Ge event relative to muon hit [sec]");  
  // events inside +/- 2 sec of muon hits (symmetrical)
  TH1D *timeCoins4 = new TH1D("timeCoins4","Time diff between muon and Ge events (+/- 2sec)",40,-2,2);
  timeCoins4->GetXaxis()->SetTitle("Time of Ge event relative to muon hit [sec]");

  // Counter histograms (filled based on unique events)
  TH1D *eventsByChannel = new TH1D("eventsByChannel","Total events, LG channels (arr. by relative height)",6,0,6);
  eventsByChannel->GetXaxis()->SetTitle("Channel number (chEasy)");
  TH1D *eventsByChannelCut2 = new TH1D("eventsByChannelCut2","Total events, LG channels (arr. by relative height)",6,0,6);
  eventsByChannelCut2->GetXaxis()->SetTitle("Channel number (chEasy)");
  TH1D *numDetectorsHit = new TH1D("numDetectorsHit","Multiple-hits in detectors",6,1,7);
  numDetectorsHit->GetXaxis()->SetTitle("Number of detectors hit");
  // Total unique events, per run (no timing cuts)
  TH1D *eventsPerRun = new TH1D("eventsPerRun","Total events per run",2000,0,2000);   // depends on range of runs in vetoBkgds
  eventsPerRun->GetXaxis()->SetTitle("45000000 + Run Number");
  // Events in each channel, per run (no timing cuts)
  char hname[100]; 
  TH1D *rateByChannel[7]; // call from 0 to 6, maps to chEasy.
  for (int i=0;i<7;i++) {
    sprintf(hname,"ch%dRate",i);
    rateByChannel[i]= new TH1D(hname,"",2000,0,2000);
    rateByChannel[i]->GetXaxis()->SetTitle("45000000 + Run Number");
  }

  // Channel (y-axis) vs. Energy (x-axis) - BEFORE muon cuts
  TH2F *chanVsEnergy = new TH2F("chanVsEnergy","Channel vs. Energy (no timing cut)",100,0,10000,6,0,6);
  chanVsEnergy->GetXaxis()->SetTitle("Offline Energy [KeV]");
  chanVsEnergy->GetYaxis()->SetTitle("channelEasy");

  // Channel vs. number of detectors hit -  BEFORE muon cuts
  TH2F *chanVsNumHit = new TH2F("chanVsNumHit","Channel vs. Num. Detectors Hit (no timing cut)",6,0,6,6,0,6);
  chanVsNumHit->GetXaxis()->SetTitle("Number of Detectors Hit");
  chanVsNumHit->GetXaxis()->SetTitle("channelEasy");

  // Channel (y-axis) vs. Energy (x-axis) - AFTER muon cuts
  TH2F *cutChanVsEnergy = new TH2F("cutChanVsEnergy","Channel vs. Energy (+/- 2sec cut)",100,0,10000,6,0,6);
  cutChanVsEnergy->GetXaxis()->SetTitle("Offline Energy [KeV]");
  cutChanVsEnergy->GetYaxis()->SetTitle("channelEasy");

  // Channel vs. number of detectors hit -  BEFORE muon cuts
  TH2F *cutChanVsNumHit = new TH2F("cutChanVsNumHit","Channel vs. Num. Detectors Hit (+/- 2sec cut)",6,0,6,6,0,6);
  cutChanVsNumHit->GetXaxis()->SetTitle("Number of Detectors Hit");
  cutChanVsNumHit->GetXaxis()->SetTitle("channelEasy");  

  // Scan over list of runs
  while( fgets(params, 200, input)!=NULL ) {
    int hits = 0;		
    int counter200 = 0;   
    int counter70 = 0;
    int counter4 = 0;
    int c[18] = {0};		  // counts how many detectors are hit per event.
    int d[18] = {0};      // counts hits in individual channels 
    int chan2[18] = {0};  // counts hits in individual channels, in the +/- 2 second window.
    int overflow = 0;     // overflow (events >10 MeV)
    
    lookHere = 0;  	
    sscanf(params, "%d %e", &runNumber, &lookHere);
    if (GoodRunCheck(runNumber,lookHere)==0) {
      cout << "Failed GoodRunCheck, run:" << runNumber << endl;
      continue;
    } 
    printf("Now scanning run %d (muon at %dsec) ...... ",runNumber,lookHere);
    
    // time between muon hits (used to prevent double counting)
    double timeDiff = lookHere-prevlookHere;

    // for the "200" histograms
    double l = (lookHere-100)*clock;
    double u = (lookHere+100)*clock;
    double pl = (prevlookHere-100)*clock;
    double pu = (prevlookHere+100)*clock;

    // for the "4" histograms
    double sl = (lookHere-2)*clock;
    double su = (lookHere+2)*clock;
    double spl = (prevlookHere-2)*clock;
    double spu = (prevlookHere+2)*clock;

    // for the "70" histograms
    double el = (lookHere-10)*clock;
    double eu = (lookHere+60)*clock;
    double epl = (prevlookHere-10)*clock;
    double epu = (prevlookHere+60)*clock;    
   
    TChain* t1 = new TChain("mjdTree");  	// put this inside loop to count up events separately.
    vector<double>* energyCal = 0;
    vector<double>* timestamp = 0;
    vector<double>* channel = 0;
    vector<double>* trap4usMax = 0;
    vector<double>* energy = 0;
    t1->SetBranchAddress("energyCal",&energyCal);
    t1->SetBranchAddress("timestamp",&timestamp);
    t1->SetBranchAddress("channel",&channel);
    t1->SetBranchAddress("trap4usMax",&trap4usMax);
    t1->SetBranchAddress("energy",&energy);
    
    char infile[200], infilename[200];
    sprintf(infile,"mjd_run%d",runNumber);  
    if (pdsf==1) sprintf(infilename,"/global/project/projectdirs/majorana/data/mjd/surfprot/data/gatified/P3END/%s.root",infile);
    else if (pdsf==0) sprintf(infilename,"%s.root",infile);
    t1->Add(infilename);   
    //t1->Add("mjd_run45001058.root");		// look at single file directly
    printf("%u entries.\n",t1->GetEntries());

    if (runNumber==prevrunNumber) 
      printf("Warning! re-analyzing same file for different muon event.\n");
           
     // Scan run 
    Int_t nentries = t1->GetEntries();
    for (Int_t i = 0; i<nentries;i++) {
      t1->GetEntry(i); 
      int n = channel->size(); 
          
      // Loop over channels with nonzero entries
      for(int j=0; j<n; j++) {
        int ch = channel->at(j);
        chEasy = mapchannel(ch);
      	double trap4 = trap4usMax->at(j);
        double e_kev = energyCal->at(j);       
      	double e_raw = energy->at(j);
      	
        if (runNumber>=45001829 && runNumber<=45002768)     // get calibrated energy (Wenqin's email "offline energy")
          e_cal = (trap4*e_kev/e_raw/0.0022255);
        else
          e_cal = (trap4*e_kev/e_raw/0.002494);
        

        double time = timestamp->at(j);
        timeSec = time/clock;
      	diffTime = time/clock-lookHere;

      	if (chEasy != 1000 && e_cal>e_cut && time>loWindow && time<hiWindow) {
          
          /* Fill overall energy spectrum (no time cuts) without double counting.
           * This only outputs unique events, which are the ones checked for double counting. */
          if(runNumber!=prevrunNumber) {  
            //printf("fullSpectrum data: chEasy:%u  e_kev:%.1f  e_cal:%.1f  time:%.1fsec  diffTime:%.1f  abs(timeDiff):%.1f \n",chEasy,e_kev,e_cal,time/clock,diffTime);
            fullSpectrum->Fill(e_cal);
            Spectrum3k->Fill(e_cal);
            chanVsEnergy->Fill(e_cal,chEasy); // x, y
            chanVsNumHit->Fill(chEasy,j);
            energyByChannel[chEasy]->Fill(e_cal);
      	    hits++; c[j]++; d[chEasy]++; 
            if (e_cal>=10000) {
              printf("Overflow event: %.1f KeV  run:%u  time:%.5f  detector:%u  diffTime:%.1f \n",e_cal,runNumber,time/clock,chEasy,diffTime);
              overflow++;
            }
            if (e_cal>=3000) {
              printf("High-energy event: %.1f KeV  run:%u  time:%.5f  detector:%u  diffTime:%.1f \n",e_cal,runNumber,time/clock,chEasy,diffTime);
              highE->Fill();
            }
          }
        
          /* Fill the "200" histograms without double counting.
           * WARNING! This code will only compare to the previous event.  If there is an overlapping event from (say) two events
           * ago, it will be missed.  Fortunately this should be extremely rare for muon events, since they're generally far apart.
           * For a more permanent solution, perhaps a TEventList should be utilized.  */
          if (runNumber!=prevrunNumber && time>l && time<u) {         // normal scenario
            //printf(" timeCoins200 filled normally\n");
            timeCoins200->Fill(diffTime);
            energyCoins200->Fill(e_cal);
            counter200++;
          }
          else if (runNumber==prevrunNumber && abs(timeDiff)>200 && time>l && time<u) {   // no overlap (events are far apart)
              //printf("timeCoins200 lookHere:%.1f  prevlH:%.1f  timeDiff:%.0fsec  time: %.0f  l:%.0f  u:%.0f  pl:%.0f  pu:%.0f, No Overlap Case filled\n", lookHere,prevlookHere,timeDiff,time/clock,l/clock,u/clock,pl/clock,pu/clock);
              timeCoins200->Fill(diffTime);               
              energyCoins200->Fill(e_cal);
              counter200++;
          }
          else if (runNumber==prevrunNumber && abs(timeDiff)<=200) {  // overlap exists, two cases possible.
            if (timeDiff > 0 && time<u && time>pu) {            // common.  prevlookHere is earlier than lookHere.
              timeCoins200->Fill(diffTime);                     // overlap region is l < time < pu, don't fill here. Fill instead pu < time < u.
              energyCoins200->Fill(e_cal);
              counter200++;
            }
            else if (timeDiff < 0 && time>l && time<pl)) {      // overlap, case a, not common (probably only a glitch.)
              timeCoins200->Fill(diffTime);                     // overlap region is pl < time < u, don't fill here. Fill instead l < time < pl.
              energyCoins200->Fill(e_cal);
              counter200++;
            }
          }

          /* Fill the "70" histograms without double counting.
           * WARNING! Same issue as above.  */
          if (runNumber!=prevrunNumber && time>el && time<eu) {         // normal scenario
            energyCoins70->Fill(e_cal);
            timeCoins70->Fill(diffTime);
            counter70++;
          }
          else if (runNumber==prevrunNumber && abs(timeDiff)>70 && time>el && time<eu) {   // no overlap (events are far apart)
              energyCoins70->Fill(e_cal);
              timeCoins70->Fill(diffTime);
              counter70++;               
          }
          else if (runNumber==prevrunNumber && abs(timeDiff)<=70) {  // overlap exists, two cases possible.
            if (timeDiff > 0 && time<eu && time>epu) {          // common.  prevlookHere is earlier than lookHere.
              energyCoins70->Fill(e_cal);                       // overlap region is el < time < epu, don't fill here.  Fill instead epu < time < eu.
              timeCoins70->Fill(diffTime);
              counter70++;
            }
            else if (timeDiff < 0 && time>el && time<epl)) {    // overlap, case a, not common.
              energyCoins70->Fill(e_cal);                       // overlap region is epl < time < eu, don't fill here.  Fill instead el < time < epl.
              timeCoins70->Fill(diffTime);
              counter70++;
            }
          }
          /* Fill the "4" histograms without double counting.  */
          if (runNumber!=prevrunNumber && time>sl && time<su) {         // normal scenario
            timeCoins4->Fill(diffTime);
            energyCoins4->Fill(e_cal);
            cutChanVsNumHit->Fill(chEasy,j);
            cutChanVsEnergy->Fill(e_cal,chEasy);
            chan2[chEasy]++;
            counter4++;
            coinEvents4->Fill();
          }
          else if (runNumber==prevrunNumber && abs(timeDiff)>4 && time>sl && time<su) {   // no overlap (events are far apart)
              timeCoins4->Fill(diffTime);               
              energyCoins4->Fill(e_cal);
              cutChanVsNumHit->Fill(chEasy,j);
              cutChanVsEnergy->Fill(e_cal,chEasy);
              chan2[chEasy]++;
              counter4++;
              coinEvents4->Fill();
          }
          else if (runNumber==prevrunNumber && abs(timeDiff)<=4) {  // overlap exists, two cases possible.
            if (timeDiff > 0 && time<su && time>spu) {            // common.  prevlookHere is earlier than lookHere.
              timeCoins4->Fill(diffTime);                     // overlap region is l < time < pu, don't fill here. Fill instead pu < time < u.
              energyCoins4->Fill(e_cal);
              cutChanVsNumHit->Fill(chEasy,j);
              cutChanVsEnergy->Fill(e_cal,chEasy);
              chan2[chEasy]++;
              counter4++;
              coinEvents4->Fill();
            }
            else if (timeDiff < 0 && time>sl && time<spl)) {      // overlap, case a, not common (probably only a glitch.)
              timeCoins4->Fill(diffTime);                     // overlap region is pl < time < u, don't fill here. Fill instead l < time < pl.
              energyCoins4->Fill(e_cal);
              cutChanVsNumHit->Fill(chEasy,j);
              cutChanVsEnergy->Fill(e_cal,chEasy);
              chan2[chEasy]++;
              counter4++;
              coinEvents4->Fill();
            }
          }
        }
      }
    }

    // Fill counter histograms
    for (int i=1;i<=7;i++) {  
      eventsByChannel->Fill(i-1,d[i-1]);  // Fill(bin,weight)  d[0]-d[6]
      eventsByChannelCut2->Fill(i-1,chan2[i-1]);
      numDetectorsHit->Fill(i,c[i-1]);    // c[1]-c[7]
      rateByChannel[i-1]->Fill(runNumber-45000000,d[i-1]);
    }
    eventsPerRun->Fill(runNumber-45000000,hits);  

    // Screen output
    printf("    Finished. Unique events in the spectrum between %.0f and %.0f seconds, above %.0f KeV: %u \n",loWindow/clock,hiWindow/clock,e_cut,hits);
    printf("      Of %u entries, %u survived (%.5f percent)\n",nentries,hits,((hits/(double)nentries)*100));
    printf("      Channel counters: ");  // watch for excessive events in a channel (noisy run)
    for (int i=0;i<=6;i++) cout << " ch" << i << "=" << d[i];
    cout << endl;
    printf("      Filled histos -- counter200: %u entries, counter70: %u entries, counter4: %u entries.\n",counter200, counter70, counter4);
    if (counter70>counter200) printf("      WARNING, 70 has more than 200!\n");
    printf("      LG Overflow Count (>10 MeV): %u \n",overflow);

    
    // File output    
    fullSpectrum->Write("", TObject::kOverwrite);
    Spectrum3k->Write("", TObject::kOverwrite);
    energyCoins4->Write("", TObject::kOverwrite);
    energyCoins70->Write("", TObject::kOverwrite);
    energyCoins200->Write("", TObject::kOverwrite);
    timeCoins4->Write("", TObject::kOverwrite);
    timeCoins70->Write("", TObject::kOverwrite);
    timeCoins200->Write("", TObject::kOverwrite);
    eventsPerRun->Write("",TObject::kOverwrite);
    eventsByChannel->Write("",TObject::kOverwrite);
    eventsByChannelCut2->Write("",TObject::kOverwrite);
    numDetectorsHit->Write("",TObject::kOverwrite);
    chanVsEnergy->Write("",TObject::kOverwrite);
    cutChanVsEnergy->Write("",TObject::kOverwrite);
    chanVsNumHit->Write("",TObject::kOverwrite);
    cutChanVsNumHit->Write("",TObject::kOverwrite);
    for (int i=0;i<7;i++) rateByChannel[i]->Write("",TObject::kOverwrite);
    for (int i=0;i<7;i++) energyByChannel[i]->Write("",TObject::kOverwrite);

    // Delete chain for memory management (required for PDSF operation)
    delete t1;              

    // remember this run number and muon hit time
    prevrunNumber=runNumber;
    prevlookHere=lookHere;
  }

  // Draw some FANCY composite plots.

  // Unique events, no timing cuts.
  TCanvas *c1 = new TCanvas("c1", "Bob Ross's Canvas",900,900);
  //gStyle->SetOptStat(0);  //false
  gStyle->SetPalette(1);  //true
  TPad *centerPad = new TPad("centerPad", "centerPad",0.0,0.0,0.65,0.6);  //xlow, ylow, xup, yup
  centerPad->Draw();
  rightPad = new TPad("rightPad", "rightPad",0.65,0.0,1.0,0.6);
  rightPad->Draw();
  botPad = new TPad("botPad", "botPad",0.0,0.55,0.65,1.0);
  botPad->Draw();
  centerPad->cd();
  chanVsEnergy->SetFillColor(kBlue+1);
  chanVsEnergy->Draw("COLZ");
  rightPad->cd();
  eventsByChannel->SetFillColor(kBlue-2);
  eventsByChannel->Draw("hbar");
  botPad->cd();
  botPad->SetLogy();
  //Spectrum3k->SetFillColor(kBlue);
  //Spectrum3k->Draw("bar");
  fullSpectrum->Draw();
  c1->Print(outputNoCuts);

  // Unique events, in the +/- 2 second window
  TCanvas *c2 = new TCanvas("c2", "Bob Ross's Other Canvas",900,900);
  //gStyle->SetOptStat(0);  //false
  gStyle->SetPalette(1);  //true
  TPad *p1 = new TPad("p1","p1",0.0,0.0,0.65,0.6);  //xlow, ylow, xup, yup
  p1->Draw();
  p2 = new TPad("p2","p2",0.65,0.0,1.0,0.6);
  p2->Draw();
  p3 = new TPad("p3","p3",0.0,0.55,0.65,1.0);
  p3->Draw();
  p4 = new TPad("p4","p4",0.65,0.6,1.0,1.0);
  p4->Draw();
  p1->cd();
  cutChanVsEnergy->SetFillColor(kBlue+1);
  cutChanVsEnergy->Draw("COLZ");
  p2->cd();
  eventsByChannelCut2->SetFillColor(kBlue-2);
  eventsByChannelCut2->Draw("hbar");
  p3->cd();
  p3->SetLogy();
  energyCoins4->SetFillColor(kBlue);
  energyCoins4->Draw("bar");
  p4->cd();
  timeCoins4->Draw();
  c2->Print(outputCut2);

  // Hits in channels vs number of detectors hit.
  TCanvas *c3 = new TCanvas("c3","Bob Ross's Other Other Canvas",1400,800);
  TPad *leftPad2 = new TPad("leftPad2","leftPad2",0.0,0.0,0.5,1.0);
  leftPad2->Draw();
  TPad *rightPad2 = new TPad("rightPad2","rightPad2",0.5,0.0,1.0,1.0);
  rightPad2->Draw();
  leftPad2->cd();
  chanVsNumHit->Draw("COLZ");
  rightPad2->cd();
  cutChanVsNumHit->Draw("COLZ");
  c3->Print(outputCut3);

  fclose(input);
  highE->Write();
  coinEvents4->Write();
  RootFile->Close();
  cout << " Wrote root file.\n";
  
  return;
  }

// Susanne's channel map function, modified by Clint.
int mapchannel(int ch)
{
    int chEasy = 1000;
    
    // Susanne's original mapping
    /*if(ch == 88) chEasy = 0;  //So1 HG
    if(ch == 89) chEasy = 1;  //So1 LG
    if(ch == 146) chEasy = 2; //P2 HG
    if(ch == 147) chEasy = 3; //P2 LG
    if(ch == 148) chEasy = 4; //P1 HG
    if(ch == 149) chEasy = 5; //P1 LG
    if(ch == 150) chEasy = 6; //D2 HG
    if(ch == 151) chEasy = 7; //D2 LG
    if(ch == 152) chEasy = 8; //D1 HG
    if(ch == 153) chEasy = 9; //D1 LG
    */
    
    // List on Vince's wall
    // See WenqinDetectorByDetector.pdf for more information
    // S3D3 is dead
    /*if(ch == 150) chEasy = 0; //S1D1  HG  (unstable gain shifts)
    if(ch == 151) chEasy = 1; //  LG  (unstable gain shifts)
    if(ch == 148) chEasy = 2; //S1D2  HG
    if(ch == 149) chEasy = 3; //  LG
    if(ch == 146) chEasy = 4; //S1D3  HG
    if(ch == 147) chEasy = 5; //  LG
    if(ch == 144) chEasy = 6; //S1D4  HG
    if(ch == 145) chEasy = 7; //  LG
    if(ch == 152) chEasy = 8; //S2D1  HG
    if(ch == 153) chEasy = 9; //  LG
    if(ch == 120) chEasy = 10;  //S3D1  HG
    if(ch == 121) chEasy = 12;  //  LG
    if(ch == 118) chEasy = 12;  //S3D2  HG
    if(ch == 119) chEasy = 13;  //  LG
    if(ch == 114) chEasy = 14;  //S3D4  HG
    if(ch == 115) chEasy = 15;  //  LG
    if(ch == 112) chEasy = 16;  //S3D5  HG
    if(ch == 113) chEasy = 17;  //  LG
    */
    /*
    // Choose only low-gain channels for muon events
    // arrange by relative height
    //if(ch == 153) chEasy = 3; // noisy! 
    //if(ch == 151) chEasy = 0; // noisy! 
    if(ch == 121) chEasy = 6;
    if(ch == 119) chEasy = 5;
    if(ch == 149) chEasy = 4; 
    if(ch == 147) chEasy = 3; 
    //if(ch == 115) chEasy = 2; // producing weird peak at ~55 keV
    if(ch == 145) chEasy = 1;
    if(ch == 113) chEasy = 0;
    */

    // try high-gain channels to see if the 115 problem goes away.
    if(ch == 120) chEasy = 6;
    if(ch == 118) chEasy = 5;
    if(ch == 148) chEasy = 4; 
    if(ch == 146) chEasy = 3; 
    if(ch == 114) chEasy = 2; 
    if(ch == 144) chEasy = 1;
    if(ch == 112) chEasy = 0;

    
    return chEasy;
}

// returns 0 for bad run, 1 for good run
int GoodRunCheck(int run, float time) {
  // Good run ranges (Wenqin's list, email: "runs in my analysis")
  int goodrun_start[100]={45000509, 45000950, 45001458, 45001597, 45002191, 45002462, 45002816, 45003847, 45004179, 45004508, 0};
  int goodrun_end[100]={45000600, 45001100, 45001571, 45001821, 45002449, 45002768, 45003830, 45004137, 45004495, 45004709, 0};
  // Bad runs 
  int badrun_start[100]={45000261, 45000455, 45000901, 45001318, 45001944, 45002184, 45002319, 45002416, 45002796, 45002990, 0};
  int badrun_end[100]={45000283, 45000458, 45000919, 45001415, 45001962, 45002190, 45002321, 45002418, 45002815, 45003001, 0};
  
  int runStatus = 0; 
  for(int i=0;i<10;i++) {
    if (run <= goodrun_end[i] && run >= goodrun_start[i] && time>=360) {
      runStatus = 1;
      break;
    }
    if (run <= badrun_end[i] && run >= badrun_start[i]) {
      runStatus = 0;
      cout << "Found bad run.\n";
      break;
    }
  }
  if (runStatus==1) return 1;
  else return 0;
}