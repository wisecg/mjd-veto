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

//GatRunPlaying(Char_t FileName[],Char_t ORFileName[], Char_t OutputFile[])
//MuonTagging(Char_t FileNames[])
int MuonTagging()
{
  	//Char_t OutputFile[128];
	//strcpy (OutputFile,FileName);
	// strip off .root
	//OutputFile[strlen(OutputFile)-5]='\0';
	//strcat (OutputFile, "_Histo.root");
	
	// No output file yet
    //TFile *RootFile = new TFile(OutputFile, "RECREATE");
	
    //Open text file to store thee muon tags
    FILE * OutputFile;
    OutputFile = fopen("MuonHits.txt","w");
    
    //Read in list of root files
    Char_t TheFile[128];
    Char_t RunNumber[128];
    Char_t BaseFile[128];
    
    ifstream FileList;
    FileList.open("FileList.txt");
    while(!FileList.eof()){
        FileList >> TheFile;
        cout << "Looking in " << TheFile << endl;
        
        // strip off .root
        strcpy (BaseFile,TheFile);
        BaseFile[strlen(BaseFile)-5]='\0';
        strcpy (RunNumber,BaseFile);
        memset (RunNumber,' ',6);
        //sprintf(RunNumber,"%s",BaseFile);
        //string Run = RunNumber.substr(5);
        cout << "Run " << RunNumber << endl;
        
        
    
    
        //Read in root files
        TChain* MyChain = new TChain("vetoTree");
        MyChain->Add(TheFile);
        //cout << "Number of entries in veto chain "<< MyChain->GetEntries() << endl;
        Int_t Entries=MyChain->GetEntries();
    
    
    
        //define vectors
        vector<double>* QDC;
    
        //variables
        UInt_t card;
        Double_t Time;
        Double_t LEDTime;
    
    
    
        //address
        MyChain->SetBranchAddress("qdc",&QDC);
        MyChain->SetBranchAddress("card",&card);
    
    
    //Read in OR root files to get run info
	//TChain* ORChain = new TChain("MGTree");
    //ORChain->Add(ORFileName);
	//cout << "Number of entries in OR chain "<< ORChain->GetEntries() << endl;
    //MGTRun* MyRun = new MGTRun();
    //ORChain->SetBranchAddress("run",&MyRun);
    
    //Friend
    //MyChain->AddFriend(ORChain);
    
    //setup Histo    TH1::AddDirectory(kFALSE);
	
    
        //loop through events
    
        Bool_t GotBot1, GotBot2, GotTop, GotFile=0;
        Int_t LEDCount;
        Int_t LEDScaler=0, LEDStamp=0;
    
    
        for(Int_t i=0; i<MyChain->GetEntries(); i++){
            MyChain->GetEntry(i);
            GotBot1=0;
            GotBot2=0;
            GotTop=0;
            LEDCount=0;
            //check for card 11
            if(!(card==11)){
                continue;
            }
            //loop through channels
            for(Int_t j=0; j<QDC->size();j++){
                //look for qdc > 500
                if(QDC->at(j) > 500){
                    LEDCount++;
                    if(j< 6){
                        //got bottom1
                        GotBot1=1;
                    }
                    if(j> 5 && j<12){
                        //got bottom2
                        GotBot2=1;
                    }
                    if (j>11 && j<14)
                        //got top
                        GotTop=1;
                }
            }
            //lets check for all three, and no more than 5 hits
            if(LEDCount > 6){
                LEDScaler++;
                continue;
            }
            if(GotBot1 && GotBot2 && GotTop){
                Time=(double(i)/double(Entries)*3600);
                printf("Entry %.0f, Roughly %.2f seconds in, After %0.f LED's \n",i,Time,LEDScaler);
                if(GotFile){
                    //already have a hit in this file. Add a \n to the output
                    fprintf(OutputFile,"\n");
                }
                fprintf(OutputFile,"%s \t %.2f \t",RunNumber,Time);
                LEDStamp=LEDScaler;
                GotFile=1;
            }
        }
        //Calculate LED Time
        LEDTime=(3600/double(LEDScaler)*double(LEDStamp));
        fprintf(OutputFile,"%.2f \n",LEDTime);
        printf("Run %s had %0.f LED's \n\n",TheFile,LEDScaler);
        
        // done looking through events
        //delete chain
        delete MyChain;
        

        
    }
    FileList.close();
    fclose (OutputFile);
    return 0;
}
