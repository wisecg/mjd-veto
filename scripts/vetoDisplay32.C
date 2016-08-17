/*
	vetoDisplay32.C
	Clint Wiseman & David Tedeschi, USC/Majorana
	August 2015.

	Illustrates a muon candidate event (input from text file) to help with 
	visualizing all the possible coincidences between 32 veto panels.

	Ver. 3: CW (Aug. 2015) -- added final 8 veto panels and revised channel map.
	Ver. 2: CW (Apr. 2015) -- streamlined original code, programmed input file and dynamic coloring
	Ver. 1: David J Tedeschi, USC (Mar. 2015)

	Usage: 
	root[0] .X vetoDisplay32.C++
*/

#include <iostream>
#include <TROOT.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoManager.h>
#include <TPaveText.h>

using namespace std;

// source: http://www.perbang.dk/rgbgradient/
// need to update this with new software thresholds that color under 500.
int coloring(int qdc,int mode){

	int v[15]={0};

	// x^2 coloring (separates low-qdc events better) 
	// 15.556 x^2 + 500, x=1,15
	if (mode == 1)
		v[0]=500,  v[1]=562,  v[2]=640,   v[3]=749,   v[4]=889,   v[5]=1060,  v[6]=1262,  v[7]=1496,
		v[8]=1760, v[9]=2056, v[10]=2382, v[11]=2740, v[12]=3129, v[13]=3549, v[14]=4000;  

	// log coloring (separates high-qdc events better)
	// range: 1477*log(i), scaled for i=15 == 4000
	if (mode == 2)
		v[0]=500,  v[1]=1024, v[2]=1623,  v[3]=2048,  v[4]=2377,  v[5]=2646,  v[6]=2874,  v[7]=3071,
		v[8]=3245, v[9]=3401, v[10]=3542, v[11]=3670, v[12]=3788, v[13]=3898, v[14]=4000;

	// linear coloring (separate qdc events linearly)
	if (mode == 3)
		v[0]=500,  v[1]=750,  v[2]=1000,  v[3]=1250,  v[4]=1500,  v[5]=1750,  v[6]=2000,  v[7]=2250,
		v[8]=2500, v[9]=2750, v[10]=3000, v[11]=3250, v[12]=3500, v[13]=3750, v[14]=4000;

	// set color table
	if (qdc < v[0])                 { TColor *color = gROOT->GetColor(1000);  color->SetRGB(0.00, 0.00, 0.00);  return 1000;}   // black (below threshold)
	if (qdc > v[0] && qdc <=v[1])   { TColor *color = gROOT->GetColor(1001);  color->SetRGB(0.02, 0.00, 0.75);  return 1001;}   // blue side
	if (qdc > v[1] && qdc <=v[2])   { TColor *color = gROOT->GetColor(1002);  color->SetRGB(0.08, 0.00, 0.69);  return 1002;}
	if (qdc > v[2] && qdc <=v[3])   { TColor *color = gROOT->GetColor(1003);  color->SetRGB(0.15, 0.00, 0.64);  return 1003;}
	if (qdc > v[3] && qdc <=v[4])   { TColor *color = gROOT->GetColor(1004);  color->SetRGB(0.21, 0.00, 0.59);  return 1004;}
	if (qdc > v[4] && qdc <=v[5])   { TColor *color = gROOT->GetColor(1005);  color->SetRGB(0.27, 0.01, 0.53);  return 1005;}
	if (qdc > v[5] && qdc <=v[6])   { TColor *color = gROOT->GetColor(1006);  color->SetRGB(0.33, 0.01, 0.48);  return 1006;}
	if (qdc > v[6] && qdc <=v[7])   { TColor *color = gROOT->GetColor(1007);  color->SetRGB(0.40, 0.01, 0.43);  return 1007;}
	if (qdc > v[7] && qdc <=v[8])   { TColor *color = gROOT->GetColor(1008);  color->SetRGB(0.46, 0.01, 0.37);  return 1008;}
	if (qdc > v[8] && qdc <=v[9])   { TColor *color = gROOT->GetColor(1009);  color->SetRGB(0.52, 0.02, 0.32);  return 1009;}
	if (qdc > v[9] && qdc <=v[10])  { TColor *color = gROOT->GetColor(1010);  color->SetRGB(0.58, 0.02, 0.27);  return 1010;}
	if (qdc > v[10] && qdc <=v[11]) { TColor *color = gROOT->GetColor(1011);  color->SetRGB(0.65, 0.02, 0.21);  return 1011;}
	if (qdc > v[11] && qdc <=v[12]) { TColor *color = gROOT->GetColor(1012);  color->SetRGB(0.71, 0.02, 0.16);  return 1012;}
	if (qdc > v[12] && qdc <=v[13]) { TColor *color = gROOT->GetColor(1013);  color->SetRGB(0.77, 0.02, 0.11);  return 1013;}
	if (qdc > v[13] && qdc <=v[14]) { TColor *color = gROOT->GetColor(1014);  color->SetRGB(0.84, 0.02, 0.05);  return 1014;}
	if (qdc > v[14])                { TColor *color = gROOT->GetColor(1015);  color->SetRGB(0.90, 0.03, 0.00);  return 1015;}    // red side
	else return 0;
}

void assembly()
{
	//--- Definition of a simple geometry
	TGeoManager *geom = new TGeoManager("Assemblies","Geometry using assemblies");

	//--- define some materials
	// (const char* name, Double_t a, Double_t z, Double_t rho, Double_t radlen = 0, Double_t intlen = 0)
	TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0,0,0);
	//TGeoMaterial *matAl = new TGeoMaterial("Al", 26.98,13,2.7);
	//TGeoMaterial *Fe = new TGeoMaterial("Fe",55.845,26,7.87);
	//TGeoMaterial *matSi = new TGeoMaterial("Si", 28.085,14,2.329);
	//TGeoMaterial *Cu = new TGeoMaterial("Cu",63.549,29,8.92);

	//--- define some tracking media
	TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1, matVacuum);
	//TGeoMedium *Al = new TGeoMedium("Aluminium",2, matAl);

	//--- create the starting point of the volume structure - 
	// half length units are in cm
	Double_t xtop=1000;
	Double_t ytop=1000;
	Double_t ztop=1000;
	TGeoVolume *top = geom->MakeBox("TOP", Vacuum, xtop,ytop,ztop);
	geom->SetTopVolume(top);
	geom->SetTopVisible(1);

	// declare panels
	TGeoVolume *panel[32];  

	// bottom veto panels-------------------------------------------------------------
	// make box for each layer and fill with 6 panels
	// panel 0-5  lower bottom
	// panel 6-11 upper bottom
	Double_t xlayer = xtop;
	Double_t ylayer = ytop;
	Double_t zlayer = 50;
	TGeoVolume *layerb1 = geom->MakeBox("layerb1", Vacuum, xlayer,ylayer,zlayer);   
	TGeoVolume *layerb2 = geom->MakeBox("layerb2", Vacuum, xlayer,ylayer,zlayer);   

	// put 6 panels in layer 1
	Double_t xpanel = xlayer;
	Double_t ypanel = ylayer/6.;    
	Double_t zpanel = zlayer;
	Char_t panel_name[50];
	Double_t ypos;

	for (Int_t i=0; i<6;i++){
		sprintf(panel_name,"panel%d",i);
		panel[i]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel,zpanel);
		ypos = -ylayer+(2*i+1)*ypanel;  
		layerb1->AddNode(panel[i],i, new TGeoTranslation(0,ypos,0));

		sprintf(panel_name,"panel%d",i);
		panel[6+i]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel,zpanel);
		ypos = -ylayer+(2*i+1)*ypanel;  
		layerb2->AddNode(panel[6+i],i, new TGeoTranslation(0,ypos,0));
	}
  
	top->AddNode(layerb1,1, new TGeoTranslation(0,0,-ztop+zlayer));  // add bottom layers to mother

	TGeoRotation *rot1 = new TGeoRotation();  
	rot1->RotateZ(90);   
	top->AddNode(layerb2,2, new TGeoCombiTrans(0,0,-ztop+3*zlayer,rot1)); //add second layer

	// top veto panels-------------------------------------------------------------
	// panel 17,18 top outer
	// panel 20,21 top inner   
	TGeoVolume *layert1 = geom->MakeBox("layert1", Vacuum, xlayer,ylayer,zlayer);   
	TGeoVolume *layert2 = geom->MakeBox("layert2", Vacuum, xlayer,ylayer,zlayer);   
	ypanel = ylayer/2.;    

	for (Int_t i=0; i<2;i++){
	    sprintf(panel_name,"panel%d",i);
	    panel[17+i]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel,zpanel);
	    ypos = -ylayer+(2*i+1)*ypanel;  
	    layert1->AddNode(panel[17+i],i, new TGeoTranslation(0,ypos,0));

	    sprintf(panel_name,"panel%d",i);
	    panel[20+i]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel,zpanel);
	    ypos = -ylayer+(2*i+1)*ypanel;  
	    layert2->AddNode(panel[20+i],i, new TGeoTranslation(0,ypos,0));
  	}
  	// add top layers to mother
  	top->AddNode(layert1,1, new TGeoTranslation(0,0,ztop-zlayer));
  	TGeoRotation *rot2 = new TGeoRotation();
  	rot2->RotateZ(-90);  
  	top->AddNode(layert2,2, new TGeoCombiTrans(0,0,ztop-3*zlayer,rot2));

	// north veto panels-------------------------------------------------------------
	// panel 15,16 north outer
	// panel 19,23 north inner
	TGeoVolume *layern1 = geom->MakeBox("layern1", Vacuum, xlayer,ylayer,zlayer);   
	TGeoVolume *layern2 = geom->MakeBox("layern2", Vacuum, xlayer,ylayer,zlayer);   

	xpanel = xlayer - 4*zlayer;    
	ypanel = ylayer/2.;    

	// 15, 16 (outer)
	for (Int_t i=0; i<2;i++){
	    sprintf(panel_name,"panel%d",i);
	    panel[15+i]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel,zpanel);
	    ypos = -ylayer+(2*i+1)*ypanel;  
	    layern1->AddNode(panel[15+i],i, new TGeoTranslation(0,ypos,0));
	 }
	//19
	sprintf(panel_name,"panel%d",19);
	panel[19]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel,zpanel);
	ypos = -ylayer+(2*0+1)*ypanel;  
	layern2->AddNode(panel[19],19, new TGeoTranslation(0,ypos,0));
	//23
	sprintf(panel_name,"panel%d",23);
	panel[23]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel,zpanel);
	ypos = -ylayer+(2*1+1)*ypanel;  
	layern2->AddNode(panel[23],23, new TGeoTranslation(0,ypos,0));

	// add north layers to mother
	TGeoRotation *rot3 = new TGeoRotation();
	rot3->RotateY(90);   
	top->AddNode(layern1,1, new TGeoCombiTrans(xtop-zlayer,0,0,rot3));
	top->AddNode(layern2,2, new TGeoCombiTrans(xtop-3*zlayer,0,0,rot3));

	// west veto panels-------------------------------------------------------------
	// panel 12,13 west inner
	// panel 14,22 west outer
	TGeoVolume *layerw1 = geom->MakeBox("layerw1", Vacuum, xlayer,ylayer,zlayer);   
	TGeoVolume *layerw2 = geom->MakeBox("layerw2", Vacuum, xlayer,ylayer,zlayer);   
	xpanel = xlayer - 4*zlayer;    
	ypanel = (ylayer-2*zlayer)/2.;    

	// 12, 13 (inner)
	for (Int_t i=0; i<2;i++){
		sprintf(panel_name,"panel%d",i);
		panel[12+i]= geom->MakeBox(panel_name, Vacuum, ypanel,xpanel,zpanel);
		ypos = -ylayer+(2*i+1)*ypanel;  
		layerw1->AddNode(panel[12+i],i, new TGeoTranslation(ypos,0,0));
	}
	//14
	sprintf(panel_name,"panel%d",14);
	panel[14]= geom->MakeBox(panel_name, Vacuum, ypanel,xpanel,zpanel);
	ypos = -ylayer+(2*0+1)*ypanel;  
	layerw2->AddNode(panel[14],14, new TGeoTranslation(ypos,0,0));
	//22
	sprintf(panel_name,"panel%d",22);
	panel[22]= geom->MakeBox(panel_name, Vacuum, ypanel,xpanel,zpanel);
	ypos = -ylayer+(2*1+1)*ypanel;  
	layerw2->AddNode(panel[22],22, new TGeoTranslation(ypos,0,0));

	// add west layers to mother
	TGeoRotation *rot4 = new TGeoRotation();
	rot4->RotateX(90);   
	top->AddNode(layerw1,1, new TGeoCombiTrans(0,ytop-zlayer,0,rot4));
	top->AddNode(layerw2,2, new TGeoCombiTrans(0,ytop-3*zlayer,0,rot4));


	///////////////////////// NEW ADDITIONS

	// EAST veto panels-------------------------------------------------------------
	// panel 28,30 EAST inner
	// panel 29,31 EAST outer

	// Z-axis points toward TOP
	// Y-axis points toward EAST
	// X-axis points toward NORTH

	TGeoVolume *layerE1 = geom->MakeBox("layerE1", Vacuum, xlayer,ylayer,zlayer);   
	TGeoVolume *layerE2 = geom->MakeBox("layerE2", Vacuum, xlayer,ylayer,zlayer);   
	xpanel = xlayer - 4*zlayer;    
	ypanel = (ylayer-2*zlayer)/2.; 

	//28
	sprintf(panel_name,"panel%d",28);
	panel[28]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel,zpanel);
	ypos = -ylayer+(2*0+1)*ypanel;  
	layerE1->AddNode(panel[28],28, new TGeoTranslation(0,ypos,0));

	//30
	sprintf(panel_name,"panel%d",30);
	panel[30]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel,zpanel);
	ypos = -ylayer+(2*1+1)*ypanel;  
	layerE1->AddNode(panel[30],30, new TGeoTranslation(0,ypos,0));

	//29
	sprintf(panel_name,"panel%d",29);
	panel[29]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel,zpanel);
	ypos = -ylayer+(2*0+1)*ypanel;  
	layerE2->AddNode(panel[29],29, new TGeoTranslation(0,ypos,0));

	//31
	sprintf(panel_name,"panel%d",31);
	panel[31]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel,zpanel);
	ypos = -ylayer+(2*1+1)*ypanel;  
	layerE2->AddNode(panel[31],31, new TGeoTranslation(0,ypos,0));

	// add EAST layers to mother
	TGeoRotation *rot5 = new TGeoRotation();
	rot5->RotateX(90);
	rot5->RotateY(90);   
	top->AddNode(layerE1,1, new TGeoCombiTrans(0,-ytop+zlayer,0,rot5));
	top->AddNode(layerE2,2, new TGeoCombiTrans(0,-ytop+3*zlayer,0,rot5));


	// SOUTH veto panels-------------------------------------------------------------
	// panel 24,26 SOUTH inner
	// panel 25,27 SOUTH outer

	TGeoVolume *layerS1 = geom->MakeBox("layerS1", Vacuum, xlayer,ylayer,zlayer);   
	TGeoVolume *layerS2 = geom->MakeBox("layerS2", Vacuum, xlayer,ylayer,zlayer);   

	xpanel = xlayer - 4*zlayer;    
	ypanel = (ylayer-4*zlayer)/2.; 

	//24
	sprintf(panel_name,"panel%d",24);
	panel[24]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel,zpanel);
	ypos = -ylayer+(2*0+1)*ypanel;  
	layerS1->AddNode(panel[24],24, new TGeoTranslation(0,ypos,0));

	//26
	sprintf(panel_name,"panel%d",26);
	panel[26]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel,zpanel);
	ypos = -ylayer+(2*1+1)*ypanel;  
	layerS1->AddNode(panel[26],26, new TGeoTranslation(0,ypos,0));

	//25
	sprintf(panel_name,"panel%d",25);
	panel[25]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel,zpanel);
	ypos = -ylayer+(2*0+1)*ypanel;  
	layerS2->AddNode(panel[25],25, new TGeoTranslation(0,ypos,0));

	//27
	sprintf(panel_name,"panel%d",27);
	panel[27]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel,zpanel);
	ypos = -ylayer+(2*1+1)*ypanel;  
	layerS2->AddNode(panel[27],27, new TGeoTranslation(0,ypos,0));

	// add SOUTH layers to mother
	TGeoRotation *rot6 = new TGeoRotation();
	rot6->RotateY(90);   
	top->AddNode(layerS1,1, new TGeoCombiTrans(-xtop+zlayer,+4*zlayer,0,rot6)); // need to move 2x the width of a panel in the -y direction
	top->AddNode(layerS2,2, new TGeoCombiTrans(-xtop+3*zlayer,+4*zlayer,0,rot6));

	// done adding panels! -------------------------------------------------------------

	// close the geometry
	geom->CloseGeometry();
	geom->SetVisLevel(4);
	geom->SetVisOption(0);

	TCanvas *ecan = new TCanvas("ecan","veto hits",0,0,700,700);
	top->Draw(); // first time makes a blank screen

	// add some geometry markers
	// these DON'T FUCKING WORK, need to fix!
	//TPaveText *pt = new TPaveText(-.8,-.8,-.5,-.7,"br");
	//pt->AddText("This is the South Side"); 
	//pt->Draw();


	// Coloring of particular panels
	// help with orientation
	panel[0]->SetLineColor(kBlue);  // lower bottom
	panel[0]->SetLineWidth(3.0); 
	panel[6]->SetLineColor(kGreen);  // upper bottom
	panel[6]->SetLineWidth(3.0); 
	panel[17]->SetLineColor(kRed); // top outer
	panel[17]->SetLineWidth(3.0); 
	panel[20]->SetLineColor(kViolet); // top inner
	panel[20]->SetLineWidth(3.0); 
	panel[15]->SetLineColor(kPink); // north outer
	panel[15]->SetLineWidth(3.0); 
	panel[19]->SetLineColor(kOrange); // north inner
	panel[19]->SetLineWidth(3.0); 
	panel[12]->SetLineColor(kTeal); // west inner
	panel[12]->SetLineWidth(3.0); 
	panel[14]->SetLineColor(kSpring); // west outer
	panel[14]->SetLineWidth(3.0); 
	panel[24]->SetLineColor(kMagenta); // south inner
	panel[24]->SetLineWidth(3.0);   
	panel[28]->SetLineColor(kCyan); // south outer
	panel[28]->SetLineWidth(3.0); 

	cout << "hit enter to draw the first time: ";
	getchar();
	cout << endl;
	top->Draw();

	/*
	// start drawing events -------------------------------------------------------------  
	// loop over hit pattern from input file and re-draw
	Int_t runNumber = 0;
	Int_t entry = 0;
	Int_t eventCount = 0;
	Double_t scalerTime = 0;
	Int_t qdcVals[32]={0};

	ifstream VetoHitsFile;
	VetoHitsFile.open("/home/wisecg/Dev/veto/txtLists/events.txt");
	while(!VetoHitsFile.eof()){

		// wait for keystroke before drawing this event
		getchar();

		// qdc range: ~150 - ~4000
		// corrupted scaler times will equal -8.58994e+07.  
		VetoHitsFile >> runNumber >> entry >> eventCount >> scalerTime;

		printf("runNumber:%i  entry:%i  eventCount:%i  scalerTime:%.5f \n qdcVals:",runNumber,entry,eventCount,scalerTime);
		for (Int_t j=0;j<24;j++){
		  VetoHitsFile >> qdcVals[j];
		  cout << qdcVals[j] << " ";
		}
		cout << endl;
		if (scalerTime<0) cout << "WARNING! Corrupted scaler, no time info available." << endl;

		// list of panels:
		// panel 0-5  lower bottom
		// panel 6-11 upper bottom
		// panel 12,13 west inner   ** did mapping change when card 18 was added? 
		// panel 15,16 north outer  ** why would west and north be on card 11?
		// panel 17,18 top outer    ** were the top panels not installed first?
		// panel 20,21 top inner    **  
		// panel 19-23 north inner
		// panel 14,22 west outer

		// color particular panels based on qdc value
		for (Int_t k = 0; k<24; k++) {
		  panel[k]->SetLineColor(coloring(qdcVals[k],1));
		  if (qdcVals[k] > 0) panel[k]->SetLineWidth(3.0);
		  else panel[k]->SetLineWidth(0.0);
		}

		top->Draw();
		} // end of loop over hits

  // end drawing events -------------------------------------------------------------  
  */
  
}



//--------------------------------------------------------------

void vetoDisplay32()  {
	assembly();
}
