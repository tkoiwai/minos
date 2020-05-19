//----------! Developed by C. Santamaria / CEA Saclay !----------
//----------!      Version date :: 2014/11/24         !----------
// Code to decode the ridf format file for MINOS
// Calibration, Analysis and Tracking
// Version for SEASTAR3 exp. (no vertex reconstruction yet or beam tracking)

// > make
// > ./AnalyzerMINOS RunNumber


#include <fstream>
#include <bitset>
#include <string>
#include <iostream>
#include <cstdlib>
#include <vector>
#include "TSystem.h"
#include "TCanvas.h"
#include "TArtStoreManager.hh"
#include "TArtEventStore.hh"
#include "TArtBigRIPSParameters.hh"
#include "TArtCalibFocalPlane.hh"
#include "TArtDALIParameters.hh"
#include "TArtMINOSParameters.hh"
#include "TArtCalibCoin.hh"
#include "TArtCalibPID.hh"
#include "TArtCalibDALI.hh"
#include "TArtCalibMINOS.hh"
#include "TArtCalibMINOSData.hh"
#include "TArtAnalyzedMINOS.hh"
#include "TArtTrackMINOS.hh"
#include "TArtTrackMINOSData.hh"
#include "TArtVertexMINOS.hh"
#include "TArtCalibPPAC.hh"
#include "TArtCalibPlastic.hh"
#include "TArtFocalPlane.hh"
#include "TArtEventInfo.hh"
#include "TArtPlastic.hh"
#include "TArtPPAC.hh"
#include "TArtRecoPID.hh"
#include "TArtRecoRIPS.hh"
#include "TArtRecoTOF.hh"
#include "TArtRecoBeam.hh"
#include "TArtPPAC.hh"
#include "TArtBeam.hh"
#include "TArtTOF.hh"
#include "TArtRIPS.hh"
#include "TTree.h"
#include "TFile.h"
#include "TClonesArray.h"
#include <vector>
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TGraph.h"
#include <TMinuit.h>
#include <TVirtualFitter.h>
#include <TRegexp.h>
#include <TFitter.h>
#include <time.h>
#include <sys/stat.h>
#include "/home/apollo/root-5.34/include/Math/Vector3D.h"
#include "TVector3.h"
#include "/home/koiwai/analysis/include/liboffline/TMinosClust.h"
#include "/home/koiwai/analysis/include/liboffline/TMinosResult.h"
#include "/home/koiwai/analysis/include/liboffline/Tracking.h"
#include "/home/koiwai/analysis/segidlist.hh"
#include "TArtRawFeminosDataObject.hh"
#include "TCutG.h"
#include "TApplication.h"

#include "/home/koiwai/analysis/include/time_tk.h"

using namespace std;
using namespace ROOT::Math;

inline Bool_t exists_test (const std::string&);
inline Bool_t exists_test (const TString&);

Double_t angleRotation = 30.7*TMath::DegToRad();
Double_t offsetx       = 1.69192;
Double_t offsety       = 0.60047;

Double_t c0 = 0.2;
Double_t c1 = 0.012;

//Definitions to use TMinuit
Tracking     *Tracking_functions;
TClonesArray data_result;
TMinosResult *minosdata_result;

Int_t NclusterFit;
void SumDistance(Int_t &, Double_t *, Double_t & sum, Double_t * par,  Int_t);

// function to exit loop at keyboard interrupt.
Bool_t stoploop = false;
void stop_interrupt()
{
  printf("keyboard interrupt\n");
  stoploop = true;
}

Int_t main(Int_t argc, char** argv)
{

  initiate_timer_tk();
  
  const char *ROOTFILEDIR = "/home/koiwai/analysis/rootfiles/";
  const char *RIDFFILEDIR = "/home/koiwai/analysis/ridf/sdaq02/";
  
  Bool_t IsChTrue(Int_t,Int_t);

  //time_t start,stop;
  //time(&start);

  // Input handling
  Int_t RunNumber; 
  if(argc < 1){
    cerr << "Missing RIDF file argument" << endl;
    exit(EXIT_FAILURE);
  }
  RunNumber = atoi(argv[1]);
  cout << "\n Run Number: " << RunNumber << "\n\n";
  if (RunNumber == 39 || RunNumber == 69 || RunNumber == 86 || RunNumber == 89 ||
      RunNumber == 118|| RunNumber == 139||RunNumber == 171|| RunNumber == 176 ||
      RunNumber == 178)
    exit(1); 
  // RIDF filename definition
  string ridffile = RIDFFILEDIR;
  ridffile.append(Form("run%04d.ridf.gz",RunNumber));
  cout << " RIDF file: " << ridffile << "\n\n";
  if (!exists_test(ridffile)){
    cerr << " ERROR - '" << ridffile << "' does not exist '" << "\n\n";
    exit(EXIT_FAILURE);
  }

  // Defining outfile
  string rootfileout; 
  char  *minosdir="/home/koiwai/analysis/rootfiles/minos/cal_new/";
  if (gSystem->OpenDirectory(minosdir) == 0){
    cerr << " ERROR - '" << minosdir << "' does not exist '" << "\n\n";
    exit(EXIT_FAILURE);
  }
  rootfileout = minosdir;
  rootfileout.append(Form("Tracks_run_%04d.root",RunNumber));
  
  cout << "\n Output ROOT file: " << rootfileout << endl;
  TFile * fout = new TFile(rootfileout.c_str(),"RECREATE");
  TTree * tree = new TTree("tree","MINOS tracks");
  
  //Open ridf
  TArtStoreManager *sman = TArtStoreManager::Instance();
  //TArtCalibCoin  *myInfo = new TArtCalibCoin();
  
  UInt_t   runnumber;
  Long64_t evenumber;
  TArtEventStore *estore = new TArtEventStore();
  estore->SetInterrupt(&stoploop);
  estore->Open(ridffile.c_str());

  // Create MINOSParameters to get ".xml"
  TArtMINOSParameters  *setup = new TArtMINOSParameters("MINOSParameters","MINOSParameters");
  setup->LoadParameters(const_cast<char *>("/home/koiwai/analysis/db/MINOS.xml"));
  TArtCalibMINOS       *CalibMINOS    = new TArtCalibMINOS();
  TArtCalibMINOSData   *minos; 
  TArtAnalyzedMINOS    *AnalyzedMINOS = new TArtAnalyzedMINOS(CalibMINOS);
  TArtTrackMINOS       *TrackMINOS    = new TArtTrackMINOS();
  
  TArtRawEventObject   *fEvent;
  TArtRawSegmentObject *seg;
  TArtRawDataObject    *d;
  
  //Load function library
  Tracking_functions = new Tracking();

  //EventInfo is important for the fBit information to know the trigger!
  TClonesArray * info_array = (TClonesArray *)sman->FindDataContainer("EventInfo");
  
  //Variables for the run
  Int_t    evtOrig;         //the same as neve
  Int_t    nch         = 7; //Number of channels in the coincidences register, i think ... 
  Double_t MINOSthresh = 15;
  Double_t TimeBinElec = 30;     // in ns
  Double_t Tshaping    = 333.3;  // in ns
  Double_t VDrift      = 0.0368; // in mm/ns
  Double_t DelayTrig   = 11330;  // in ns
  Double_t StopT;

  TClonesArray fitdata;
  fitdata.SetClass("TMinosClust");
  data_result.SetClass("TMinosResult"); // A TClonesArray defined before
  Int_t         trackNbr;
  Int_t         trackNbr_FINAL;
  Double_t      Pulser_charge;
  vector<Int_t> EventInfo_fBit;
  // Variables for the fits.. none of them is initialized
  Double_t pStart[4];
  Double_t parFit_temp[4], err_temp[4];
  Double_t chi[2];
  Int_t    fitStatus[2];
  Double_t arglist[10];
  Int_t    iflag;
  Int_t    nvpar,nparx;
  Double_t amin,edm, errdef;
  // Variables only filled when trackNbr_FINAL>=1 // Needed at all?
  vector<Double_t> length, chargeTot;
  vector<Double_t> parFit1, parFit2, parFit3, parFit4;
  
  Double_t maxCharge = 0.;
  Int_t    filled    = 0;
  vector<Double_t> Xpad, Ypad, Qpad, XpadNew, YpadNew, QpadNew, ZpadNew;
  vector<Double_t> XpadTemp, YpadTemp, QpadTemp; 
  vector<Int_t> clusterringboolTemp; 
  vector<Int_t> clusterringbool;
  vector<Int_t> clusternbr;
  vector<Int_t> clusterpads;
  
  Int_t Iteration     = 0;
  Int_t filter_result = 0;
  Int_t indexfill     = 0;
  Bool_t fitbool      = false; 
  Int_t fit2DStatus   = 0;
  Double_t Chi2       = 0.; //not initialized
  Double_t x_mm,y_mm,z_mm,q_pad,t_pad, r_mm;
  
  TF1 *fit_function = new TF1("fit_function",
			      Tracking_functions,&Tracking::conv_fit,
			      0, 511, 3,
			      "Tracking","conv_fit"); // only 3 param. because baseline is fixed. not ini
  
  Double_t hfit_max, hfit_max_T, T_min, T_max;
  
  TH1F *    hfit  = new TH1F("hfit","hfit",512,0,512);
  TGraph* grxztmp = new TGraph(); //not init
  TGraph *gryztmp = new TGraph(); //not init
  //2nd step variables
  Int_t npoint            = 0;
  Int_t cluster_temp      = 0; 
  Int_t ringsum           = 0;
  Int_t ringtouch[18]     = {0};
  Double_t zmax           = 0.;
  Int_t allevt_2pfiltered = 0;
  Int_t array_final       = 0;
  
  vector<TGraph >  gryz;
  vector<TGraph >  grxz;
  vector<Double_t> xin, yin, zin, qin, xout, yout, zout, xoutprime, youtprime, zoutprime, qout; //not ini
  
  TMinuit *min;   //not ini
  TVector3 point; //not ini
  
  //Branches definition
  tree->Branch("fitdata",         &fitdata);
  tree->Branch("EventNumber",     &evtOrig,       "EventNumber/I");
  tree->Branch("Ridf_EventNumber",&evenumber,     "Ridf_EventNumber/L");
  tree->Branch("RunNumber",       &runnumber,     "RunNumber/I");
  tree->Branch("NumberTracks",    &trackNbr_FINAL,"NumberTracks/I");
  tree->Branch("VDrift",          &VDrift,        "VDrift/D");
  tree->Branch("DelayTrig",       &DelayTrig,     "DelayTrig/D");
  
  tree->Branch("chargeTot","vector<Double_t>",&chargeTot);  
  tree->Branch("length",   "vector<Double_t>",&length);
  tree->Branch("parFit1",  "vector<Double_t>",&parFit1);
  tree->Branch("parFit2",  "vector<Double_t>",&parFit2);
  tree->Branch("parFit3",  "vector<Double_t>",&parFit3);
  tree->Branch("parFit4",  "vector<Double_t>",&parFit4);

  tree->Branch("Pulser_charge",&Pulser_charge,"Pulser_charge/D");

  tree->Branch("MyEventInfo_fBit",&EventInfo_fBit);

  // Get parameters for the MINOS ANALYSIS
  ifstream ConfigFile;
  ConfigFile.open("/home/koiwai/analysis/db/config_MINOSdrift.txt");
  Int_t HeaderFile;

  while(ConfigFile.is_open()) {
    ConfigFile >> HeaderFile >> DelayTrig >> StopT >> VDrift;
      cout << "HEADER: " << HeaderFile << endl;    
    if(HeaderFile == RunNumber)
      break;
  }
  ConfigFile.close();
  
  cout << endl;
  cout << " *** MINOS Configuration Parameters *** " << endl;
  cout << "MINOSthresh = " << MINOSthresh << " (bins)"  << endl ;
  cout << "TimeBinElec = " << TimeBinElec << " (ns)"    << endl ;    
  cout << "   Tshaping = " << Tshaping    << " (ns)"    << endl;  
  cout << "  DelayTrig = " << DelayTrig   << " (ns)"    << endl ;   
  cout << "     VDrift = " << VDrift      << " (mm/ns)" << endl;
  cout << endl;
  cout<< "-----------------------------------------------"<<endl;
  cout << "Conversion started (Wait a few seconds for display) " << endl;
  cout<< "-----------------------------------------------"<<endl;

  prepare_timer_tk();
  
  Int_t neve = 0;

  cout << argc << endl;

  int nEntry;
  if(argc == 3)
    nEntry = atoi(argv[2]);
  else if(argc == 2)
    nEntry = 1000000;
  cout << "You will process " << nEntry << " events.\n"  << endl;

  while(estore->GetNextEvent()&&neve<nEntry){
  //while(estore->GetNextEvent() ) {

    start_timer_tk(neve,nEntry,10);
    
    //if(neve%1000==0)
    //  cout << "Event " << neve << "\r" << flush;

    
    evtOrig = neve;

    //cout << "ok1" << endl;


    //myInfo->LoadData(); //coin register
    runnumber = ((TArtEventInfo *)info_array->At(0))->GetRunNumber();
    evenumber = ((TArtEventInfo *)info_array->At(0))->GetEventNumber();

    //cout << "ok2"  << endl;
      
    //get the trigger distribution
    fEvent = (TArtRawEventObject *)sman->FindDataContainer("RawEvent"); //TArtRawEventObject
    EventInfo_fBit.clear();  //vector

    for(Int_t i=0; i<fEvent->GetNumSeg(); i++){ // loop over the raw event segment
      seg = fEvent->GetSegment(i);              // TArtRawSegmentObject
      Int_t detector = seg->GetDetector();
      if(COIN == detector){                     // if the data comes from the coincidences register
	d = seg->GetData(0);                    // TArtRawDataObject
	Int_t val = d->GetVal();
	if (seg->GetFP()==F13){                 // fp is 13 or 63
	  for (Int_t id=1; id<nch+1; ++id){     // loop over 1 to 7
	    if(IsChTrue(id,val)) 
	      EventInfo_fBit.push_back(id-1);
	  }
	}
      }
    }

    //Clear & Reset variables
    fitdata.Clear();
    data_result.Clear();
    trackNbr       = 0;
    trackNbr_FINAL = 0;
    Pulser_charge  = 0.;

    length.clear();
    chargeTot.clear();
    parFit1.clear();
    parFit2.clear();
    parFit3.clear();
    parFit4.clear();
       
    maxCharge = 0.;
    filled    = 0;
    Xpad.clear();
    Ypad.clear();
    Qpad.clear();
    XpadNew.clear();
    YpadNew.clear();
    ZpadNew.clear();
    QpadNew.clear();
    clusterringbool.clear();
    clusternbr.clear();
    clusterpads.clear();
    Iteration     = 0;
    filter_result = 0;
    indexfill     = 0;
    fitbool       = false;
    fit2DStatus   = 0;
    x_mm = 0.; y_mm = 0.; z_mm = 0.; q_pad = 0.; t_pad = 0.;
    hfit->Reset();

    npoint      = 0;
    ringsum     = 0;
    zmax        = 0.;
    array_final = 0;
    grxz.clear();
    gryz.clear();

    // Making MINOS Reconstruction
    CalibMINOS->ClearData();
    CalibMINOS->ReconstructData();
    // Get the pulser charge
    for(Int_t i=0; i<CalibMINOS->GetNumCalibMINOS(); i++){ // Loop over pads
      minos = CalibMINOS->GetCalibMINOS(i); 
      Double_t maxCharge = 0.;
      if(minos->GetDetID() == 3 && minos->GetAsic()==0 && minos->GetChannel()==64){ // There is never detID==3
	for(Int_t j=0; j<minos->GetNData(); j++){ // Loop over samples
	  if(minos->GetCalibValue(j)>maxCharge) maxCharge = minos->GetCalibValue(j);
	}
	Pulser_charge = maxCharge;
	cout << "MAX_charge 1 :" << maxCharge << endl;
      }
    }

    //OFFLINE tracking algorithm 
    //----------------------------------------------------------------------
    TMinosClust *minosfitdata;
    //----------------------------------------------------------------------------------------------
    //  IN EVENT LOOP  : STEP1
    //  Filling vectors with (x,y,q) information for TPC and count the number of filled (valid pads)
    //-----------------------------------------------------------------------------------------------
    for(Int_t i=0; i<CalibMINOS->GetNumCalibMINOS(); i++){ // Loop over pads
      minos     = CalibMINOS->GetCalibMINOS(i);            // Get pad
      maxCharge = 0.;
      x_mm      = minos->GetX()+offsetx;                   // Get positions
      y_mm      = minos->GetY()+offsety;
      if(minos->GetDetID() != 0)                           //2 and 3 are dssd or something else
	continue; 
      else {
	if(!(abs(x_mm)<0.0001 && abs(y_mm)<0.0001) ){  // NON connected MINOS (TPC) channels...
	  for(Int_t j=0; j<minos->GetNData(); j++){    // Loop over samples
	    if(minos->GetCalibValue(j)>maxCharge){
	      maxCharge = minos->GetCalibValue(j);     //get max charge of each pad
	      Pulser_charge = maxCharge;
	      //cout << "MAX_charge 2 :" << maxCharge << endl;
	    }
	  }
	  if(maxCharge>=MINOSthresh){ //if above threshold
	    Xpad.push_back(x_mm);
	    Ypad.push_back(y_mm);
	    Qpad.push_back(maxCharge);
	    filled++; //number of valid pads in the event
	  }
	}
      }
    }
    //------------------------------------------------------------------
    //  IN EVENT LOOP  : STEP 2
    //  Modified (xy) Hough transform (only for TPC)
    //------------------------------------------------------------------
    if(filled>0){                             // If there are good pads
      while(Xpad.size()>=10 && Iteration<20){ // Take events until there are less than 10 pads or got 20 clusters
	filter_result = 0;
	Iteration++;
	XpadTemp.clear();                
	YpadTemp.clear();                
	QpadTemp.clear();                
	clusterringboolTemp.clear();

	filter_result = Tracking_functions->Hough_modified(&Xpad, &Ypad, &Qpad, 
							   &XpadTemp, &YpadTemp, &QpadTemp, 
							   &clusterringboolTemp);

	if(filter_result<0) break; 
	if(filter_result>10 && clusterringboolTemp.back()==1){ // More than 10 pads and more than 2 in the first ring
	  trackNbr++; //proper track!
	  for(Int_t ik=0; ik<filter_result; ik++){ // Loop over the number of good pads
	    XpadNew.push_back(XpadTemp[ik]); // Fill the vectors containing the pads that form the cluster  
	    YpadNew.push_back(YpadTemp[ik]);                
	    ZpadNew.push_back(-10000);       // No z pad here..
	    QpadNew.push_back(QpadTemp[ik]);  
	    clusterringbool.push_back(clusterringboolTemp[ik]);                
	    clusternbr.push_back(Iteration); //cluster number
	    clusterpads.push_back(filter_result);
	  }
	}
      }
    }
    for(UInt_t il=0; il<XpadNew.size(); il++){ // Loop over the number of new pads
      minosfitdata = (TMinosClust*)fitdata.ConstructedAt(il);
      minosfitdata->Set(XpadNew[il], YpadNew[il],
			-10000, -10000, QpadNew[il], 
			clusternbr[il],clusterpads[il], 0.);
      ZpadNew.push_back(-10000.);
    }
	
    //------------------------------------------------------------------
    //  IN EVENT LOOP  : STEP 3
    //  Following analysis only if #tracks = 1||2||3||4
    //------------------------------------------------------------------
    if(trackNbr>0 && trackNbr<5) {
	if(filled==0) cerr << "Error !!!" << endl;
	//-------------------------------------------------------
	//  STEP 3.1: Fitting taken pads for Qmax and Ttrig info
	//-------------------------------------------------------
	for(Int_t i=0; i<CalibMINOS->GetNumCalibMINOS(); i++) { // Loop over all pads
	  minos = CalibMINOS->GetCalibMINOS(i);
	  hfit->Reset();
	  fitbool = false;
	  if(minos->GetDetID() != 0) continue; 
	  x_mm = minos->GetX()+offsetx;
	  y_mm = minos->GetY()+offsety;
	  r_mm = sqrt(x_mm*x_mm+y_mm*y_mm);
			       
	  for(UInt_t jj=0; jj<XpadNew.size(); jj++){ // Get the pads that were good in the previous step...
	    if( abs(XpadNew[jj]-x_mm)<0.0001 && abs(YpadNew[jj]-y_mm)<0.0001){
	      fitbool   = true;
	      indexfill = jj;
	      break;
	    }
	  }
	  // Check if New channel is of interest
	  //(if so, we read the Q(t) signal, and we should 
	  // fill the vectors w/ t & q information after fitting E(T))
	  if( fitbool==true ) {
	    for(Int_t j=0; j<minos->GetNData(); j++){ // Loop over samples
	      if(minos->GetCalibValue(j)>=0)
		hfit->SetBinContent(hfit->FindBin(minos->GetCalibTime(j)), minos->GetCalibValue(j)+250);
	    }
	    // Fitting the hfit histogram of last channel if not empty
	    if(hfit->GetSumOfWeights()>0){
	      hfit->GetXaxis()->SetRange(0,510);
	      hfit_max = hfit->GetMaximum();
	      hfit_max_T = hfit->GetMaximumBin();
	      T_min=-1;
	      T_max=-1;
	      // Find the T_min & T_max limits of the signal non zero
	      for(Int_t h=hfit_max_T;h>0;h--) {
		if(T_min == -1 && (hfit->GetBinContent(h))<=250 ){
		  T_min = h;
		  break;
		}
	      }
	      for(Int_t h=hfit_max_T;h<510;h++){
		if(T_max == -1 && (hfit->GetBinContent(h))==0 ){
		  T_max = h;
		  break;
		}
	      }
	      // Take only 1.5*Tshaping before the max if other signals before...
	      if((hfit_max_T-3.5*(Tshaping/TimeBinElec)) > T_min)
		T_min = hfit_max_T-2*Tshaping/TimeBinElec;
	      if((hfit_max_T+10) < T_max || T_max==-1)
		T_max = hfit_max_T+10.;
	      T_min = max(T_min,0.);
	      if(T_max>510) T_max = 510;
	      // Set fit parameters
	      fit_function->SetParameter(0, hfit_max-250.);
	      fit_function->SetParameter(1,hfit_max_T - Tshaping/TimeBinElec);
	      fit_function->SetParameter(2, Tshaping/TimeBinElec);
	      fit_function->SetParLimits(0,0,100000);
	      fit_function->SetParLimits(1,-20,512);
	      fit_function->SetParLimits(2,0,512);
	      // Fit of the signal within the range defined: [T_min, T_max]
	      fit2DStatus = hfit->Fit(fit_function,"QN","",T_min,T_max);
	      Double_t fit_function_max = 0., fit_function_Tpad = 0.;
	      if(fit2DStatus==0) {
		Chi2 = fit_function->GetChisquare();               // Get chi2 of the fit
		fit_function_max = fit_function->GetMaximum();     // Q max
		fit_function_Tpad = fit_function->GetParameter(1); // Tpad
	      }
			
	      // Attribute q_pad and z_mm value (z is inside the tpc)
	      if(fit2DStatus!=0                      ||
		 fit_function_max<=20.               || fit_function_max>=100000. ||
		 fit_function_Tpad<=0.15             || fit_function_Tpad>=513.   ||
		 fit_function->GetParameter(2)<=0.15 || fit_function->GetParameter(2)>=513.){
		  q_pad = hfit_max-250.;
		  z_mm = -10000;
	      }  else {
		// Add to the variables the fit parameters
		t_pad = fit_function_Tpad;//trigger time
			       
		Double_t t1 = t_pad*TimeBinElec;
		Double_t t0 = DelayTrig;
		z_mm = VDrift*(t1-t0); //time bin*(ns/us)*vdrift(mm/ns) =>mm

		Double_t phi    = atan2(y_mm,x_mm);
		Double_t zprime = (z_mm/300)-1.0;
		Double_t yprime = x_mm*cos(angleRotation)-y_mm*sin(angleRotation);

		Double_t Er         = 81-r_mm;
		Double_t IntegralEz = sqrt(1-zprime*zprime);
			       
		Double_t deltaR = c0*Er*IntegralEz*(1+c1/c0*sin(phi+angleRotation));

		Double_t kz=0.06;
		Double_t k=0.2*(1.+kz*(x_mm/r_mm*cos(angleRotation)-y_mm/r_mm*sin(angleRotation)));
		Double_t r_new = (r_mm-81.95)*(1.-k*sqrt(1-pow(z_mm/300.-1,2)) )+81.95;
		//		Double_t r_new = r_mm+(81.95-r_mm)*0.2;
			       			       
		Double_t ratio = r_new/r_mm;
		Double_t x_new = x_mm*ratio;
		Double_t y_new = y_mm*ratio;
			       
		XpadNew[indexfill]=x_new;
		YpadNew[indexfill]=y_new;
			       
		q_pad = fit_function_max-250.;  // Max charge/event
			       
	      }
	      ZpadNew[indexfill] = z_mm; 

	      minosfitdata = (TMinosClust*)fitdata.ConstructedAt(indexfill);
	      minosfitdata->Set(XpadNew[indexfill], YpadNew[indexfill], t_pad*TimeBinElec, 
				ZpadNew[indexfill], q_pad, clusternbr[indexfill], clusterpads[indexfill], Chi2);
	      //Fill the z and q information for next steps (3D Hough filter & 3D fit weighted by charge)
	      QpadNew[indexfill] = q_pad;
	    }//end if histogram not empty
	  }// end if fitbool==true
	  else
	    continue;
	}//END of entries in tclonesarray for the event
	       
	//-------------------------------------------------------
	//  STEP 3.2:  Filtering the tracks off possible noise 
	//             with Hough3D (3*2D planes)
	//-------------------------------------------------------
	for(UInt_t i=0;i<(XpadNew.size());i++) {
	  if(xin.size()>0 && ( (cluster_temp!=Int_t(clusternbr[i]) && i!=0) || i==(XpadNew.size() - 1) ) ) {
	    Tracking_functions->Hough_3D(&xin, &yin, &zin, &qin, &xout, &yout, &zout, &qout);
	    // xout is the vector of valid pads after teh second Hough transformation
	    for(UInt_t ij=0; ij<xout.size();ij++) {
	      if(zout[ij]>zmax) zmax = zout[ij];
	      ringtouch[Int_t((sqrt(xout[ij]*xout[ij]+yout[ij]*yout[ij])-45.2)/2.1)]++; // Add the pad to the ring
	    }
	    for(Int_t ko=0; ko<18; ko++) {
	      if(ringtouch[ko]>0) ringsum++; // Count how many rings are hit
	    }
	    if(zmax>290) ringsum = 16; // Dunnoo, me neither
	    // Tracks of interest: >10 pads and >=12 rings hit
	    if(xout.size()>10 && ringsum>=12) {
	      npoint=0;
	      trackNbr_FINAL++; //number of tracks Double_t filtered
	      Double_t charge_temp = 0.;
	      Double_t length_temp = 0;
	      Double_t zintarget   = 300;
	      Double_t zouttarget  = 0;
	      Int_t indexin=0, indexout=0;
	      for(UInt_t ij=0; ij<xout.size(); ij++){ //loop over good pads
		point.SetXYZ(xout[ij],yout[ij],zout[ij]);  //vector to the pad
		point.RotateZ(angleRotation);
			       
		xoutprime.push_back(point.X()); //rotated points
		youtprime.push_back(point.Y());
		zoutprime.push_back(point.Z());

		grxztmp->SetPoint(npoint,zoutprime[ij],xoutprime[ij]); //create graphs of the xz and yz planes
		gryztmp->SetPoint(npoint,zoutprime[ij],youtprime[ij]);
		
		charge_temp += qout[ij]; //total charge
		
		minosdata_result = (TMinosResult*)data_result.ConstructedAt(array_final);
		minosdata_result->Set(xoutprime[ij], youtprime[ij], zoutprime[ij], qout[ij], trackNbr_FINAL, xout.size(), zmax);
		array_final++; 
		npoint++;
		if(zoutprime[ij]<zintarget) {
		  zintarget=zoutprime[ij];
		  indexin=ij;
		}
		if(zoutprime[ij]>zouttarget){
		  zouttarget=zoutprime[ij];
		  indexout=ij;
		}
	      }
	      //distance inside the tpc
	      if(xout.size()>0)
		length_temp=sqrt(pow(zoutprime[indexin]-zoutprime[indexout],2)
				 +pow(youtprime[indexin]-youtprime[indexout],2)
				 +pow(xoutprime[indexin]-xoutprime[indexout],2));

	      grxz.push_back(*grxztmp);
	      gryz.push_back(*gryztmp);
	      grxztmp->Set(0);
	      gryztmp->Set(0);
	      chargeTot.push_back(charge_temp);
	      length.push_back(length_temp);
	    }

	    xin.clear();
	    yin.clear();
	    zin.clear();
	    qin.clear();
	    xout.clear();
	    yout.clear();
	    zout.clear();
	    xoutprime.clear();
	    youtprime.clear();
	    zoutprime.clear();
	    qout.clear();

	    ringsum=0;// angle between 1st track and z axis in 3D in degrees
	    zmax=0.;
	    for(Int_t ko=0; ko<18; ko++) ringtouch[ko] = 0;

	  }

	  cluster_temp = clusternbr[i];

	  if(!(clusterpads[i]>=10 && clusterringbool[i]==1 && ZpadNew[i]>-10000 && ZpadNew[i]<=320)) continue;
	  else {
	    xin.push_back(XpadNew[i]);
	    yin.push_back(YpadNew[i]);
	    zin.push_back(ZpadNew[i]);
	    qin.push_back(QpadNew[i]);
	  }

	}//end of loop on pads

	//-------------------------------------------------------
	//  STEP 3.3:  Fitting the filtered tracks in 3D 
	//             (weight by charge, TMinuit)
	//-------------------------------------------------------
	//For 1 track found or more (less than 5 in total)
	if(trackNbr_FINAL>=1) {
	  //////////Minimization in 3D to reconstruct track lines
	  allevt_2pfiltered++;
	  for(Int_t itr=0; itr<trackNbr_FINAL; itr++){
	    pStart[0]=0; pStart[2]=0; pStart[1]=1; pStart[3]=3;
	    min = new TMinuit(4);
	    min->SetPrintLevel(-1);
	    arglist[0] = 3;
	    //Fit the 2D lines in xz and yz and get parameters
	    //pStart[0] and [1] are cut and slope of the xz line
	    //pStart[2] and [3] are cut and slope of the yz line
	    Tracking_functions->FindStart(pStart,chi,fitStatus, &grxz.at(itr),&gryz.at(itr));
	    NclusterFit = itr+1;
	    min->SetFCN(SumDistance);
	    // Set starting values and step sizes for parameters
	    min->mnparm(0,"x0",pStart[0],0.1,-500,500,iflag);
	    min->mnparm(1,"Ax",pStart[1],0.1,-10,10,iflag);
	    min->mnparm(2,"y0",pStart[2],0.1,-500,500,iflag);
	    min->mnparm(3,"Ay",pStart[3],0.1,-10,10,iflag);
	    arglist[0] = 100; // number of function calls
	    arglist[1] = 0.000001; // tolerance
	    min->mnexcm("MIGRAD",arglist,2,iflag); // minimization with MIGRAD
		       
	    min->mnstat(amin,edm,errdef,nvpar,nparx,iflag);  //returns current status of the minimization
	    // get fit parameters
	    for(Int_t i = 0; i <4; i++) min->GetParameter(i,parFit_temp[i],err_temp[i]);
		       
	    parFit1.push_back(parFit_temp[0]);
	    parFit2.push_back(parFit_temp[1]);
	    parFit3.push_back(parFit_temp[2]);
	    parFit4.push_back(parFit_temp[3]);
	    delete min;
	  }
		   
		   
	}// end if trackNbr_FINAL>=1
	       
	// tree->Fill();
	xin.clear();
	yin.clear();
	zin.clear();
	qin.clear();
	xout.clear();
	yout.clear();
	zout.clear();
	xoutprime.clear();
	youtprime.clear();
	zoutprime.clear();
	qout.clear();
	//Hough_canvas.clear();

      }//loop for E(T) fits for less than 5 track evts
    tree->Fill();  // fill the tree in ALL cases to obtain the same number of evts for DALI2 analysis
    estore->ClearData();
    neve ++;
  }//loop on events

  //     tree->Print();
  cout<<"Write..."<<endl;
  fout->Write();
  cout<<"Close..."<<endl;
  fout->Close();
  cout<<"Conversion to Root done!"<<endl;

  //time(&stop);
  //printf("Elapsed time: %.1f seconds\n",difftime(stop,start));

  stop_timer_tk(nEntry);
  cout << "AnalyzerMINOS: run " << RunNumber << " finished!: tkok" << endl;
  
  return 0;

}
   
/// Functions to be minimized
void SumDistance(Int_t &, Double_t *, Double_t & sum, Double_t * par,  Int_t)
{
  Int_t nused=0;
  Double_t qtot=0;
  sum = 0;
  for(Int_t i=0; i<data_result.GetEntriesFast(); i++) {
    minosdata_result = (TMinosResult*)data_result.At(i);
    if(minosdata_result->n_Cluster==NclusterFit){
      float x=minosdata_result->x_mm;
      float y=minosdata_result->y_mm;
      float z=minosdata_result->z_mm;
      float q=minosdata_result->Chargemax;
      //if(nused<2)cout<<minosdata_result->n_Cluster<<" "<<x<<" "<<y<<" "<<z<<" "<<q<<endl;
      Double_t d = Tracking_functions->distance2(x, y, z, par);
      sum += d*q;
      nused++;
      qtot+=q;
    }
  }
  //sum/=nused;
  sum/=qtot;
  return;
}

Bool_t IsChTrue(Int_t id, Int_t val)
{
  //     cout<<"hallo "<<std::bitset<8>(val)<<endl;
  return ((1 << (id - 1)) & val) ? true : false;
}

inline Bool_t exists_test (const std::string& name) {
  return ( access( name.c_str(), F_OK ) != -1 );
}
inline Bool_t exists_test (const TString& name) {
  return ( access( name.Data(), F_OK ) != -1 );
}
