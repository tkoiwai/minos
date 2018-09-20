//----------! Developed by C. Santamaria / CEA Saclay !----------
//----------!      Version date :: 2014/11/24         !----------
// Code to decode the ridf format file for MINOS
// Calibration, Analysis and Tracking
// Version for SEASTAR3 exp. (no vertex reconstruction yet or beam tracking)

// > make
// > ./AnalyzerMINOS ../ridf/FileName.ridf   0      -> Creates a root file: ../rootfiles/FileName_minos.root 0:full 1:short
// > ./AnalyzerMINOS ../ridf/FileName.ridf   1      -> Uses the online MINOS code
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
//#include "/home/liliana/Packages/root_5.34.36/math/genvector/inc/Math/Vector3D.h"
#include "TVector3.h"
#include "/home/koiwai/analysis/include/liboffline/TMinosClust.h"
#include "/home/koiwai/analysis/include/liboffline/TMinosResult.h"
#include "/home/koiwai/analysis/include/liboffline/Tracking.h"
#include "/home/koiwai/analysis/segidlist.hh"
#include "TArtRawFeminosDataObject.hh"
#include "TCutG.h"
#include "TApplication.h"
using namespace std;
using namespace ROOT::Math;

bool plots=false;
double angleRotation=30.7*TMath::DegToRad(); //### should be confirmed by myself
double offsetx=1.69192;
double offsety=0.60047;

double c0=0.2;
double c1=0.012;


char ROOTFILEDIR[] = "/home/koiwai/analysis/rootfiles/";

//===== Definitions to use TMinuit =====
Tracking *Tracking_functions;
TClonesArray data_result;
TMinosResult *minosdata_result;
int NclusterFit;
void SumDistance(int &, double *, double & sum, double * par,  int);
bool MINOSOnline = false;

//===== function to exit loop at keyboard interrupt. =====
bool stoploop = false;
void stop_interrupt()
{
  printf("keyboard interrupt\n");
  stoploop = true;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int main(int argc, char** argv)
{
  bool IsChTrue(int,int);
  //===== Variables to calculate the elapsed time of the process =====
  time_t start,stop;
  time(&start);

  TApplication *a = new TApplication("a",0,0);
    // Get ridf  file names
  //------------------------------
  /*
  int RidfRunNumber;
  char* ridffile;
  ridffile = argv[1];
  RidfRunNumber = atoi(ridffile);
  
  int online= atof(argv[2]);
  //  angleRotation= atoi(argv[2]);
  if(online == 1) MINOSOnline = true;
  TString RidfRunPath = Form("../../ridf/run%04d.ridf.gz",RidfRunNumber);
  cout << " *** RIDF file: " << RidfRunPath << endl;
  */
  Int_t FileNum = TString(argv[1]).Atoi();
  TString ridffile = Form("/home/koiwai/analysis/ridf/sdaq02/run%04d.ridf.gz",FileNum);

  
  /*
  //PID Cuts
  //----------------------
  //Open merged file to get PID
  TFile *ffrag=new TFile(Form("../../rootfiles/run%d/run%d_MergedAll.root",RidfRunNumber,RidfRunNumber));
  TTree *tfrag=(TTree*) ffrag->Get("tree");
  double aoqSamurai_cut, zetSamurai_cut;
  double aoqSamurai, zetSamurai;
  double aoq37_313, zet313;
  double aoqBR_cut, zetBR_cut;
  tfrag->SetBranchAddress("aoqSamurai", &aoqSamurai) ;
  tfrag->SetBranchAddress("zetSamurai", &zetSamurai) ;
  tfrag->SetBranchAddress("aoq37_313", &aoq37_313) ;
  tfrag->SetBranchAddress("zet313", &zet313) ;
  */
  TString infnameb = Form("/home/koiwai/analysis/rootfiles/ana/beam/ana_beam%04d.root",FileNum);
  TFile *infileb = TFile::Open(infnameb);
  TTree *trb;
  infileb->GetObject("anatrB",trb);
  double aoqBR, zetBR;
  trb->SetBranchAddress("aoqBR",&aoqBR);
  trb->SetBranchAddress("zetBR",&zetBR);
  TString infnames = Form("/home/koiwai/analysis/rootfiles/ana/smri/ana_smri%04d.root",FileNum);
  TFile *infiles = TFile::Open(infnames);
  TTree *trs;
  infiles->GetObject("anatrS",trs);
  double aoqSA, zetSA;
  trs->SetBranchAddress("aoqSA",&aoqSA);
  trs->SetBranchAddress("zetSA",&zetSA);

  trb->AddFriend(trs);

  TFile *brcuts  = new TFile("/home/koiwai/analysis/cutfiles/BRpid.root","");
  TCutG *cbr56ca = (TCutG*)brcuts->Get("br56ca");
  TCutG *cbr56sc = (TCutG*)brcuts->Get("br56sc");
  TCutG *cbr54ca = (TCutG*)brcuts->Get("br54ca");
  TFile *sacuts  = new TFile("/home/koiwai/analysis/cutfiles/SApid.root","");
  TCutG *csa55ca = (TCutG*)sacuts->Get("sa55ca");
  TCutG *csa53ca = (TCutG*)sacuts->Get("sa53ca");
  TCutG *csa55k  = (TCutG*)sacuts->Get("sa55k");
  //-------------
  /*
  TFile *BRcut = new TFile("/home/liliana/Documents/SEASTAR3_Analysis/liliana/Cuts/BRCut_63V.root","READ");
  TCutG *brcut;
  BRcut->GetObject("mycut",brcut);
  BRcut->Close();
  
  TFile *SAMcut1 = new TFile("/home/liliana/Documents/SEASTAR3_Analysis/liliana/Cuts/SamuraiCut_62Ti.root","READ");
  TCutG *samcut1;
  SAMcut1->GetObject("mycut",samcut1);
  SAMcut1->Close();
  TFile *SAMcut2 = new TFile("/home/liliana/Documents/SEASTAR3_Analysis/liliana/Cuts/SamuraiCut_50Ar.root","READ");
  TCutG *samcut2;
  SAMcut2->GetObject("mycut",samcut2);
  SAMcut2->Close();
  TFile *SAMcut3 = new TFile("/home/liliana/Documents/SEASTAR3_Analysis/liliana/Cuts/SamuraiCut_60Ti.root","READ");
  TCutG *samcut3;
  SAMcut3->GetObject("mycut",samcut3);
  SAMcut3->Close();
  TFile *SAMcut4 = new TFile("/home/liliana/Documents/SEASTAR3_Analysis/liliana/Cuts/SamuraiCut_52Ca.root","READ");
  TCutG *samcut4;
  SAMcut4->GetObject("mycut",samcut4);
  SAMcut4->Close();
  TFile *SAMcut5 = new TFile("/home/liliana/Documents/SEASTAR3_Analysis/liliana/Cuts/SamuraiCut_53Ca.root","READ");
  TCutG *samcut5;
  SAMcut5->GetObject("mycut",samcut5);
  SAMcut5->Close();
   TFile *SAMcut6 = new TFile("/home/liliana/Documents/SEASTAR3_Analysis/liliana/Cuts/SamuraiCut_57Sc.root","READ");
  TCutG *samcut6;
  SAMcut6->GetObject("mycut",samcut6);
  SAMcut6->Close();
  */
  //Get root file name
  //------------------------
  /*
  TString rootfilesubpath = Form("../../rootfiles/run%d/",RidfRunNumber);
  TString OutputFilePath;
  if(MINOSOnline) OutputFilePath = rootfilesubpath + Form("run%d_MINOS_Online.root",RidfRunNumber);
  else OutputFilePath = rootfilesubpath + Form("run%d_MINOS_test6.root",RidfRunNumber);
  struct stat infosubpath;
  cout << " *** ROOT file: " << OutputFilePath  << endl;
  if(stat( rootfilesubpath, &infosubpath ) != 0)
    {
      cout <<  "\e[33m  ANALYSIS-Info : Create directory '" << rootfilesubpath <<  "'\e[0m" << endl;
      int goodfile=system("mkdir "+rootfilesubpath);
      cout<<"the value returned by the system was "<<goodfile<<endl;
     }
  
  TFile *fout = new TFile(OutputFilePath.Data(),"RECREATE");
  TTree * tree = new TTree("tree","ridf tree");
  */
  TString rootfile = Form("/home/koiwai/analysis/rootfiles/minos/cal/cal_minos%04d.root",FileNum);
  //TString rootfile = Form("/home/koiwai/analysis/rootfiles/minos/test/cal_minos%04d.root",FileNum);
  TFile *fout = new TFile(rootfile,"RECREATE");
  TTree *tree = new TTree("caltrM","caltrM");
  
  //===== Open ridf =====
  //------------
  TArtStoreManager *sman = TArtStoreManager::Instance();
  //TArtCalibCoin *myInfo = new TArtCalibCoin();
  UInt_t runnumber;
  Long64_t evenumber ;
  TArtEventStore *estore = new TArtEventStore();
  estore->SetInterrupt(&stoploop);
  estore->Open(ridffile);

  //===== Create MINOSParameters to get ".xml" =====
  //------------------------------------
  TArtMINOSParameters *setup = new TArtMINOSParameters("MINOSParameters","MINOSParameters");
  setup->LoadParameters(const_cast<char *>("/home/koiwai/analysis/db/MINOS.xml"));
  //setup->PrintListOfMINOSPara();
  TArtCalibMINOS *CalibMINOS = new TArtCalibMINOS();
  TArtCalibMINOSData *minos ; 
  TArtAnalyzedMINOS *AnalyzedMINOS = new TArtAnalyzedMINOS(CalibMINOS);
  TArtTrackMINOS *TrackMINOS = new TArtTrackMINOS();
  TArtTrackMINOSData *trackdata = NULL;
  TArtTrackMINOSData *trackdata2 = NULL;
  //TArtVertexMINOS *VertexMINOS = new TArtVertexMINOS();
  //TArtRawEventObject *fEvent;// = (TArtRawEventObject *)sman->FindDataContainer("RawEvent");
  TArtRawEventObject *fEvent;
  TArtRawSegmentObject *seg;
  TArtRawDataObject *d;
  
  //===== Load function library =====
  //----------------------
  Tracking_functions = new Tracking();

  //EventInfo is important for the fBit information to know the trigger!
  TClonesArray * info_array = (TClonesArray *)sman->FindDataContainer("EventInfo");
  
  //===== Variables for the run =====
  //------------------------
  int evtOrig; //the same as neve
  int nch = 7; //Number of channels in the coincidences register, i think ... 
  double MINOSthresh=25;
  double TimeBinElec=30; //in ns
  double Tshaping=333.3; // in ns
  double VDrift; //in mm/ns
  double DelayTrig; // in ns
  double StopT;
  TClonesArray fitdata;
  fitdata.SetClass("TMinosClust");
  data_result.SetClass("TMinosResult"); //a TClonesArray defined before
  int trackNbr;
  int trackNbr_FINAL;
  double  Pulser_charge;
  vector<Int_t> EventInfo_fBit;
  //===== Variables for the fits.. none of them is initialized =====
  Double_t pStart[4];
  double parFit_temp[4], err_temp[4];
  Double_t chi[2];
  Int_t fitStatus[2];
  Double_t arglist[10];
  Int_t iflag;
  int nvpar,nparx;
  double amin,edm, errdef;
  //===== Variables only filled when trackNbr_FINAL>=1 =====
  vector<double> lenght, chargeTot;
  vector<double> parFit1, parFit2, parFit3, parFit4, errFit1_local, errFit2_local, errFit3_local, errFit4_local;
  vector<double> parFit1_global, parFit2_global, parFit3_global, parFit4_global, errFit1_global, errFit2_global, errFit3_global, errFit4_global; //not used
  //===== variables for all the MINOS part... don't understand what they are.. yet =====
  //----------------------------------
  vector<TCanvas*> Filter_canvas;
  vector<TCanvas*> Hough_canvas;
  
  double ChargeBin=0.,maxCharge=0.;
  int filled = 0;
  vector<double> Xpad, Ypad, Qpad, XpadNew, YpadNew, QpadNew, ZpadNew;
  vector<double> XpadTemp, YpadTemp, QpadTemp, ZpadTemp; //zpadtemp still not ini
  vector<int> clusterringboolTemp; //not initialized
  vector<int> clusterringbool;
  vector<int> clusternbr;
  vector<int> clusterpads;
  int Iteration=0;
  int filter_result=0;
  int indexfill=0;
  bool fitbool = false; 
  int fit2DStatus = 0;
  double Chi2=0.; //not initialized
  double x_mm,y_mm,z_mm,q_pad,t_pad, r_mm;
  TF1 *fit_function = new TF1("fit_function",Tracking_functions,&Tracking::conv_fit, 0, 511, 3,"Tracking","conv_fit"); // only 3 param. because baseline is fixed. not ini
  double hfit_max, hfit_max_T, T_min, T_max; //not initialized
  TH1F *hfit = new TH1F("hfit","hfit",512,0,512);
  TGraph* grxztmp=new TGraph(); //not init
  TGraph *gryztmp=new TGraph(); //not init
  //2nd step variables
  int npoint=0;
  int npoint_temp=0, cluster_temp=0; //not ini
  int ringsum=0;
  int ringtouch[18]={0}; //not init
  double zmax=0.;
  int allevt_2pfiltered=0; //not ini
  int array_final=0;
  
  bool padsnbr1=false, padsnbr2=false;
  vector<TGraph > gryz;
  vector<TGraph > grxz;
  vector<double> xin, yin, zin, qin, xout, yout, zout, xoutprime, youtprime, zoutprime, qout; //not ini
  TMinuit *min ;  //not ini
  //  TF1 *fit_result;
  TVector3 point; //not ini
  //TVector3 offset(-1.4,-1.4,-4507.3-75-7);
  TVector3 offset(-1.4,-1.4,-75-7);

  Int_t br56sc, br56ca, br54ca, sa55ca, sa55k, sa53ca;
  
  //Branches definition
  //-----------------------
  tree->Branch("fitdata",&fitdata);
  tree->Branch("MINOSOnline",&MINOSOnline,"MINOSOnline/O");
  tree->Branch("EventNumber",&evtOrig,"EventNumber/I");
  tree->Branch("Ridf_EventNumber",&evenumber,"Ridf_EventNumber/L");
  tree->Branch("RunNumber",&runnumber,"RunNumber/I");
  tree->Branch("NumberTracks",&trackNbr_FINAL,"NumberTracks/I");
  //tree->Branch("data_result",&data_result);
  tree->Branch("VDrift",&VDrift,"VDrift/D");
  tree->Branch("DelayTrig",&DelayTrig,"DelayTrig/D");
  
  tree->Branch("chargeTot","vector<double>",&chargeTot);  
  tree->Branch("lenght","vector<double>",&lenght);
  tree->Branch("parFit1","vector<double>",&parFit1);
  tree->Branch("parFit2","vector<double>",&parFit2);
  tree->Branch("parFit3","vector<double>",&parFit3);
  tree->Branch("parFit4","vector<double>",&parFit4);

  tree->Branch("Pulser_charge",&Pulser_charge,"Pulser_charge/D");

  tree->Branch("MyEventInfo_fBit",&EventInfo_fBit);
  /*
  tree->Branch("aoqSamurai_cut",&aoqSamurai_cut);
  tree->Branch("zetSamurai_cut",&zetSamurai_cut);
  tree->Branch("aoqBR_cut",&aoqBR_cut);
  tree->Branch("zetBR_cut",&zetBR_cut);
  */
  tree->Branch("br56ca",&br56ca);
  tree->Branch("br56sc",&br56sc);
  tree->Branch("br54ca",&br54ca);
  tree->Branch("sa55ca",&sa55ca);
  tree->Branch("sa55k",&sa55k);
  tree->Branch("sa53ca",&sa53ca);
  
  //===== Get parameters for the MINOS ANALYSIS =====
  //-------------------------------------------
  ifstream ConfigFile;
  //ConfigFile.open("../../parameters/MinosVDrift.par");
  ConfigFile.open("/home/koiwai/analysis/db/config_MINOSdrift.txt");
  //int HeaderFile;
  int header;
  while(ConfigFile.is_open())
    {
      //ConfigFile >> HeaderFile >> DelayTrig >> StopT >> VDrift;
      ConfigFile >> header >> DelayTrig >> StopT >> VDrift;
      // cout<<"searching: "<<RidfRunNumber<<" among "<<HeaderFile<<endl;
      //if(HeaderFile == RidfRunNumber)
      if(header == FileNum)
	break;
    }
  ConfigFile.close();

  cout << endl;
  cout << " *** MINOS Configuration Parameters *** " << endl;
  cout << "MINOSthresh = " << MINOSthresh << " (bins)" << endl ;
  cout << "TimeBinElec = " << TimeBinElec << " (ns) "<< endl ;    
  cout << "Tshaping = " << Tshaping << " (ns) " << endl;  
  cout << "DelayTrig = " << DelayTrig << " (ns) "<< endl ;   
  cout << "VDrift = " << VDrift << " (mm/ns) " << endl;
  cout << endl;
  cout<< "-----------------------------------------------"<<endl;
  cout << "Conversion started (Wait a few seconds for display) " << endl;
  cout<< "-----------------------------------------------"<<endl;
  int neve = 0;

  //Main loop over the events
  //----------------------------------------------------------
  while(estore->GetNextEvent()&&neve<1000000)
    //while(estore->GetNextEvent() && neve<100000)
    {
      //if(neve%10000==0)
      //cout << "Event " << neve << "\r" << flush;
      if(neve%100==0) clog << neve/1000 << "k events treated..." << "\r";
      evtOrig = neve;
      //if(neve%1000)tree->AutoSave();
      
      //get the event and run number
      //myInfo->LoadData(); //coin register
      runnumber = ((TArtEventInfo *)info_array->At(0))->GetRunNumber();
      evenumber = ((TArtEventInfo *)info_array->At(0))->GetEventNumber();
      
      //get the trigger distribution
      fEvent = (TArtRawEventObject *)sman->FindDataContainer("RawEvent"); //TArtRawEventObject
      EventInfo_fBit.clear();  //vector
      for(Int_t i=0;i<fEvent->GetNumSeg();i++) //loop over the raw event segment
	{
	  seg = fEvent->GetSegment(i); //TArtRawSegmentObject
	  Int_t detector = seg->GetDetector();
	  if(COIN == detector) //if the data comes from the coincidences register
	    {
	      d = seg->GetData(0); //TArtRawDataObject
	      Int_t val = d->GetVal();
	      if (seg->GetFP()==F13) //if focal plane is 3.. it is in this moment that the condition in never fulfilled.. fp is 13 or 63
		{
		  for (int id=1;id<nch+1;++id) //loop over 1 to 7
		    {
		      if(IsChTrue(id,val)) 
			EventInfo_fBit.push_back(id-1);
		    }
		}
	    }
	}

       //---------------------------------------------------------------------
      //Clear & Reset variables
       fitdata.Clear();
       data_result.Clear();
       trackNbr=0;
       trackNbr_FINAL=0;
       Pulser_charge = 0.;
       /*
       aoqSamurai=-999;
       aoqSamurai_cut=-999;
       zetSamurai=-999;
       zetSamurai_cut=-999;
       aoq37_313=-999;
       zet313=-999;
       aoqBR_cut=-999;
       zetBR_cut=-999;
       */
       lenght.clear();
       chargeTot.clear();
       parFit1.clear();
       parFit2.clear();
       parFit3.clear();
       parFit4.clear();
       errFit1_local.clear();
       errFit2_local.clear();
       errFit3_local.clear();
       errFit4_local.clear();
       
      // Filter_canvas.clear();
      // Hough_canvas.clear();
       ChargeBin = 0.;
       maxCharge = 0.;
       filled=0;
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
       Iteration=0;
       filter_result=0;
       indexfill = 0;
       fitbool = false;
       fit2DStatus=0;
       x_mm = 0.; y_mm = 0.; z_mm = 0.; q_pad = 0.; t_pad = 0.;
       hfit->Reset();

       npoint=0;
       ringsum=0;
       zmax=0.;
       array_final=0;
       padsnbr1=false;
       padsnbr2=false;
       grxz.clear();
       gryz.clear();

       br56ca = 0;
       br56sc = 0;
       br54ca = 0;
       sa55ca = 0;
       sa55k = 0;
       sa53ca = 0;
       
       //Making MINOS Reconstruction
       CalibMINOS->ClearData();
       CalibMINOS->ReconstructData();
       //Get the pulser charge
       for(Int_t i=0;i<CalibMINOS->GetNumCalibMINOS();i++) //loop over pads
	 {
	   minos = CalibMINOS->GetCalibMINOS(i); 
	   double maxCharge = 0.;
	   if(minos->GetDetID() == 3 && minos->GetAsic()==0 && minos->GetChannel()==64) //there is never detID==3
	     {
	       for(Int_t j=0; j<minos->GetNData(); j++) //loop over samples
		 {
		   if(minos->GetCalibValue(j)>maxCharge) maxCharge = minos->GetCalibValue(j);//cout<<maxCharge<<endl;
		 }
	       Pulser_charge = maxCharge;
	     }
	 }
       
       //Cut in BR and ZD
       //----------------
       bool SACutBool = false;
       bool BRCutBool = false;
       //tfrag->GetEntry(neve);
       trb->GetEntry(neve);
       if(cbr56ca->IsInside(aoqBR,zetBR)){
	 BRCutBool = true;
	 br56ca = 1;
       }
       else if(cbr56sc->IsInside(aoqBR,zetBR)){
	 BRCutBool = true;
	 br56sc = 1;
       }
       else if(cbr54ca->IsInside(aoqBR,zetBR)){
	 BRCutBool = true;
	 br54ca = 1;
       }
       if(csa55ca->IsInside(aoqSA,zetSA)){
	 SACutBool = true;
	 sa55ca;
       }
       else if(csa55k->IsInside(aoqSA,zetSA)){
	 SACutBool = true;
	 sa55k;
       }
       else if(csa53ca->IsInside(aoqSA,zetSA)){
	 SACutBool = true;
	 sa53ca = 1;
       }
//       if(samcut1->IsInside(aoqSamurai,zetSamurai)) samCutBool = true;
//       else if(samcut2->IsInside(aoqSamurai,zetSamurai)) samCutBool = true;
//       else if(samcut3->IsInside(aoqSamurai,zetSamurai)) samCutBool = true;
//       else if(samcut4->IsInside(aoqSamurai,zetSamurai)) samCutBool = true;
//       else if(samcut5->IsInside(aoqSamurai,zetSamurai)) samCutBool = true;
//       else if(samcut6->IsInside(aoqSamurai,zetSamurai)) samCutBool = true;
       //if(aoqSamurai>2.5 && aoqSamurai<3.1 && zetSamurai>15.5 && zetSamurai<25 ) samCutBool = true;
       
       if((SACutBool==false)||(BRCutBool==false)){
	   tree->Fill();
	   neve++;
	   continue;
	 }
       //aoqSamurai_cut=aoqSamurai;
       //zetSamurai_cut=zetSamurai;
       //aoqBR_cut=aoq37_313;
       //zetBR_cut=zet313;
       //gcout<<"passed"<<endl;
       
       //ONLINE  tracking algorithm
       //----------------------------
       if(MINOSOnline==true)
	 {
	   AnalyzedMINOS->ClearData();
	   TrackMINOS->ClearData();
	   AnalyzedMINOS->SetConfig(VDrift, TimeBinElec, DelayTrig); //, MINOSthresh);
	   AnalyzedMINOS->ReconstructData();
	   if(AnalyzedMINOS->GetNumAnalyzedMINOS()>10) {
	     TrackMINOS->ReconstructData();
	     //TArtTrackMINOSData *trackdata = NULL;
	     //TArtTrackMINOSData *trackdata2 = NULL;
	     trackNbr_FINAL = TrackMINOS->GetTrackNumMINOS();
	     if(TrackMINOS->GetTrackNumMINOS()>=1) {
               trackdata = TrackMINOS->GetTrackMINOS(0);
               parFit1.push_back(trackdata->GetPar_x0());
               parFit2.push_back(trackdata->GetPar_Ax());
               parFit3.push_back(trackdata->GetPar_y0());
               parFit4.push_back(trackdata->GetPar_Ay());
               if(TrackMINOS->GetTrackNumMINOS()>1) {
		 trackdata2 = TrackMINOS->GetTrackMINOS(1);
		 parFit1.push_back(trackdata2->GetPar_x0());
		 parFit2.push_back(trackdata2->GetPar_Ax());
		 parFit3.push_back(trackdata2->GetPar_y0());
		 parFit4.push_back(trackdata2->GetPar_Ay());
               }
	     }// end if trackMINOS->GetTrackNumMINOS()>=1
	     else {
               trackNbr_FINAL = -1;
	     }                
	   }//loop for AnalyzedMINOS > 10 pads
	 }//end of MINOS Online code (if MINOSOnline == true)   // ONLINE algorithm end
       //--------------------------
       //OFFLINE tracking algorithm 
       //----------------------------------------------------------------------
       else 
	 {
	   TMinosClust *minosfitdata;
	   //----------------------------------------------------------------------------------------------
	   //  IN EVENT LOOP  : STEP1
	   //  Filling vectors with (x,y,q) information for TPC and count the number of filled (valid pads)
	   //-----------------------------------------------------------------------------------------------
	   for(Int_t i=0;i<CalibMINOS->GetNumCalibMINOS();i++) //loop over pads
	     {
	       minos = CalibMINOS->GetCalibMINOS(i); //get pad
	       ChargeBin = 0.;
	       maxCharge = 0.;
	       x_mm = minos->GetX()+offsetx; //get positions
	       y_mm = minos->GetY()+offsety;
	       if(minos->GetDetID() != 0) continue; //2 and 3 are dssd or something else
	       else
		 {
		   if(!(abs(x_mm)<0.0001 && abs(y_mm)<0.0001) )  // NON connected MINOS (TPC) channels...
		     {
		       for(Int_t j=0; j<minos->GetNData(); j++) //loop over samples
			 {
			   if(minos->GetCalibValue(j)>maxCharge)
			     maxCharge = minos->GetCalibValue(j); //get max charge of each pad
			 }
		       if(maxCharge>=MINOSthresh) //if above threshold
			 {
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
	   if(filled>0) //if there are good pads
	     {
	       while(Xpad.size()>=10 && Iteration<20)  //take events until there are less than 10 pads or got 20 clusters
		 {
		   filter_result = 0;
		   Iteration++;
		   XpadTemp.clear();                
		   YpadTemp.clear();                
		   QpadTemp.clear();                
		   clusterringboolTemp.clear();
		   /*
		   if(plots)
		     {
		       Filter_canvas.push_back(new TCanvas(Form("Event%d_cluster%d", neve, Iteration), 
							       Form("Event%d_cluster%d", neve, Iteration)));
		       //filter result tells me the number of pads belonging to a clusters. 
		       filter_result = Tracking_functions->Hough_modified(Filter_canvas.back(), 
									      &Xpad, &Ypad, &Qpad, 
									      &XpadTemp, &YpadTemp, &QpadTemp, 
									      &clusterringboolTemp);
		       
		       Filter_canvas.back()->Write();
		     }
		   */
		   //else //without plots
		       //{
		   filter_result = Tracking_functions->Hough_modified(&Xpad, &Ypad, &Qpad, 
								      &XpadTemp, &YpadTemp, &QpadTemp, 
								      &clusterringboolTemp);
		   //}
		   
		   
		   if(filter_result<0) break; 
		   if(filter_result>10 && clusterringboolTemp.back()==1) //more than 10 pads and more than 2 in the first ring
		     {
		       trackNbr++; //proper track!
		       for(int ik=0; ik<filter_result; ik++) //loop over the number of good pads
			 {
			   XpadNew.push_back(XpadTemp[ik]); //fill the vectors containing the pads that form the cluster  
			   YpadNew.push_back(YpadTemp[ik]);                
			   ZpadNew.push_back(-10000);       // no z pad here..
			   QpadNew.push_back(QpadTemp[ik]);  
			   clusterringbool.push_back(clusterringboolTemp[ik]);                
			   clusternbr.push_back(Iteration); //cluster number
			   clusterpads.push_back(filter_result);
			 }
		     }
		 }
	     }
	   for(unsigned int il=0; il<XpadNew.size(); il++) //loop over the number of new pads
	     {
	       minosfitdata = (TMinosClust*)fitdata.ConstructedAt(il);
	       //Set(Double_t xmm, Double_t ymm,
	       //    Double_t tns, Double_t zmm, Double_t chargemax,
	       //    Double_t ncluster, Double_t npads, Double_t chi2)
	       minosfitdata->Set(XpadNew[il], YpadNew[il],
				 -10000, -10000, QpadNew[il], 
				 clusternbr[il],clusterpads[il], 0.);
	       ZpadNew.push_back(-10000.);
	     }
	
	   //------------------------------------------------------------------
	   //  IN EVENT LOOP  : STEP 3
	   //  Following analysis only if #tracks = 1||2||3||4
	   //------------------------------------------------------------------
	   if(trackNbr>0 && trackNbr<5) 
	     {
	       if(filled==0) cerr << "Error !!!" << endl;
	       //-------------------------------------------------------
	       //  STEP 3.1: Fitting taken pads for Qmax and Ttrig info
	       //-------------------------------------------------------
	       for(Int_t i=0;i<CalibMINOS->GetNumCalibMINOS();i++) //loop over all pads
		 {
		   minos = CalibMINOS->GetCalibMINOS(i);
		   hfit->Reset();
		   fitbool = false;
		   if(minos->GetDetID() != 0) continue; 
		   x_mm = minos->GetX()+offsetx;
		   y_mm = minos->GetY()+offsety;
		   r_mm=sqrt(x_mm*x_mm+y_mm*y_mm);
			       
		   for(unsigned int jj=0; jj<XpadNew.size(); jj++) //get the pads that were good in the previous step...
		     {
		       if( abs(XpadNew[jj]-x_mm)<0.0001 && abs(YpadNew[jj]-y_mm)<0.0001)
			 {
			   fitbool = true;
			   indexfill=jj;
			   break;
			 }
		     }
		   // Check if New channel is of interest
		   //(if so, we read the Q(t) signal, and we should 
		   // fill the vectors w/ t & q information after fitting E(T))
		   if( fitbool==true ) 
		     {
		       for(Int_t j=0; j<minos->GetNData(); j++) //loop over samples
			 {
			   if(minos->GetCalibValue(j)>=0)
			     hfit->SetBinContent(hfit->FindBin(minos->GetCalibTime(j)), minos->GetCalibValue(j)+250);
			 }
		       // Fitting the hfit histogram of last channel if not empty
		       if(hfit->GetSumOfWeights()>0)
			 {
			   hfit->GetXaxis()->SetRange(0,510);
			   hfit_max = hfit->GetMaximum();
			   hfit_max_T = hfit->GetMaximumBin();
			   T_min=-1;
			   T_max=-1;
			   // Find the T_min & T_max limits of the signal non zero
			   for(int h=hfit_max_T;h>0;h--)
			     {
			       if(T_min == -1 && (hfit->GetBinContent(h))<=250 )
				 {
				   T_min = h;
				   break;
				 }
			     }
			   for(int h=hfit_max_T;h<510;h++)
			     {
			       if(T_max == -1 && (hfit->GetBinContent(h))==0 )
				 {
				   T_max = h;
				   break;
				 }
			     }
			   //Take only 1.5*Tshaping before the max if other signals before...
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
			   //gStyle->SetOptFit(1111);
			   double fit_function_max = 0., fit_function_Tpad = 0.;
			   if(fit2DStatus==0)
			     {
			       //fit_result = hfit->GetFunction("fit_function");
			       Chi2 = fit_function->GetChisquare(); //get chi2 of the fit
			       fit_function_max = fit_function->GetMaximum(); //q max
			       fit_function_Tpad = fit_function->GetParameter(1); //tpad
			     }
			
			   //attribute q_pad and z_mm value (z is inside the tpc)
			   if(fit2DStatus!=0 || fit_function_max<=20. || 
			      fit_function_max>=100000. || fit_function_Tpad<=0.15 || 
			      fit_function_Tpad>=513. || fit_function->GetParameter(2)<=0.15 || 
			      fit_function->GetParameter(2)>=513.)
			     {
			       //cout << "NOT CORRECTLY FITTED !!!!!!! chi2 = " << vectout->chi2.back() << endl;
			       q_pad = hfit_max-250.;
			       z_mm = -10000;
			     }
			   else 
			     {
			       // Add to the variables the fit parameters
			       t_pad = fit_function_Tpad;//trigger time
			       
			       double t1=t_pad*TimeBinElec;
			       double t0=DelayTrig;
			       //z_mm = ((t_pad*TimeBinElec-DelayTrig)*VDrift); //time bin*(ns/us)*vdrift(mm/ns) =>mm
			       //-----------------------------
			       z_mm = VDrift*(t1-t0); //time bin*(ns/us)*vdrift(mm/ns) =>mm
			       //one more try
			       double phi=atan2(y_mm,x_mm);
			       double zprime=(z_mm/300)-1.0;
			       double yprime=x_mm*cos(angleRotation)-y_mm*sin(angleRotation);

			       double Er=81-r_mm;
			       double IntegralEz=sqrt(1-zprime*zprime);
			       
			       double deltaR=c0*Er*IntegralEz*(1+c1/c0*sin(phi+angleRotation));
			       //double r_new= r_mm + deltaR;

			       //Chen
			       double kz=0.06;
			       double k=0.2*(1.+kz*(x_mm/r_mm*cos(angleRotation)-y_mm/r_mm*sin(angleRotation)));
			       double r_new = (r_mm-81.95)*(1.-k*sqrt(1-pow(z_mm/300.-1,2)) )+81.95;
			       
			       //double r_new=r_mm+(81.95-r_mm)*0.2;
			       
			       			       
			       double ratio=r_new/r_mm;
			       double x_new=x_mm*ratio;
			       double y_new=y_mm*ratio;
			       
			       XpadNew[indexfill]=x_new;
			       YpadNew[indexfill]=y_new;
			       //XpadNew[indexfill] = x_mm;
			       //YpadNew[indexfill] = y_mm;			       
			       
			       q_pad = fit_function_max-250.;  // Max charge/event
			       
			     }
			   ZpadNew[indexfill] = z_mm; 
			   // Comment out for saving histograms of Q(t)
			   // TH1F *hfit_clone = (TH1F*)hfit->Clone(Form("E(T)_%d_%f_%f",neve,x_mm,y_mm));
			   //hfit_clone->Write();
			   minosfitdata = (TMinosClust*)fitdata.ConstructedAt(indexfill);
			   minosfitdata->Set(XpadNew[indexfill], YpadNew[indexfill], t_pad*TimeBinElec, 
					     ZpadNew[indexfill], q_pad, clusternbr[indexfill], clusterpads[indexfill], Chi2);
			   //Fill the z and q information for next steps (3D Hough filter & 3D fit weighted by charge)
			   QpadNew[indexfill] = q_pad;
			
			 }//end if histogram not empty
		    
		     }// end if fitbool==true
		   else continue;
		 }//END of entries in tclonesarray for the event
	       
	       
	       //-------------------------------------------------------
	       //  STEP 3.2:  Filtering the tracks off possible noise 
	       //             with Hough3D (3*2D planes)
	       //-------------------------------------------------------
	       
	       for(unsigned int i=0;i<(XpadNew.size());i++)
		 {
		   //cerr << "Event " << neve << ", xpadnew=" << i << endl;
		   if(xin.size()>0 && ((cluster_temp!=int(clusternbr[i]) && i!=0) || i==(XpadNew.size() - 1)))
		     {
		       /*
		       if(plots)
			 {
			   Hough_canvas.push_back(new TCanvas(Form("HEvent%d_cluster%d", neve, cluster_temp), 
							      Form("HEvent%d_cluster%d", neve, cluster_temp)));
			   Tracking_functions->Hough_3D(Hough_canvas.back(), &xin, &yin, &zin, &qin,
							&xout, &yout, &zout, &qout);
			   
			   Hough_canvas.back()->Write();
			 }
		       */
		       //else //without plots
			 Tracking_functions->Hough_3D(&xin, &yin, &zin, &qin, &xout, &yout, &zout, &qout);
		       
		       //xout is the vector of valid pads after the second Hough transformation
		       for(unsigned int ij=0; ij<xout.size();ij++)
			 {
			   if(zout[ij]>zmax) zmax = zout[ij];
			   ringtouch[int((sqrt(xout[ij]*xout[ij]+yout[ij]*yout[ij])-45.2)/2.1)]++; //ad the pad to the ring
			 }
		       for(int ko=0; ko<18; ko++)
			 {
			   if(ringtouch[ko]>0) ringsum++; //count how many rings are hit
			 }
		       if(zmax>290) ringsum=16; //dunnoo
		       // Tracks of interest: >10 pads and >=12 rings hit
		       if(xout.size()>10 && ringsum>=12)
			 {
			   npoint=0;
			   trackNbr_FINAL++; //number of tracks double filtered
			   //grxz.push_back(new TGraph());
			   //gryz.push_back(new TGraph());
			   double charge_temp=0.;
			   double lenght_temp=0;
			   double zintarget=300;
			   double zouttarget=0;
			   int indexin=0; int indexout=0;
			   //cout<<"xout size = "<<xout.size()<<endl;
			   for(unsigned int ij=0; ij<xout.size(); ij++) //loop over good pads
			     {
			       //grxz.back()->SetPoint(npoint,zout[ij],xout[ij]);
			       //  gryz.back()->SetPoint(npoint,zout[ij],yout[ij]);
			       //cout<<"before "<<xout[ij]<<" "<<yout[ij]<<" "<<zout[ij]<<endl;
			       point.SetXYZ(xout[ij],yout[ij],zout[ij]);  //vector to the pad
			       //point.RotateY(0.0064);
			       //point.RotateX(0.004);
			       //point=point+offset;
			       //-------------------------ANGLE!!!
			       ///point.RotateZ(15.6*TMath::DegToRad());//15.6 ROTATION!!!!
			       
			       point.RotateZ(angleRotation);
			       
			       xoutprime.push_back(point.X()); //rotated points
			       youtprime.push_back(point.Y());
			       zoutprime.push_back(point.Z());
			       //cout<<"after "<<point.X()<<" "<<point.Y()<<" "<<point.Z()<<endl<<endl;
			       grxztmp->SetPoint(npoint,zoutprime[ij],xoutprime[ij]); //create graphs of the xz and yz planes
			       gryztmp->SetPoint(npoint,zoutprime[ij],youtprime[ij]);
			       charge_temp += qout[ij]; //total charge 
			       minosdata_result = (TMinosResult*)data_result.ConstructedAt(array_final);
			       minosdata_result->Set(xoutprime[ij], youtprime[ij], zoutprime[ij], qout[ij], trackNbr_FINAL, xout.size(), zmax);
			       array_final++; 
			       npoint++;
			       if(zoutprime[ij]<zintarget)
				 {
				   zintarget=zoutprime[ij];
				   indexin=ij;
				 }
			       if(zoutprime[ij]>zouttarget)
				 {
				   zouttarget=zoutprime[ij];
				   indexout=ij;
				 }
			     }
			   //distance inside the tpc
			   if(xout.size()>0)
			     lenght_temp=sqrt(pow(zoutprime[indexin]-zoutprime[indexout],2)+pow(youtprime[indexin]-youtprime[indexout],2)+pow(xoutprime[indexin]-xoutprime[indexout],2));

			   grxz.push_back(*grxztmp);
			   gryz.push_back(*gryztmp);
			   grxztmp->Set(0);
			   gryztmp->Set(0);
			   //cout<<"L = "<<lenght_temp<<" --- Q = "<<charge_temp<<endl;
			   chargeTot.push_back(charge_temp);
			   lenght.push_back(lenght_temp);
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
		       npoint_temp=0;
		       ringsum=0;// angle between 1st track and z axis in 3D in degrees
		       zmax=0.;
		       for(int ko=0; ko<18; ko++) ringtouch[ko] = 0;

		     }

		   cluster_temp = clusternbr[i];

		   if(!(clusterpads[i]>=10 && clusterringbool[i]==1 && ZpadNew[i]>-10000 && ZpadNew[i]<=320)) continue;
		   else
		     {
		       xin.push_back(XpadNew[i]);
		       yin.push_back(YpadNew[i]);
		       zin.push_back(ZpadNew[i]);
		       qin.push_back(QpadNew[i]);
		       npoint_temp++;
		     }

		 }//end of loop on pads


	       //-------------------------------------------------------
	       //  STEP 3.3:  Fitting the filtered tracks in 3D 
	       //             (weight by charge, TMinuit)
	       //-------------------------------------------------------
	       
	       //For 1 track found or more (less than 5 in total)
	       if(trackNbr_FINAL>=1)
		 {
		   //////////Minimization in 3D to reconstruct track lines
		   allevt_2pfiltered++;
		   for(int itr=0; itr<trackNbr_FINAL; itr++)
		     {
		       pStart[0]=0; pStart[2]=0; pStart[1]=1; pStart[3]=3;
		       //cout<<"start Minuit "<<endl;
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
		       for(int i = 0; i <4; i++) min->GetParameter(i,parFit_temp[i],err_temp[i]);
		       
		       parFit1.push_back(parFit_temp[0]);
		       parFit2.push_back(parFit_temp[1]);
		       parFit3.push_back(parFit_temp[2]);
		       parFit4.push_back(parFit_temp[3]);
		       //errFit1_local.push_back(err_temp[0]);
		       //errFit2_local.push_back(err_temp[1]);
		       //errFit3_local.push_back(err_temp[2]);
		       //errFit4_local.push_back(err_temp[3]);
		       //offset SAMURAI-mid target+mid target-target entrance+target entrance+TPC MM plane (z=0 in local coord)
		       /*parFit1_global.push_back(parFit_temp[0]+parFit_temp[1]*(4507+75+7)); 
			 parFit2_global.push_back(parFit_temp[1]);
			 parFit3_global.push_back(parFit_temp[2]+parFit_temp[3]*(4507+75+7));
			 parFit4_global.push_back(parFit_temp[3]);
			 errFit1_global.push_back(err_temp[0]);
			 errFit2_global.push_back(err_temp[1]);
			 errFit3_global.push_back(err_temp[2]);
			 errFit4_global.push_back(err_temp[3]);
		       */
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
	 }// end of MINOS Offline code (if MINOSOnline = false)
       //---------------------------------------------------------------------------
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

     time(&stop);
     printf("Elapsed time: %.1f seconds\n",difftime(stop,start));

     return 0;

   }
   
/// Functions to be minimized
void SumDistance(int &, double *, double & sum, double * par,  int)
{
  int nused=0;
  double qtot=0;
  sum = 0;
  //cout<<"sum "<<sum<<" over "<<npoints<<endl;
  //double factor;
  //cout<<"*************after fit "<<endl;
  for(int i=0; i<data_result.GetEntriesFast(); i++)
    {
      minosdata_result = (TMinosResult*)data_result.At(i);
      if(minosdata_result->n_Cluster==NclusterFit)
	{
	  float x=minosdata_result->x_mm;
	  float y=minosdata_result->y_mm;
	  float z=minosdata_result->z_mm;
	  float q=minosdata_result->Chargemax;
	  //if(nused<2)cout<<minosdata_result->n_Cluster<<" "<<x<<" "<<y<<" "<<z<<" "<<q<<endl;
	  double d = Tracking_functions->distance2(x, y, z, par);
	  sum += d*q;
	  nused++;
	  qtot+=q;
	}
    }
  //sum/=nused;
  sum/=qtot;
  return;
}

bool IsChTrue(Int_t id, Int_t val)
{
  //     cout<<"hallo "<<std::bitset<8>(val)<<endl;
  return ((1 << (id - 1)) & val) ? true : false;
}

