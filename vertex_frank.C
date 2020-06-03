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
#include "TArtMINOSParameters.hh"
#include "TArtCalibMINOS.hh"
#include "TArtCalibMINOSData.hh"
#include "TArtAnalyzedMINOS.hh"
#include "TArtTrackMINOS.hh"
#include "TArtTrackMINOSData.hh"
#include "TArtVertexMINOS.hh"
#include "TArtEventInfo.hh"
#include "TTree.h"
#include "TFile.h"
#include "TClonesArray.h"
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
#include "TVector3.h"
#include "TArtRawFeminosDataObject.hh"
#include "TCutG.h"
#include "/home/koiwai/analysis/include/liboffline/Tracking.h"
#include "/home/koiwai/analysis/include/liboffline/TMinosClust.h"
#include "/home/koiwai/analysis/include/liboffline/TMinosResult.h"

#include "/home/koiwai/analysis/include/time_tk.h"

using namespace std;
using namespace TMath;

void Vertex(double *p, double *pp, double &xv,double &yv,double &zv, double &min_dist);

int main(int argc, char *argv[]){
  
  initiate_timer_tk();
  
  Int_t filenum = TString(argv[1]).Atoi();

  //---------------------------------
  //===== Load input MINOS file =====
  //TString infnameM = Form("/home/koiwai/analysis/rootfiles/minos/cal/cal_minos%04d.root",filenum);
  TString infnameM = Form("/home/koiwai/analysis/minos/Tracks_run_%04d.root",filenum);
  TFile   *infileM = TFile::Open(infnameM);
  TTree   *caltrM;
  infileM->GetObject("tree",caltrM);

  //===== Input tree variables =====
  Int_t EventNumber_minos, RunNumber_minos;
  vector<double> *parFit1 = 0; // x = p0 + p1*z
  vector<double> *parFit2 = 0;
  vector<double> *parFit3 = 0; // y = p2 + p3*z
  vector<double> *parFit4 = 0;
  Int_t NumberTracks;

  //===== SetBranchAddress =====
  caltrM->SetBranchAddress("EventNumber",&EventNumber_minos);
  caltrM->SetBranchAddress("RunNumber",&RunNumber_minos);

  //caltrM->SetBranchAddress("parFit1",&p0);
  caltrM->SetBranchAddress("parFit1",&parFit1);
  caltrM->SetBranchAddress("parFit2",&parFit2);
  caltrM->SetBranchAddress("parFit3",&parFit3);
  caltrM->SetBranchAddress("parFit4",&parFit4);
  caltrM->SetBranchAddress("NumberTracks",&NumberTracks);

  //------------------------------
  //===== Load input DC file =====
  TString infnameDC = Form("/home/koiwai/analysis/rootfiles/ana/mwdc/anaDC%04d.root",filenum);
  TFile   *infileDC = TFile::Open(infnameDC);
  TTree   *anatrDC;
  infileDC->GetObject("anatrDC",anatrDC);

  //===== Input tree variables =====
  Long64_t EventNumber_dc;
  Int_t RunNumber_dc;

  Double_t BDC1_X, BDC1_Y;
  Double_t BDC2_X, BDC2_Y;
  Double_t Target_X, Target_Y;
  
  //===== SetBranchAddress =====
  anatrDC->SetBranchAddress("EventNumber",&EventNumber_dc);
  anatrDC->SetBranchAddress("RunNumber",&RunNumber_dc);

  anatrDC->SetBranchAddress("BDC1_X",&BDC1_X);
  anatrDC->SetBranchAddress("BDC1_Y",&BDC1_Y);
  anatrDC->SetBranchAddress("BDC2_X",&BDC2_X);
  anatrDC->SetBranchAddress("BDC2_Y",&BDC2_Y);
  anatrDC->SetBranchAddress("Target_X",&Target_X);
  anatrDC->SetBranchAddress("Target_Y",&Target_Y);


  //===== AddFriend =====
  caltrM->AddFriend(anatrDC);

  //===== Load cut files =====

  //===== Load .dat files =====
  
  //===== Create output file/tree =====
  //TString ofname = Form("/home/koiwai/analysis/rootfiles/minos/vertex/vertex%04d.root",filenum);
  TString ofname = Form("/home/koiwai/analysis/minos/vertex%04dtest.root",filenum);
  TFile   *outf  = new TFile(ofname,"RECREATE");
  TTree   *tr    = new TTree("tr","tr");

  //===== Declare const.s =====
  Double_t d = 2310.31; // [mm]: BDC1 - MINOS entrance.
  Double_t Dist_BDC1TPC = 2310.31; // [mm]: BDC1 - MINOS entrance.
  Double_t dBDC = 999.53 ; // [mm]: BDC1 - BDC2.
  Double_t Dist_BDC1BDC2 = 999.53 ; // [mm]: BDC1 - BDC2.
  Double_t PI = TMath::Pi();

  //===== Declare variables =====
  Double_t bdc_dx, bdc_dy;
  Double_t p0beam, p1beam, p2beam, p3beam;

  //===== Declare tree variables =====
  Int_t RunNumber, EventNumber;

  Double_t vertexX, vertexY, vertexZ;
  Double_t vertexR, vertexTheta, vertexPhi;
  
  Double_t p0a, p1a, p2a, p3a;
  Double_t p0b, p1b, p2b, p3b;
  Double_t dp0, dp2, beta, alpha, B, C;
  Double_t xa, xb, ya, yb, za, zb;
  
  Double_t theta2p;

  //===== Declare variables for Frank's part =====

  Double_t parTrack1[4], parTrack2[4], parTrackBDC[4];
  Double_t x_vertexBDC, y_vertexBDC, z_vertexBDC;
  Double_t x_vertex, y_vertex, z_vertex;
  Double_t MINOS_X_BDC, MINOS_Y_BDC, MINOS_Z_BDC;
  Double_t MINOS_X, MINOS_Y, MINOS_Z, MINOS_D_min;
  Double_t Dist_minBDC, Dist_min;
  Double_t MINOS_tr1_phi, MINOS_tr1_theta, MINOS_tr2_phi, MINOS_tr2_theta, MINOS_tr_phi, MINOS_tr_theta;
  Double_t MINOS_Radius;
  Int_t MINOS_NumberTracks;
  Double_t Target_R, TargetLength;


  //===== Create tree Branch =====
  tr->Branch("RunNumber",&RunNumber);
  tr->Branch("EventNumber",&EventNumber);

  tr->Branch("vertexX",&vertexX);
  tr->Branch("vertexY",&vertexY);
  tr->Branch("vertexZ",&vertexZ);
  tr->Branch("vertexR",&vertexR);
  tr->Branch("vertexTheta",&vertexTheta);
  tr->Branch("vertexPhi",&vertexPhi);
  
  tr->Branch("theta2p",&theta2p);

//===== Begin LOOP =====

  prepare_timer_tk();

  int nEntry = caltrM->GetEntries();
  for(int iEntry=0;iEntry<nEntry;iEntry++){
   //for(int iEntry=0;iEntry<10000;iEntry++){

    start_timer_tk(iEntry,nEntry,1000);
    
    caltrM->GetEntry(iEntry);
    
    RunNumber   = RunNumber_minos;
    EventNumber = EventNumber_minos;

    //=== Initialization =====
    vertexX     = Sqrt(-1);
    vertexY     = Sqrt(-1);
    vertexZ     = Sqrt(-1);
    vertexR     = Sqrt(-1);
    vertexTheta = Sqrt(-1);
    vertexPhi   = Sqrt(-1);

    bdc_dx = BDC2_X - BDC1_X;
    bdc_dy = BDC2_Y - BDC1_Y;

    p0beam = Sqrt(-1);
    p1beam = Sqrt(-1);    
    p2beam = Sqrt(-1);
    p3beam = Sqrt(-1);

    p0a = Sqrt(-1);
    p1a = Sqrt(-1);
    p2a = Sqrt(-1);
    p3a = Sqrt(-1);
    p0b = Sqrt(-1);
    p1b = Sqrt(-1);
    p2b = Sqrt(-1);
    p3b = Sqrt(-1);

    dp0   = Sqrt(-1);
    dp2   = Sqrt(-1);
    beta  = Sqrt(-1);
    alpha = Sqrt(-1);
    B     = Sqrt(-1);
    C     = Sqrt(-1);

    xa = Sqrt(-1);
    xb = Sqrt(-1);
    ya = Sqrt(-1);
    yb = Sqrt(-1);
    za = Sqrt(-1);
    zb = Sqrt(-1);

    theta2p = Sqrt(-1);
    
    //=== Calc ===
    p0beam = BDC1_X + bdc_dx/dBDC*d;
    p1beam = bdc_dx/dBDC;
    p2beam = BDC1_Y + bdc_dy/dBDC*d;
    p3beam = bdc_dy/dBDC;

    double p[4], pp[4];
    
/*
    if(tracknum==1){  
      p0a = p0beam;
      p1a = p1beam;
      p2a = p2beam;
      p3a = p3beam;
      p0b = p0->at(0);
      p1b = p1->at(0);
      p2b = p2->at(0);
      p3b = p3->at(0);  
    }
    else if(tracknum>=2){
      p0a = p0->at(1);
      p1a = p1->at(1);
      p2a = p2->at(1);
      p3a = p3->at(1);
      p0b = p0->at(0);
      p1b = p1->at(0);
      p2b = p2->at(0);
      p3b = p3->at(0);  
    }
    else{
      tr->Fill();
      continue;
    }
    
    dp0 = p0b - p0a;
    dp2 = p2b - p2a;
    alpha = -(p1b*dp0+p3b*dp2)/(p1b*p1b+p3b*p3b+1);
    beta  = (p1a*p1b+p3a*p3b+1)/(p1b*p1b+p3b*p3b+1);
    B = (p1a*p1a+p3a*p3a+1) - beta*(p1a*p1b+p3a*p3b+1);
    C = beta*(p1b*dp0+p3b*dp2) - p1a*dp0 - p3a*dp2;

    za = -C/B;
    zb = beta*za + alpha;
    xa = p0a + p1a*za;
    xb = p0b + p1b*zb;
    ya = p2a + p3a*za;
    yb = p2b + p3b*zb;

    vertexX = (xa + xb)/2.;
    vertexY = (ya + yb)/2.;
    vertexZ = (za + zb)/2.;

    vertexR     = Sqrt(vertexX*vertexX + vertexY*vertexY + vertexZ*vertexZ);
    vertexTheta = ACos(vertexZ/vertexR) *180./PI;
    vertexPhi   = ACos(vertexX/vertexR/Sin(vertexTheta)) *180./PI;

    theta2p = ACos((p1a*p1b + p3a*p3b +1)/(Sqrt(p1a*p1a+p3a*p3a+1)*Sqrt(p1b*p1b+p3b*p3b+1)))*180./PI;
*/


    //Frank

     // MINOS 
    if ( parFit1->size() > 0 ){
      parTrack1[0]=parFit1->at(0); // first track in MINOS
      parTrack1[1]=parFit2->at(0); // first track in MINOS
      parTrack1[2]=parFit3->at(0); // first track in MINOS
      parTrack1[3]=parFit4->at(0); // first track in MINOS
      // calc track from BDC
      parTrackBDC[0] = BDC1_X + Dist_BDC1TPC / Dist_BDC1BDC2 * (BDC2_X - BDC1_X);
      parTrackBDC[2] = BDC1_Y + Dist_BDC1TPC / Dist_BDC1BDC2 * (BDC2_Y - BDC1_Y);
      parTrackBDC[1] = (BDC2_X - BDC1_X) / Dist_BDC1BDC2;
      parTrackBDC[3] = (BDC2_Y - BDC1_Y) / Dist_BDC1BDC2;
      Vertex(parTrack1,parTrackBDC,x_vertexBDC,y_vertexBDC,z_vertexBDC,Dist_minBDC);
      MINOS_X_BDC = x_vertexBDC;
      MINOS_Y_BDC = y_vertexBDC;
      MINOS_Z_BDC = z_vertexBDC;
      x_vertex = x_vertexBDC ;
      y_vertex = y_vertexBDC ;
      z_vertex = z_vertexBDC ;
      if(parFit1->size() > 1){
	parTrack2[0]=parFit1->at(1); // second track in MINOS
	parTrack2[1]=parFit2->at(1); // second track in MINOS
	parTrack2[2]=parFit3->at(1); // second track in MINOS
	parTrack2[3]=parFit4->at(1); // second track in MINOS
	Vertex(parTrack1,parTrack2,x_vertex,y_vertex,z_vertex,Dist_min);
	TVector3 ptr1(parTrack1[1],parTrack1[3],1);
	ptr1 = ptr1.Unit();
	MINOS_tr1_theta = ptr1.Theta()/3.1415927*180.;
	MINOS_tr1_phi   = ptr1.Phi()/3.1415927*180.;
	TVector3 ptr2(parTrack2[1],parTrack2[3],1);
	ptr2 = ptr2.Unit();
	MINOS_tr2_theta = ptr2.Theta()/3.1415927*180.;
	MINOS_tr2_phi   = ptr2.Phi()/3.1415927*180.;
	MINOS_tr_theta  = TMath::ACos( ptr1.Dot(ptr2) )/3.1415927*180.;
	MINOS_tr_phi    = TMath::Abs(MINOS_tr1_phi - MINOS_tr2_phi);
	MINOS_tr_phi    = MINOS_tr_phi>180 ? 360-MINOS_tr_phi : MINOS_tr_phi;
      }
      MINOS_X = x_vertex;
      MINOS_Y = y_vertex;
      MINOS_Z = z_vertex;
      MINOS_D_min = Dist_min;
      MINOS_Radius = TMath::Sqrt(x_vertex*x_vertex + y_vertex*y_vertex);
      MINOS_NumberTracks = NumberTracks;
    }
    double TargetLength0 = 150; // in mm
    double tmpk=4./400.;
    Target_R = sqrt(Target_X*Target_X+Target_Y*Target_Y);
    TargetLength = TargetLength0 - tmpk*Target_R*Target_R;





    
    tr->Fill();
  }//for LOOP end
  outf->cd();
  tr->Write();
  outf->Close();

  stop_timer_tk(nEntry);
}//main()



void Vertex(double *p, double *pp, double &xv,double &yv,double &zv, double &min_dist)
{
  double a1 = p[0];
  double a2 = p[2];
  double b1 = p[1];
  double b2 = p[3];
  double ap1 = pp[0];
  double ap2 = pp[2];
  double bp1 = pp[1];
  double bp2 = pp[3];
  double alpha, beta, A, B, C;
  alpha = (bp1*(a1-ap1)+bp2*(a2-ap2))/(bp1*bp1 + bp2*bp2 + 1);
  beta = (bp1*b1+bp2*b2+1)/(bp1*bp1 + bp2*bp2 + 1);
  A = beta*(bp1*bp1 + bp2*bp2 + 1) - (bp1*b1 + bp2*b2 + 1);
  B = (b1*b1 + b2*b2 + 1) - beta*(bp1*b1+bp2*b2+1);
  C = beta*(bp1*(ap1-a1) + bp2*(ap2-a2)) - (b1*(ap1-a1) + b2*(ap2-a2));
  double sol1, solf1;
  double x,y,z,xp,yp,zp;
  sol1 = -(A*alpha + C)/(A*beta + B);
  solf1 = alpha + beta* sol1;
  x = a1 + b1*sol1;
  y = a2 + b2*sol1;
  z = sol1;
  xp = ap1 + bp1*solf1;
  yp = ap2 + bp2*solf1;
  zp = solf1;
  xv = (x+xp)/2.;
  yv = (y+yp)/2.;
  zv = (z+zp)/2.;
  min_dist = sqrt(pow((x-xp),2) + pow((y-yp),2) + pow((z-zp),2));
}
