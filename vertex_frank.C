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
using namespace std;
using namespace TMath;

int main(int argc, char *argv[]){
  
  time_t start, stop;
       time(&start);
  
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
  vector<double> *p0 = 0; // x = p0 + p1*z
  vector<double> *p1 = 0;
  vector<double> *p2 = 0; // y = p2 + p3*z
  vector<double> *p3 = 0;
  Int_t tracknum;

  //===== SetBranchAddress =====
  caltrM->SetBranchAddress("EventNumber",&EventNumber_minos);
  caltrM->SetBranchAddress("RunNumber",&RunNumber_minos);

  caltrM->SetBranchAddress("parFit1",&p0);
  caltrM->SetBranchAddress("parFit2",&p1);
  caltrM->SetBranchAddress("parFit3",&p2);
  caltrM->SetBranchAddress("parFit4",&p3);
  caltrM->SetBranchAddress("NumberTracks",&tracknum);

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

  //===== Load input Beam file =====
  TString infnameB = Form("/home/koiwai/analysis/rootfiles/ana/beam/ana_beam%04d.root",filenum);
  TFile   *infileB = TFile::Open(infnameB);
  TTree   *anatrB;
  infileB->GetObject("anatrB",anatrB);


  Double_t zetBR, aoqBR;

  anatrB->SetBranchAddress("zetBR",&zetBR);
  anatrB->SetBranchAddress("aoqBR",&aoqBR);


  //===== Load input smri file =====
  TString infnameS = Form("/home/koiwai/analysis/rootfiles/ana/smri/ana_smri%04d.root",filenum);
  TFile   *infileS = TFile::Open(infnameS);
  TTree   *anatrS;
  infileS->GetObject("anatrS",anatrS);


  Double_t zetSA, aoqSA;

  anatrS->SetBranchAddress("zetSA",&zetSA);
  anatrS->SetBranchAddress("aoqSA",&aoqSA);
  
  
  //===== AddFriend =====
  caltrM->AddFriend(anatrDC);
  caltrM->AddFriend(anatrB);
  caltrM->AddFriend(anatrS);
  

  //===== Load cut files =====

  //===== Load .dat files =====
  
  //===== Create output file/tree =====
  TString ofname = Form("/home/koiwai/analysis/rootfiles/minos/vertex/vertex%04d.root",filenum);
  //TString ofname = Form("/home/koiwai/analysis/minos/vertex%04dtest.root",filenum);
  TFile   *outf  = new TFile(ofname,"RECREATE");
  TTree   *tr    = new TTree("tr","tr");

  //===== Declare const.s =====
  Double_t d = 2310.31; // [mm]: BDC1 - MINOS entrance.
  Double_t dBDC = 999.53 ; // [mm]: BDC1 - BDC2.
  Double_t PI = TMath::Pi();

  //===== Declare variables =====
  Double_t bdc_dx, bdc_dy;

  Double_t p0beam, p1beam, p2beam, p3beam;

  

  //===== Declare tree variables =====
  Int_t runnum, eventnum;

  Double_t vertexX, vertexY, vertexZ;
  Double_t vertexR, vertexTheta, vertexPhi;
  
  Double_t p0a, p1a, p2a, p3a;
  Double_t p0b, p1b, p2b, p3b;
  Double_t dp0, dp2, beta, alpha, B, C;
  Double_t xa, xb, ya, yb, za, zb;
  
  Double_t theta2p;

  //===== Create tree Branch =====
  tr->Branch("runnum",&runnum);
  tr->Branch("eventnum",&eventnum);

  //tr->Branch("NumberTracks","tracknum");

  tr->Branch("vertexX",&vertexX);
  tr->Branch("vertexY",&vertexY);
  tr->Branch("vertexZ",&vertexZ);
  tr->Branch("vertexR",&vertexR);
  tr->Branch("vertexTheta",&vertexTheta);
  tr->Branch("vertexPhi",&vertexPhi);

  
  tr->Branch("theta2p",&theta2p);

  //===== Begin LOOP =====
  int nEntry = caltrM->GetEntries();
 for(int iEntry=0;iEntry<nEntry;iEntry++){
  //for(int iEntry=0;iEntry<10000;iEntry++){

    if(iEntry%100==0) clog << iEntry/1000 << "k events treated..." << "\r";
    caltrM->GetEntry(iEntry);
    
    runnum   = RunNumber_minos;
    eventnum = EventNumber_minos;

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

    tr->Fill();
  }//for LOOP end
  outf->cd();
  tr->Write();
  outf->Close();

  time(&stop);
  printf("Elapsed time: %.1f seconds\n",difftime(stop,start));
}//main()
