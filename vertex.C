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
using namespace std;
using namespace TMath;

int main(int argc, char *argv[]){
  
  time_t start, stop;
  time(&start);
  
  Int_t filenum = TString(argv[1]).Atoi();

  //---------------------------------
  //===== Load input MINOS file =====
  TString infnameM = Form("/home/koiwai/analysis/rootfiles/minos/cal/cal_minos%04d.root",filenum);
  TFile   *infileM = TFile::Open(infnameM);
  TTree   *caltrM;
  infileM->GetObject("caltrM",caltrM);
  //infileM->Close();

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
  TString infnameDC = Form("/home/koiwai/analysis/rootfiles/ana/mwdc/ana_mwdc%04d.root",filenum);
  TFile   *infileDC = TFile::Open(infnameDC);
  TTree   *anatrDC;
  infileDC->GetObject("anatrDC",anatrDC);
  //infileDC->Close();

  //===== Input tree variables =====
  Int_t EventNumber_dc, RunNumber_dc;

  Double_t BDC1_X, BDC1_Y;
  Double_t BDC2_X, BDC2_Y;
  Double_t BDC_X, BDC_Y, BDC_A, BDC_B;
  Double_t Target_X, Target_Y;

  Int_t BG_flag_dc;

  
  //===== SetBranchAddress =====
  anatrDC->SetBranchAddress("EventNum",&EventNumber_dc);
  anatrDC->SetBranchAddress("RunNum",&RunNumber_dc);

  anatrDC->SetBranchAddress("BDC1_X",&BDC1_X);
  anatrDC->SetBranchAddress("BDC1_Y",&BDC1_Y);
  anatrDC->SetBranchAddress("BDC2_X",&BDC2_X);
  anatrDC->SetBranchAddress("BDC2_Y",&BDC2_Y);
  anatrDC->SetBranchAddress("BDC_X",&BDC_X);
  anatrDC->SetBranchAddress("BDC_Y",&BDC_Y);
  anatrDC->SetBranchAddress("BDC_A",&BDC_A);
  anatrDC->SetBranchAddress("BDC_B",&BDC_B);
  anatrDC->SetBranchAddress("Target_X",&Target_X);
  anatrDC->SetBranchAddress("Target_Y",&Target_Y);
  anatrDC->SetBranchAddress("BG_flag",&BG_flag_dc);

  //===== AddFriend =====
  caltrM->AddFriend(anatrDC);

  //===== Load cut files =====

  
  //===== Load .dat files =====
  
  //===== Create output file/tree =====
  TString ofname = Form("/home/koiwai/analysis/rootfiles/minos/vertex/vertex%04d.root",filenum);
  TFile   *outf  = new TFile(ofname,"RECREATE");
  TTree   *tr    = new TTree("tr","tr");

  //===== Declare const.s =====
  Double_t d = 2233.78; // [mm]: BDC1 - MINOS entrance.
  Double_t dBDC = 998.; // [mm]: BDC1 - BDC2.

  //===== Declare variables =====
  Double_t bdc_dx, bdc_dy;
  vector<Double_t> tmpx, tmpy, tmpz;
  Double_t tmpx_ave, tmpy_ave;

  Double_t p0beam, p1beam, p2beam, p3beam;

  TVector3 a[5], m[5];
  Double_t s[5];

  Double_t p0a, p1a, p2a, p3a;
  Double_t p0b, p1b, p2b, p3b;
  Double_t dp0, dp2, beta, alpha, A, B, C;
  Double_t xa, xb, ya, yb, za, zb;

  
  //===== Declare tree variables =====
  Int_t runnum, eventnum;

  Double_t vertexX, vertexY, vertexZ;

  Double_t lmin;
  
  //===== Create tree Branch =====
  tr->Branch("runnum",&runnum);
  tr->Branch("eventnum",&eventnum);

  tr->Branch("vertexX",&vertexX);
  tr->Branch("vertexY",&vertexY);
  tr->Branch("vertexZ",&vertexZ);

  tr->Branch("lmin",&lmin);
  
  //===== Begin LOOP =====
  int nEntry = caltrM->GetEntries();
  for(int iEntry=0;iEntry<nEntry;iEntry++){

    if(iEntry%100==0) clog << iEntry/1000 << "k events treated..." << "\r";
    //cout << "ok" << endl;
    //sleep(0.1);
    caltrM->GetEntry(iEntry);
    
    runnum   = RunNumber_minos;
    eventnum = EventNumber_minos;

    //=== Initialization =====
    vertexX = Sqrt(-1);
    vertexY = Sqrt(-1);
    vertexZ = Sqrt(-1);
    tmpx_ave = 0;
    tmpy_ave = 0;
    
    tmpx.clear();
    tmpy.clear();
    tmpz.clear();

    bdc_dx = BDC2_X - BDC1_X;
    bdc_dy = BDC2_Y - BDC1_Y;

    p0beam = Sqrt(-1);
    p1beam = Sqrt(-1);    
    p2beam = Sqrt(-1);
    p3beam = Sqrt(-1);

    lmin = Sqrt(-1);

    for(int i=0;i<5;i++) s[i] = Sqrt(-1);

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
    A     = Sqrt(-1);
    B     = Sqrt(-1);
    C     = Sqrt(-1);

    xa = Sqrt(-1);
    xb = Sqrt(-1);
    ya = Sqrt(-1);
    yb = Sqrt(-1);
    za = Sqrt(-1);
    zb = Sqrt(-1);

    
    //=== Calc ===
    p0beam = -(BDC1_X + bdc_dx/dBDC*d);
    p1beam = bdc_dx/dBDC;
    p2beam = -(BDC1_Y + bdc_dy/dBDC*d);
    p3beam = bdc_dy/dBDC;

    /*
    a[0].SetXYZ(-p0beam/p1beam,-q0beam/q1beam,0);
    m[0].SetXYZ(1/p1beam,1/q1beam,1);
    
    for(int i=1;i<tracknum+1;i++){
      a[i].SetXYZ(-p0->at(i-1)/p1->at(i-1),-q0->at(i-1)/q1->at(i-1),0);
      m[i].SetXYZ(1/p1->at(i-1),1/q1->at(i-1),1);
    }
    
    if(tracknum==0){
      tr->Fill();
      continue;
    }
   
    else if(tracknum==1){
      TVector3 b = a[1] - a[0];
      s[0] = (m[1].Mag2()*(b*m[0]) - (m[0]*m[1])*(b*m[0]))/(m[0].Mag2()*m[1].Mag2() - (m[0]*m[1])*(m[0]*m[1]));
      s[1] = (m[0].Mag2()*(b*m[1]) - (m[0]*m[1])*(b*m[1]))/(m[0].Mag2()*m[1].Mag2() - (m[0]*m[1])*(m[0]*m[1]));

      TVector3 l(9999.,9999.,9999.);
      l = b + s[1]*m[1] - s[0]*m[0];
      lmin = l.Mag2();

      TVector3 tmp(9999.,9999.,9999.);
      tmp = a[0] + a[1] + s[0]*m[0] + s[1]*m[1];
      vertexX = 0.5*tmp.X();
      vertexY = 0.5*tmp.Y();
      vertexZ = 0.5*tmp.Z();
    }
    else{
      tr->Fill();
      continue;
    }
    */


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
    else if(tracknum==2){
      p0a = p0->at(1);
      p1a = p1->at(1);
      p2a = p2->at(1);
      p3a = p3->at(1);
      p0b = p0->at(0);
      p1b = p1->at(0);
      p2b = p2->at(0);
      p3b = p3->at(0);  
    }
    else continue;
    
    dp0 = p0b - p0a;
    dp2 = p2b - p2a;
    alpha = -(p1b*dp0+p3b*dp2)/(p1b*p1b+p3b*p3b+1);
    beta  = (p1a*p1b+p3a*p3b+1)/(p1b*p1b+p3b*p3b+1);
    A = beta*(p1b*p1b+p3b*p3b+1) - (p1a*p1b+p3a*p3b+1);
    B = (p1b*p1b+p3b*p3b+1) - beta*(p1a*p1b+p3a*p3b+1);
    C = beta*(p1b*dp0+p3b*dp2) - p1a*dp0 - p3a*dp2;

    za = -(A*alpha+C)/(A*beta+B);
    zb = beta*za + alpha;
    xa = p0a + p1a+za;
    xb = p0b + p1b+zb;
    ya = p2a + p3a+za;
    yb = p2b + p3b+zb;

    vertexX = 0.5*(xa + xb);
    vertexY = 0.5*(ya + yb);
    vertexZ = 0.5*(za + zb);



    






   
    //===== BG cut =====


    tr->Fill();
  }//for LOOP end
  outf->cd();
  tr->Write();
  outf->Close();

  time(&stop);
  printf("Elapsed time: %.1f seconds\n",difftime(stop,start));
}//main()
