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
  vector<double> *p0xz;
  vector<double> *p1xz;
  vector<double> *p0yz;
  vector<double> *p1yz;
  Int_t tracknum;
  
  //===== SetBranchAddress =====
  caltrM->SetBranchAddress("EventNumber",&EventNumber_minos);
  caltrM->SetBranchAddress("RunNumber",&RunNumber_minos);

  caltrM->SetBranchAddress("parFit1",&p0xz);
  caltrM->SetBranchAddress("parFit2",&p1xz);
  caltrM->SetBranchAddress("parFit3",&p0yz);
  caltrM->SetBranchAddress("parFit4",&p1yz);
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

  //===== Declare variables =====
  Double_t bdc_dx, bdc_dy;
  //vector<Double_t> tmpx, tmpy, tmpz;
  Double_t tmpx[4], tmpy[4], tmpz[4];
  Double_t tmpx_ave, tmpy_ave;

  //===== Declare tree variables =====
  Int_t runnum, eventnum;

  Double_t vertexX, vertexY, vertexZ;

  //===== Create tree Branch =====
  tr->Branch("runnum",&runnum);
  tr->Branch("eventnum",&eventnum);

  tr->Branch("vertexX",&vertexX);
  tr->Branch("vertexY",&vertexY);
  tr->Branch("vertexZ",&vertexZ);

//===== Begin LOOP =====
  int nEntry = caltrM->GetEntries();
  for(int iEntry=0;iEntry<nEntry;iEntry++){

    if(iEntry%100==0) clog << iEntry/1000 << "k events treated..." << "\r";

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
    
    //=== Calc ===
    if(tracknum==0) continue;
    else{
      for(Int_t i=0;i<tracknum;i++){
	tmpx.push_back((p0xz->at(i)-p0yz->at(i)+p1yz->at(i)*(bdc_dy/bdc_dx*BDC_X-BDC_Y))/(p1xz->at(i)+bdc_dy/bdc_dx*p1yz->at(i)));
	tmpx_ave += tmpx.at(i);
	tmpy.push_back((bdc_dy/bdc_dx*(p0xz->at(i)-p0yz->at(i))+p1xz->at(i)*(BDC_Y-bdc_dy/bdc_dx*BDC_X))/(p1xz->at(i)+bdc_dy/bdc_dx*p1yz->at(i)));	
	tmpy_ave += tmpy.at(i);
      }
      vertexX = tmpx_ave/tracknum;
      vertexY = tmpy_ave/tracknum;

      vertexZ = p0xz->at(0) + p1xz->at(0)*vertexX; 



    }

    








    
    //===== BG cut =====


    tr->Fill();
  }//for LOOP end
  outf->cd();
  tr->Write();
  outf->Close();

  time(&stop);
  printf("Elapsed time: %.1f seconds\n",difftime(stop,start));
}//main()
