#include <stdio.h>
#include <iostream>
#include <fstream>
#include "TH1.h"
#include "TH1F.h"
#include "TSystem.h"
#include "TROOT.h"


void histo_differential(int runnum){

  TFile *f = TFile::Open(Form("/home/koiwai/analysis/rootfiles/minos/vertex/vertex%04d.root",runnum));
  tr->Draw("vertexZ>>h(100,-100,200)");
  //TH1F *h = (TH1F*)gROOT->FindObject(histname);
  /*
  int numbin = h->GetNbinsX();
  //cout << numbin << endl;
  int xmin   = h->GetXaxis()->GetXmin();
  int xmax   = h->GetXaxis()->GetXmax();
  */
  int numbin = 100;
  int xmin   = -100;
  int xmax   = 200;
  
  int bin_content[numbin] = {0};

  TH1F *h_diff = new TH1F("h_diff","h_diff",numbin,xmin,xmax);
  
  for(int i=0;i<numbin;i++){
    bin_content[i] = abs(h->GetBinContent(i+1) - h->GetBinContent(i));
    h_diff->SetBinContent(i,bin_content[i]);
  }
  h_diff->Draw();
  
  TF1 *func_min = new TF1("func_min","gaus",-30,-10);
  func_min->SetParameter(0,400);
  func_min->SetParamater(1,-20);
  func_min->SetParameter(2,5);
  func_min->SetParLimits(0,-25,-15);
  func_min->SetParLimits(1,)
  

			 
}
