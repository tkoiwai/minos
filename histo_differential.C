#include <stdio.h>
#include <iostream>
#include <fstream>
#include "TH1.h"
#include "TH1F.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TCanvas.h"


void histo_differential(const char* histname){
  TH1F *h = (TH1F*)gROOT->FindObject(histname);
  int numbin = h->GetNbinsX();
  //cout << numbin << endl;
  int xmin   = h->GetXaxis()->GetXmin();
  int xmax   = h->GetXaxis()->GetXmax();
  
  int bin_content[numbin] = {0};

  TH1F *h_diff = new TH1F("h_diff","h_diff",numbin,xmin,xmax);
  
  for(int i=0;i<numbin;i++){
    bin_content[i] = abs(h->GetBinContent(i+1) - h->GetBinContent(i));
    h_diff->SetBinContent(i,bin_content[i]);
  }

  TCanvas *cc;
  cc = new TCanvas("cc","cc",1200,400);
  cc->Divide(2,1);
  cc->cd(1);
  h->Draw();
  cc->cd(2);
  h_diff->Draw();
  
  
  

			 
}
