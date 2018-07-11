void vcalc(int runnum){
  TFile *driftV = TFile::Open(Form("/home/koiwai/analysis/rootfiles/minos/driftV/driftVminos%04d.root",runnum));
  dvtr->Draw("t_ns>>h(400,1100,13000)");
  double y[400] = {0};
  double t[400] = {0};
  for(int i=0;i<400;i++){
    y[i] = abs(h->GetBinContent(i+1) - h->GetBinContent(i));
    t[i] = 1100 + 30*i + 15;
  }
  t[0] = 1100 + 15;
  y[0] = 0;
  int n = 400;
  TGraph *gr = new TGraph(n,t,y);
  gr->SetName("gr");
  gr->SetTitle("differencial of Tpad");
  TFile *f = new TFile("test.root","recreate");
  gr->Write();
  gr->Draw("");
  

}
