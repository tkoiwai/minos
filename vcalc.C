void vcalc(int runnum){
  TFile *driftV = TFile::Open(Form("/home/koiwai/analysis/rootfiles/minos/driftV/driftVminos%04d.root",runnum));
  dvtr->Draw("t_ns>>h(400,1100,13000)");
  double y[400] = {0};
  double t[400] = {0};
  for(int i=0;i<400;i++){
    y[i] = abs(h->GetBinContent(i+1) - h->GetBinContent(i));
    t[i] = 1100 + 30*i + 15;
  }
  int tmin_trig, tmax_trig;

  for(int i=1;i<400;i++){
    if(y[i]>200){
      for(int j=i;j<i+20;j++){
	if((y[j]-y[j-1]<0)){
	  tmin_trig = t[j-1];
	  break;
	}
	else continue;
      }
      break;
    }
    else continue;
  }
  for(int i=1;i<400;i++){
    if(y[i]>500){
      for(int j=i;j<i+20;j++){
	if((y[j]-y[j-1]<0)){
	  tmax_trig = t[j-1];
	  break;
	}
	else continue;
      }
      break;
    }
    else continue;
  }

  //cout << "tmin_trig " << tmin_trig << endl;
  //cout << "tmax_trig " << tmax_trig << endl;
  
  t[0] = 1100 + 15;
  y[0] = 0;
  int n = 400;
  TGraph *gr = new TGraph(n,t,y);
  gr->SetName("gr");
  gr->SetTitle("differencial of Tpad");
  TFile *f = new TFile("test.root","recreate");
  gr->Write();
  gr->GetXaxis()->SetRangeUser(tmin_trig-75,tmin_trig+75);
  gr->Draw("apl");
  TF1 *func_min = new TF1("func_min","gaus",tmin_trig-75,tmin_trig+75);
  func_min->SetParameter(0,300);
  func_min->SetParameter(1,tmin_trig);
  func_min->SetParameter(2,30);
  func_min->SetParLimits(0,200,400);
  func_min->SetParLimits(1,tmin_trig-75,tmin_trig+75);
  func_min->SetParLimits(2,0,60);
  gr->Fit("func_min");
  
  gr->GetXaxis()->SetRangeUser(tmax_trig-120,tmax_trig+120);
  gr->Draw("apl");
  TF1 *func_max = new TF1("func_max","gaus",tmax_trig-75,tmax_trig+75);
  func_max->SetParameter(0,700);
  func_max->SetParameter(1,tmax_trig);
  func_max->SetParameter(2,60);
  func_max->SetParLimits(0,500,850);
  func_max->SetParLimits(1,tmax_trig-50,tmax_trig+50);
  func_max->SetParLimits(2,30,90);
  gr->Fit("func_max");
  
  //cout << "Tmin " << func_min->GetParameter(1) << endl;
  //cout << "Tmax " << func_max->GetParameter(1) << endl;
  //cout << "Drift V " << 300/(func_max->GetParameter(1)-func_min->GetParameter(1)) << endl;

  ofstream fout("/home/koiwai/analysis/db/config_MINOSdrift.txt",ios::app);
  fout << runnum << "  " << func_min->GetParameter(1) <<  " " << func_max->GetParameter(1)
       << " " << 300./(func_max->GetParameter(1)-func_min->GetParameter(1)) << endl;


  
}
