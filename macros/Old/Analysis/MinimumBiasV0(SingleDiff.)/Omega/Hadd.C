{
 
  TFile *Lambda = new TFile("PtSpectraOmega.root");
  TFile *AntiLambda = new TFile("PtSpectraAntiOmega.root");
 const int ctr = 5;
  const int eff = 4;
  Double_t cmin [ctr]= {0.0, 5.0, 15.0, 50.0 ,0.0};
  Double_t cmax [ctr]= {5.0, 15.0, 50.0, 100.0, 100.0};
  

  TH1D* LPT[ctr];
  TH1D* ALPT[ctr];

  for (int j = 0 ; j < ctr ; j++) {    
    LPT[j] = (TH1D*) Lambda->Get(Form("hPT_%f-%f",cmin[j],cmax[j]));
    ALPT[j] = (TH1D*) AntiLambda->Get(Form("hPT_%f-%f",cmin[j],cmax[j]));
    //ALPT[j]->SetMarkerStyle(24);
    LPT[j]->Add(ALPT[j]);
  }

  TCanvas* c = new TCanvas();
  c->SetLogy();
  c->SetGridx();
  c->SetGridy();
  LPT[0]->Draw();
  TFile *W = new TFile("OmegaAdd.root","recreate");
  //ALPT[0]->Draw("SAME");
  for (int j = 0 ; j < ctr ; j++) {
    LPT[j]->Draw("SAME");
    LPT[j]->Write();
    //ALPT[j]->Draw("SAME");
  }
 
}
