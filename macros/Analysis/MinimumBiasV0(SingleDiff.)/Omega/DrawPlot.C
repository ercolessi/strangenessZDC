{
  TFile *fP = new TFile("../AntiOmega/AntiOmegaLevyTsallis.root");
  TFile *fAP = new TFile("OmegaLevyTsallis.root");

  TGraphErrors *YieldsNchP = (TGraphErrors *) fP->Get("YieldsNch");
  TGraphErrors *NormYieldsNchP = (TGraphErrors *) fP->Get("NormYieldsNch");

  TGraphErrors *YieldsNchAP = (TGraphErrors *) fAP->Get("YieldsNch");
  TGraphErrors *NormYieldsNchAP = (TGraphErrors *) fAP->Get("NormYieldsNch");
  
  TCanvas * c = new TCanvas("c","c",700,600);
  YieldsNchP->SetName("NormYieldsNch");
  YieldsNchP->SetTitle("#Omega");
  YieldsNchAP->SetTitle("#bar{#Omega}");
  YieldsNchP->GetXaxis()->SetTitle("< n_{ch} >");
  YieldsNchP->GetXaxis()->SetTitleSize(0.05);
  YieldsNchAP->GetYaxis()->SetTitleSize(0.05);
  YieldsNchAP->GetYaxis()->SetTitleOffset(1.4);
  YieldsNchP->GetXaxis()->SetRangeUser(0,25);
  //YieldsNchP->GetYaxis()->SetRangeUser();
  YieldsNchP->SetMarkerStyle(21);
  YieldsNchP->SetMarkerColor(4);
  YieldsNchAP->SetName("NormYieldsNch");
  YieldsNchAP->GetXaxis()->SetTitle("< n_{ch} >");
  YieldsNchAP->GetXaxis()->SetTitleSize(0.05);
  YieldsNchAP->GetYaxis()->SetTitle("dN/dy");
  YieldsNchAP->SetMarkerStyle(25);
  YieldsNchP->SetMarkerSize(0.8);
  YieldsNchAP->SetMarkerSize(1.1);
  YieldsNchAP->SetMarkerColor(4);
  YieldsNchAP->SetLineColor(4);
  YieldsNchP->SetLineColor(4);
  TLegend* l = new TLegend(0.1,0.7,0.48,0.9);
  l->AddEntry(YieldsNchP,"","P");
  l->AddEntry(YieldsNchAP,"","P");
  YieldsNchAP->Draw("AP");
  YieldsNchP->Draw("SAME");
  l->Draw("SAME");



}
