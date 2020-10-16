
#include <TFile.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TClonesArray.h>
#include <TLatex.h>
#include <TMath.h>
#include <TColor.h>

// STL includes
#include <iostream>
#include <fstream>
using std::cout;
using std::endl;

TF1 * LevyTsallis(const Char_t *name, Double_t mass, Double_t n = 5., Double_t C = 0.1, Double_t norm = 1.);

void FinalRatios() {

  const char *outputfile = "XiLevyTsallis.root";
  const char *spectrafile = "PtSpectraXi.root"; 

  const int ctr = 9;
  const int eff = 4;
  Double_t cmin [ctr]= {0.0, 1.0, 5.0, 10.0, 15.0 , 30.0, 50.0 , 70.0 ,0.0};
  Double_t cmax [ctr]= {1.0, 5.0, 10.0 , 15.0, 30.0 , 50.0 , 70.0 , 100.0 , 100.0};
  
  double mult[ctr];
  int nevt[ctr];  
  TH1D* hPT[ctr];
  TH1D* hPTClone[ctr];
  double n[ctr];
  double C[ctr];
  double norm[ctr];
  Double_t mass = 1.32131;
  TF1* ffit[ctr];
  TF1* LevyT[ctr];
  Double_t yields[ctr], eyields[ctr];
  TH1D* hmult;
  TH1D* hevt;
  
  //Input File 
  TFile *e = new TFile(spectrafile);
  hmult = (TH1D*) e->Get("hmult");
  hevt = (TH1D*) e->Get("hevt");
  for (int j = 0 ; j < ctr ; j++) { 
    hPT[j] = (TH1D*) e->Get(Form("hPT_%f-%f",cmin[j],cmax[j]));
    hPTClone[j]= (TH1D*)hPT[j]->Clone("hPTClone");
    hPTClone[j]->Reset();
    hPTClone[j]->SetTitle(Form("%f - %f",cmin[j],cmax[j]));
    //inizializzazione variabili
    mult[j]=0;
    nevt[j]=0;
    n[j]=0;
    C[j]=0;
    norm[j]=0;
    ffit[j] = LevyTsallis("LevyTsallis",mass);
    yields[j]=0;
    eyields[j]=0;

    for (int i = 8; i <= hPT[j]->GetNbinsX(); i++){
      double binc = hPT[j]->GetBinContent(i);
      double ebinc = hPT[j]->GetBinError(i);
      hPTClone[j]->SetBinContent(i,binc);
      hPTClone[j]->SetBinError(i,ebinc);
      yields[j] += (hPT[j]->GetBinContent(i)*hPT[j]->GetBinWidth(i));
      eyields[j] += (hPT[j]->GetBinError(i)*hPT[j]->GetBinWidth(i)*hPT[j]->GetBinWidth(i)*hPT[j]->GetBinError(i)); 
    }
    eyields[j] = TMath::Sqrt(eyields[j]);
    hPTClone[j]->SetMarkerStyle(21);
  }
   
  for (int j = 0 ; j < (ctr) ; j++) {   

    hPT[j]->Fit(ffit[j],"","",0.6,10.);  
    n[j] = ffit[j]->GetParameter(1);
    C[j] = ffit[j]->GetParameter(2);
    norm[j] = ffit[j]->GetParameter(3);    
    cout << " Chi Square = " << ffit[j]->GetChisquare()/ffit[j]->GetNDF() << endl;
    cout << " Integral = " << ffit[j]->Integral(0.6,10.)/ffit[j]->Integral(0.,100.) << endl;     

    yields[j] += ffit[j]->Integral(0.0,0.6) + ffit[j]->Integral(10.0,100.0);

    n[j] = hevt->GetBinContent(j+1);
    mult[j] = hmult->GetBinContent(j+1);
        
  } 
  
  Double_t normyields[ctr];
  Double_t enormyields[ctr];
  normyields[ctr-1] = yields[ctr-1]/mult[ctr-1];
  enormyields[ctr-1] = eyields[ctr-1]/mult[ctr-1];    
  Double_t  multScaled[ctr-1];
  Double_t emultScaled[ctr-1];
  Double_t prodScaled[ctr-1];
  Double_t eprodScaled[ctr-1];    
  
  for(int j=0;j < (ctr-1);j++){
    normyields[j] = yields[j]/mult[j]; 
    enormyields[j] = eyields[j]/mult[j];   
    multScaled[j] = mult[j]*1./mult[ctr-1];
    emultScaled[j] = 0;
    prodScaled[j] = normyields[j]/normyields[ctr-1];
    eprodScaled[j] = enormyields[j]/normyields[ctr-1];
  }

  const int npoint = ctr-1;
  TGraphErrors *NormYieldsNormNch = new TGraphErrors(npoint,multScaled,prodScaled,emultScaled,eprodScaled);
  TGraphErrors *NormYieldsNch = new TGraphErrors(npoint,mult,prodScaled,emultScaled,eprodScaled);
  TGraphErrors *YieldsNch = new TGraphErrors(npoint,mult,yields,emultScaled,eyields);

  hPTClone[0]->SetLineColor(kRed+1);
  hPTClone[0]->SetMarkerColor(kRed+1);
  hPTClone[1]->SetLineColor(kRed-4);
  hPTClone[1]->SetMarkerColor(kRed-4);
  hPTClone[2]->SetLineColor(kOrange+1);
  hPTClone[2]->SetMarkerColor(kOrange+1);
  hPTClone[3]->SetLineColor(kOrange);
  hPTClone[3]->SetMarkerColor(kOrange);
  hPTClone[4]->SetLineColor(kSpring+8);
  hPTClone[4]->SetMarkerColor(kSpring+8);
  hPTClone[5]->SetLineColor(kAzure-4);
  hPTClone[5]->SetMarkerColor(kAzure-4);
  hPTClone[6]->SetLineColor(kBlue);
  hPTClone[6]->SetMarkerColor(kBlue);
  hPTClone[7]->SetLineColor(kBlue+3);
  hPTClone[7]->SetMarkerColor(kBlue+3);
  hPTClone[8]->SetLineColor(kBlack);
  hPTClone[8]->SetMarkerColor(kBlack);
  
  ffit[0]->SetLineColor(kRed+1);
  ffit[0]->SetLineStyle(9);
  ffit[1]->SetLineColor(kRed-4);
  ffit[1]->SetLineStyle(9);
  ffit[2]->SetLineColor(kOrange+1);
  ffit[2]->SetLineStyle(9);
  ffit[3]->SetLineColor(kOrange);
  ffit[3]->SetLineStyle(9);
  ffit[4]->SetLineColor(kSpring+8);
  ffit[4]->SetLineStyle(9);
  ffit[5]->SetLineColor(kAzure-4);
  ffit[5]->SetLineStyle(9);
  ffit[6]->SetLineColor(kBlue);
  ffit[6]->SetLineStyle(9);
  ffit[7]->SetLineColor(kBlue+3);
  ffit[7]->SetLineStyle(9);
  ffit[8]->SetLineColor(kBlack);
  ffit[8]->SetLineStyle(9);

  hPTClone[0]->SetTitle("0.0-1.0");
  hPTClone[1]->SetTitle("1.0-5.0");
  hPTClone[2]->SetTitle("5.0-10.0");
  hPTClone[3]->SetTitle("10.0-15.0");
  hPTClone[4]->SetTitle("15.0-30.0");
  hPTClone[5]->SetTitle("30.0-50.0");
  hPTClone[6]->SetTitle("50.0-70.0");
  hPTClone[7]->SetTitle("70.0-100.0");
  hPTClone[8]->SetTitle("0.0-100.0");
  hPTClone[8]->SetMarkerStyle(25);

  TCanvas* c2 = new  TCanvas();
  c2->SetLogy();
  c2->SetGridx();
  c2->SetGridy();
  TLegend* l = new TLegend(0.1,0.7,0.48,0.9);
  l->AddEntry(hPTClone[0],"","LEP");
  l->SetHeader("V0 selection (percentiles):");
  hPTClone[0]->SetAxisRange(0.0,10.,"X");
  hPTClone[0]->SetAxisRange(1E-7,0.1,"Y");
  hPTClone[0]->Draw();
  ffit[0]->Draw("SAME");
  for(int j=1;j < (ctr-1);j++){
    l->AddEntry(hPTClone[j],"","LEP");
    hPTClone[j]->Draw("SAME");
    ffit[j]->Draw("SAME");   
    }
  l->AddEntry(ffit[ctr-1],"","L");
  l->Draw("SAME");

  new TCanvas;
  NormYieldsNormNch->SetName("NormYieldsNormNch");
  NormYieldsNormNch->SetTitle("Normalized Integrated Yields vs #frac{< n_{ch} >}{(< n_{ch} >)_{MB}} ");
  NormYieldsNormNch->GetXaxis()->SetTitle("#frac{< n_{ch} >}{(< n_{ch} >)_{MB}}");
  NormYieldsNormNch->GetYaxis()->SetTitle("#frac{1/< n_{ch} > dN/dy}{(1/< n_{ch} > dN/dy)_{MB}}");
  NormYieldsNormNch->Draw();

  new TCanvas;
  NormYieldsNch->SetName("NormYieldsNch");
  NormYieldsNch->SetTitle("Normalized Integrated Yields vs < n_{ch}> ");
  NormYieldsNch->GetXaxis()->SetTitle("< n_{ch} >");
  NormYieldsNch->GetYaxis()->SetTitle("#frac{1/< n_{ch} > dN/dy}{(1/< n_{ch} > dN/dy)_{MB}}");
  NormYieldsNch->Draw();

  new TCanvas;
  YieldsNch->SetName("YieldsNch");
  YieldsNch->SetTitle("Integrated Yields vs < n_{ch}> ");
  YieldsNch->GetXaxis()->SetTitle("< n_{ch} >");
  YieldsNch->GetYaxis()->SetTitle("dN/dy");
  YieldsNch->Draw();

  TFile* Write = new TFile (outputfile, "recreate");
  YieldsNch->Write();
  NormYieldsNch->Write();
  NormYieldsNormNch->Write();
  for(int j=0;j < (ctr);j++){
    hPTClone[j]->Write();
    ffit[j]->Write();
  }

  Write->Close();
  
 }


Double_t LevyTsallis_Func(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */
  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t n = p[1];
  Double_t C = p[2];
  Double_t norm = p[3];

  Double_t part1 = (n - 1.) * (n - 2.);
  Double_t part2 = n * C * (n * C + mass * (n - 2.));
  Double_t part3 = part1 / part2;
  Double_t part4 = 1. + (mt - mass) / n / C;
  Double_t part5 = TMath::Power(part4, -n);
  return pt * norm * part3 * part5;
}

TF1 * LevyTsallis(const Char_t *name, Double_t mass, Double_t n = 5., Double_t C = 0.1, Double_t norm = 1.)
{
  
  TF1 *fLevyTsallis = new TF1(name, LevyTsallis_Func, 0., 10., 4);
  fLevyTsallis->SetParameters(mass, n, C, norm);
  fLevyTsallis->SetParNames("mass", "n", "C", "norm");
  fLevyTsallis->FixParameter(0, mass);
  fLevyTsallis->SetParLimits(1, 1.e-3, 1.e3);
  fLevyTsallis->SetParLimits(2, 1.e-3, 1.e3);
  fLevyTsallis->SetParLimits(3, 1.e-5, 1.e5);
  return fLevyTsallis;
}
  
//------------------------------------------------------
  

