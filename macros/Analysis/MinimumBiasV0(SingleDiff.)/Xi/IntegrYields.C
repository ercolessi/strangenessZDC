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

void IntegrYields() {
  
  const char *outputname = "A_APXiLevyTsallis.root";
  const char *ParticleSpectraFile = "PtSpectraXi.root";
  const char *AntiParticleSpectraFile = "../AntiXi/PtSpectraAntiXi.root";

  //General var
  Double_t mass = 1.32131;
  const int ctr = 9;
  Double_t cmin [ctr]= {0.0, 1.0, 5.0, 10.0, 15.0 , 30.0, 50.0 , 70.0 ,0.0};
  Double_t cmax [ctr]= {1.0, 5.0, 10.0 , 15.0, 30.0 , 50.0 , 70.0 , 100.0 , 100.0};

  //Particle variables: P ----------------------------------------------------------------------------------------------
  double mult_P[ctr];
  int nevt_P[ctr];  
  TH1D* hPT_P[ctr];
  TH1D* hPTClone_P[ctr];
  double n_P[ctr];
  double C_P[ctr];
  double norm_P[ctr];  
  TF1* ffit_P[ctr];
  Double_t yieldsPart[ctr], eyieldsPart[ctr];
  TH1D* hmult_P;
  TH1D* hevt_P;

  //Input Files
  TFile *P = new TFile(ParticleSpectraFile);

  hmult_P = (TH1D*) P->Get("hmult");
  hevt_P = (TH1D*) P->Get("hevt");
  
  for (int j = 0 ; j < (ctr) ; j++) {   
    hPT_P[j] = (TH1D*) P->Get(Form("hPT_%f-%f",cmin[j],cmax[j]));
    hPTClone_P[j]= (TH1D*)hPT_P[j]->Clone("hPTClone");
    hPTClone_P[j]->Reset();
    hPTClone_P[j]->SetTitle(Form("%f - %f",cmin[j],cmax[j]));
    
    mult_P[j]=0;
    nevt_P[j]=0;
    n_P[j]=0;
    C_P[j]=0;
    norm_P[j]=0;
    ffit_P[j] = LevyTsallis("LevyTsallis",mass);
    yieldsPart[j]=0;
    eyieldsPart[j]=0;
  
    for (int i = 7; i <= hPT_P[j]->GetNbinsX(); i++){
      double binc_P = hPT_P[j]->GetBinContent(i);
      double ebinc_P = hPT_P[j]->GetBinError(i);
      hPTClone_P[j]->SetBinContent(i,binc_P);
      hPTClone_P[j]->SetBinError(i,ebinc_P);
      yieldsPart[j] += (hPT_P[j]->GetBinContent(i)*hPT_P[j]->GetBinWidth(i));
      eyieldsPart[j] += (hPT_P[j]->GetBinError(i)*hPT_P[j]->GetBinWidth(i)*hPT_P[j]->GetBinWidth(i)*hPT_P[j]->GetBinError(i));
    }  
    eyieldsPart[j] = TMath::Sqrt(eyieldsPart[j]);
    hPTClone_P[j]->SetMarkerStyle(21);
  }

  //loop over multiplicities
  for (int j = 0 ; j < (ctr) ; j++) {   

    hPT_P[j]->Fit(ffit_P[j],"","",0.6,10.);  
    n_P[j] = ffit_P[j]->GetParameter(1);
    C_P[j] = ffit_P[j]->GetParameter(2);
    norm_P[j] = ffit_P[j]->GetParameter(3);    
    cout << " Chi Square = " << ffit_P[j]->GetChisquare()/ffit_P[j]->GetNDF() << endl;
    cout << " Integral = " << ffit_P[j]->Integral(0.6,10.)/ffit_P[j]->Integral(0.,100.) << endl;     

    yieldsPart[j] += ffit_P[j]->Integral(0.0,0.6) + ffit_P[j]->Integral(10.0,100.0);

    n_P[j] = hevt_P->GetBinContent(j+1);
    mult_P[j] = hmult_P->GetBinContent(j+1);

  }

 
  //AntiParticle variables: AP ----------------------------------------------------------------------------------------------
  double mult_AP[ctr];
  int nevt_AP[ctr];  
  TH1D* hPT_AP[ctr];
  TH1D* hPTClone_AP[ctr];
  double n_AP[ctr];
  double C_AP[ctr];
  double norm_AP[ctr];  
  TF1* ffit_AP[ctr];
  Double_t yieldsAntiPart[ctr], eyieldsAntiPart[ctr];
  TH1D* hmult_AP;
  TH1D* hevt_AP;

  //Input Files
  TFile *AP = new TFile(AntiParticleSpectraFile);

  hmult_AP = (TH1D*) AP->Get("hmult");
  hevt_AP = (TH1D*) AP->Get("hevt");
  
  for (int j = 0 ; j < (ctr) ; j++) {   
    hPT_AP[j] = (TH1D*) AP->Get(Form("hPT_%f-%f",cmin[j],cmax[j]));
    hPTClone_AP[j]= (TH1D*)hPT_AP[j]->Clone("hPTClone");
    hPTClone_AP[j]->Reset();
    hPTClone_AP[j]->SetTitle(Form("%f - %f",cmin[j],cmax[j]));
    
    mult_AP[j]=0;
    nevt_AP[j]=0;
    n_AP[j]=0;
    C_AP[j]=0;
    norm_AP[j]=0;
    ffit_AP[j] = LevyTsallis("LevyTsallis",mass);
    yieldsAntiPart[j]=0;
    eyieldsAntiPart[j]=0;
    
    for (int i = 7; i <= hPT_AP[j]->GetNbinsX(); i++){
      double binc_AP = hPT_AP[j]->GetBinContent(i);
      double ebinc_AP = hPT_AP[j]->GetBinError(i);
      hPTClone_AP[j]->SetBinContent(i,binc_AP);
      hPTClone_AP[j]->SetBinError(i,ebinc_AP);
      yieldsAntiPart[j] += (hPT_AP[j]->GetBinContent(i)*hPT_AP[j]->GetBinWidth(i));
      eyieldsAntiPart[j] += (hPT_AP[j]->GetBinError(i)*hPT_AP[j]->GetBinWidth(i)*hPT_AP[j]->GetBinWidth(i)*hPT_AP[j]->GetBinError(i));
    }  
    eyieldsAntiPart[j] = TMath::Sqrt(eyieldsAntiPart[j]);
    hPTClone_AP[j]->SetMarkerStyle(21);
  }

  //loop over multiplicities
  for (int j = 0 ; j < (ctr) ; j++) {   
    
    hPT_AP[j]->Fit(ffit_AP[j],"","",0.6,10.);  
    n_AP[j] = ffit_AP[j]->GetParameter(1);
    C_AP[j] = ffit_AP[j]->GetParameter(2);
    norm_AP[j] = ffit_AP[j]->GetParameter(3);    
    cout << " Chi Square = " << ffit_AP[j]->GetChisquare()/ffit_AP[j]->GetNDF() << endl;
    cout << " Integral = " << ffit_AP[j]->Integral(0.6,10.)/ffit_AP[j]->Integral(0.,100.) << endl;     

    yieldsAntiPart[j] += ffit_AP[j]->Integral(0.0,0.6) + ffit_AP[j]->Integral(10.0,100.0);

    n_AP[j] = hevt_AP->GetBinContent(j+1);
    mult_AP[j] = hmult_AP->GetBinContent(j+1);

  }

  //Calculate yields 
  Double_t normyields[ctr];
  Double_t enormyields[ctr];
  Double_t multTOT[ctr];
  Double_t yieldsTOT[ctr];
  yieldsTOT[ctr-1]=yieldsAntiPart[ctr-1]+yieldsPart[ctr-1];
  multTOT[ctr-1] = mult_AP[ctr-1]+mult_P[ctr-1];
  normyields[ctr-1] = yieldsTOT[ctr-1]/multTOT[ctr-1];
  enormyields[ctr-1] = TMath::Sqrt((eyieldsAntiPart[ctr-1]*eyieldsAntiPart[ctr-1]))/multTOT[ctr-1];
  
  Double_t  multScaled[ctr-1];
  Double_t emultScaled[ctr-1];
  Double_t prodScaled[ctr-1];
  Double_t eprodScaled[ctr-1];
  Double_t eyieldsTOT[ctr-1];
     
  for(int j=0;j < (ctr-1);j++){

    multTOT[j] = mult_AP[j]+mult_P[j];
    yieldsTOT[j]=yieldsAntiPart[j]+yieldsPart[j];
    normyields[j] = yieldsTOT[j]/multTOT[j];
    eyieldsTOT[j]=TMath::Sqrt((eyieldsAntiPart[j]*eyieldsAntiPart[j])+(eyieldsPart[j]*eyieldsPart[j]));
    enormyields[j] = eyieldsTOT[j]*1./multTOT[j];
    
    multScaled[j] = multTOT[j]*1./multTOT[ctr-1];
    emultScaled[j] = 0;
    prodScaled[j] = normyields[j]/normyields[ctr-1];
    eprodScaled[j] = enormyields[j]/normyields[ctr-1];
  }

  const int npoint = ctr-1;
  TGraphErrors *NormYieldsNormNch = new TGraphErrors(npoint,multScaled,prodScaled,emultScaled,eprodScaled);
  TGraphErrors *NormYieldsNch = new TGraphErrors(npoint,multTOT,prodScaled,emultScaled,eprodScaled);
  TGraphErrors *YieldsNch = new TGraphErrors(npoint,multTOT,yieldsTOT,emultScaled,eyieldsTOT);
  
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

  TFile* Write = new TFile (outputname, "recreate");
  YieldsNch->Write();
  NormYieldsNch->Write();
  NormYieldsNormNch->Write();
  
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
  fLevyTsallis->SetParLimits(3, 1.e-7, 1.e7);
  return fLevyTsallis;
}
  
//------------------------------------------------------
  
