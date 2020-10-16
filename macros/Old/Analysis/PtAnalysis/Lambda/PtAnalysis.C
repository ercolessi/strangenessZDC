#include <TFile.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TClonesArray.h>
#include <TLatex.h>
#include <TMath.h>

// STL includes
#include <iostream>
#include <fstream>
using std::cout;
using std::endl;


void BinLogX(TH2 *h);
void SetBinX(TH2 *h);


void PtAnalysis() {

  const char *outputname = "yielsPt.root";

  //Efficiencies 
  TFile *e = new TFile("efficienzeLambda.root");
  TH1D* hEffpt = (TH1D*) e->Get("hEffptAll"); 
  hEffpt->Sumw2();
  new TCanvas;
  hEffpt->Draw();
   
  Float_t pt;  
 
  //Fit function for invariant mass
  TF1 *ffit = new TF1("ffit","gaus(0)+pol2(3)",1.09,1.14);
  ffit->SetParameter(0,0.0002);
  ffit->SetParameter(1,1.1159);
  ffit->SetParLimits(1,1.1155,1.1165);
  ffit->SetParameter(2,1.393e-03);
  ffit->SetParLimits(2,1.39E-3,1.4E-3);
  ffit->SetNpx(1000);

  int centrV0Min = 0;
  int centrV0Max = 100;

  //Open data files
  TFile *f = new TFile("../TreeCreation/StrangenessAnalysis13TeV/Strangeness15h.root");
  f->ls();
  TTree *te = (TTree *) f->Get("PWGLF_StrVsMult/fTreeEvent");
  te->ls();
  TTree *tc = (TTree *) f->Get("PWGLF_StrVsMult/fTreeV0");
  tc->ls();

  //TH2 definitions
  TH2F *hMinvPt = new TH2F("MinvPt"," ",55,0.1,10.,50,1.09,1.14);
  SetBinX(hMinvPt);
  
  TH2F *hMinvPt2 = new TH2F("MinvPt2"," ",55,0.1,10.,50,1.09,1.14);
  SetBinX(hMinvPt2);

  //Loop over V0 Tree entries 
  for(int i=0;i < tc->GetEntries();i++) {
    
    tc->GetEvent(i);    
    
    if (tc->GetLeaf("fTreeVariableRapLambda")->GetValue() > 0.5) continue;
    if (tc->GetLeaf("fTreeVariableRapLambda")->GetValue() < -0.5) continue;

    pt = tc->GetLeaf("fTreeVariablePt")->GetValue();
    double Lambdas = (tc->GetLeaf("fTreeVariableInvMassLambda")->GetValue());
    
    hMinvPt->Fill(pt,Lambdas);
    hMinvPt2->Fill(pt,Lambdas);
  }
     
  TCanvas* c = new TCanvas ();
  c->SetLogz();
  hMinvPt->SetXTitle("p_{T} (GeV/c)");
  hMinvPt->SetYTitle("M_{inv} (GeV/c^2)");  
  hMinvPt->Draw("colz");

  TH1D* hPT = hMinvPt->ProjectionX(); //su x hai pt così
  hPT->Reset();

  TH1D* hBinWidth = hMinvPt2->ProjectionX(); //su x hai pt così
  hBinWidth->Reset();

  const int nbins = hPT->GetNbinsX();
  float yield[nbins], eyield[nbins];
  
  new TCanvas;
  //loop over bins in PT
  for (int i=0; i<= nbins; i++) {

      TH1D* hProjMinv = hMinvPt->ProjectionY("hProjMinv",i,i);
      double binwidthMinv = hProjMinv->GetBinWidth(i);

      hProjMinv->Fit(ffit,"","",1.1,1.13);
      ffit->SetParameter(3,0);
      ffit->SetParameter(4,0);
      ffit->SetParameter(5,0);
     
      yield[i] = (ffit->Integral(1.09,1.14))/binwidthMinv/(hPT->GetBinWidth(i));
      eyield[i] = TMath::Sqrt(TMath::Abs(yield[i]));

      //Fill PT Spectra
      hPT->SetBinContent(i,yield[i]);
      hPT->SetBinError(i,eyield[i]);
      
      hBinWidth->SetBinContent(i,hBinWidth->GetBinWidth(i));

      if (i==60) {
	hProjMinv->Draw();
	new TCanvas;
      }
  }

  //new TCanvas;
  hPT->SetLineColor(4);
  hPT->SetMarkerStyle(20);
  hPT->SetMarkerSize(0.8);
  hPT->SetMarkerColor(4);
  hPT->SetYTitle("yields");
  //hPT->Draw();

  new TCanvas;
  hBinWidth->SetYTitle("Bin Width");
  hBinWidth->Draw();

  new TCanvas;
  hPT->Divide(hPT,hEffpt,1,1);
  hPT->Draw();

}


//---------------------------------------------------------------------------------------------
void BinLogX(TH2 *h) {

  // Method for the correct logarithmic binning of histograms

  TAxis *axis = h->GetXaxis();
  int bins = axis->GetNbins();

  Double_t from = axis->GetXmin();
  Double_t to = axis->GetXmax();
  Double_t *new_bins = new Double_t[bins + 1];

  new_bins[0] = from;
  Double_t factor = pow(to/from, 1./bins);

  for (int i = 1; i <= bins; i++) {
   new_bins[i] = factor * new_bins[i-1];
  }

  axis->Set(bins, new_bins);
  delete [] new_bins;

}

//-----------------------------------------------------------------------------------------------
void SetBinX(TH2 *h) {

  TAxis *axis = h->GetXaxis();
  int bins = axis->GetNbins();
  
  Double_t from = axis->GetXmin();
  Double_t to = axis->GetXmax();

  Double_t *new_bins = new Double_t[bins + 1];

  new_bins[0] = 0;
  
  for (int i = 1; i < 31; i++) {new_bins[i] = 0.1 + new_bins[i-1];}
  for (int i = 31; i < 41; i++) {new_bins[i] = 0.2 + new_bins[i-1];}
  for (int i = 41; i < 51; i++) {new_bins[i] = 0.3 + new_bins[i-1];}
  for (int i = 51; i < 56; i++) {new_bins[i] = 0.4 + new_bins[i-1];}
  
  axis->Set(bins, new_bins);
  delete [] new_bins;
}
