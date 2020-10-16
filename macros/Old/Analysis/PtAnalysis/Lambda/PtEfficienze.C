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


void SetBinX(TH1 *h);


void PtEfficienze() {
  //Definitions
  int pdg = 3122;
  float cmin = 0;
  float cmax = 100;
  ;
  int icmin = cmin+1;
  int icmax = cmax;
  double const ptmax = 6;
  
  //Open file
  TFile *f = new TFile("MonteCarloMerged15f.root");

  //Get Branches
  TTree *te = (TTree *) f->Get("PWGLF_StrVsMult/fTreeEvent");
  TTree *t = (TTree *) f->Get("PWGLF_StrVsMult_MC/fTreeV0");
  TList *l = (TList *) f->Get("PWGLF_StrVsMult_MC/cList");
  TH3D *hGen = (TH3D *) l->FindObject("fHistGeneratedPtVsYVsCentralityLambda");
  TH3D *hReco = new TH3D("hReco"," ",55,0.,10.,20,-0.5,0.5,100,0.,100.);
  SetBinX(hReco);
  int n=t->GetEntries();

  float ptMC,yMC;
  float centrality;

  //Loop over events
  for(int i=0;i < n;i++){

    t->GetEvent(i);

    if(t->GetLeaf("fTreeVariablePID")->GetValue() != pdg) continue;

    if (t->GetLeaf("fTreeVariableRapMC")->GetValue() > 0.5) continue;
    if (t->GetLeaf("fTreeVariableRapMC")->GetValue() < -0.5) continue;

    ptMC = t->GetLeaf("fTreeVariablePtMC")->GetValue();
    yMC = t->GetLeaf("fTreeVariableRapMC")->GetValue();
    centrality = t->GetLeaf("fTreeVariableCentrality")->GetValue();

    hReco->Fill(ptMC,yMC,centrality);
  }

  TH1F *hPt = new TH1F("hPt"," ",55,0.,10.);
  SetBinX(hPt);
 
  TH1D *hGenptAll = hGen->ProjectionX("hGenptAll",6,15);
  TH1D *hRecoptAll = hReco->ProjectionX();
  
  hRecoptAll->SetName("hRecoptAll");

  hGenptAll->Sumw2();
  //hRecoptAll->Sumw2();

  for (int i = 1; i <= hPt->GetNbinsX();i++) {
    
    double min = hPt->GetBinCenter(i)-(hPt->GetBinWidth(i)/2);
    double max = min + hPt->GetBinWidth(i);
    
    for (int k = 1; k <= hGenptAll->GetNbinsX();k++) {
      
      double entries = hGenptAll->GetBinContent(k);
      if (hGenptAll->GetBinCenter(k) < max && hGenptAll->GetBinCenter(k) >= min) 
	hPt->AddBinContent(i,entries);
    }
  } 
  
  TH1D *hEffptAll = hReco->ProjectionX();
  SetBinX(hEffptAll);
  hEffptAll->SetName("hEffptAll");
  //hEffptAll->Reset();

  new TCanvas;
  hPt->SetLineColor(4);
  hPt->Draw();

  hGenptAll->SetLineColor(1);
  hGenptAll->Draw("SAME");

  new TCanvas;
  hRecoptAll->SetLineColor(2);
  hRecoptAll->Draw();

  new TCanvas;   
  hEffptAll->Divide(hEffptAll,hPt,1,1,"B");
  hEffptAll->SetTitle("#varepsilon for 0 % < centr < 100 %");
  hEffptAll->SetLineColor(2);
  hEffptAll->SetMarkerStyle(21);
  hEffptAll->SetMarkerColor(2);
  hEffptAll->SetMarkerSize(0.8);
  hEffptAll->Draw("SAME");
  
  //output file
  TFile* Write = new TFile ("efficienzeLambda.root", "recreate");
  hEffptAll->Write();
  
  Write->Close();
  
}



//-----------------------------------------------------------------------------------------------
void SetBinX(TH1 *h) {

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
