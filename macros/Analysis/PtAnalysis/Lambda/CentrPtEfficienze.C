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


void CentrPtEfficienze() {
  //Definitions
  int pdg = 3122;
 
  //Open file
  TFile *f = new TFile("/alistorage/fercolessi/MCAnalysis/MonteCarloMerged15f.root");

  const int ctr = 11;
  int cmin [ctr]= {0,10,20,30,40,50,60,70,80,90,0};
  int cmax [ctr]= {10,20,30,40,50,60,70,80,90,100,100};

  //Get Branches
  TTree *te = (TTree *) f->Get("PWGLF_StrVsMult/fTreeEvent");
  TTree *t = (TTree *) f->Get("PWGLF_StrVsMult_MC/fTreeV0");
  TList *l = (TList *) f->Get("PWGLF_StrVsMult_MC/cList");
  //V0 entries
  int nV0=t->GetEntries();

  TH3D *hGen[ctr];
  TH3D *hReco[ctr];
  TH1F *hPt[ctr];
  TH1D *hGenptAll[ctr];
  TH1D *hRecoptAll[ctr];
  TH1D* hEffptAll[ctr];

  for (int j = 0 ; j < ctr ; j++) {

   hGen[j] = (TH3D *) l->FindObject("fHistGeneratedPtVsYVsCentralityLambda");
   hReco[j] = new TH3D(Form("hReco%i",j)," ",55,0.,10.,20,-0.5,0.5,100,cmin[j],cmax[j]);
   SetBinX(hReco[j]);
   hPt[j] = new TH1F(Form("hPt%i",j)," ",55,0.,10.);
   SetBinX(hPt[j]);
   
  }

   for (int j = 0 ; j < ctr ; j++) {
  
   float ptMC,yMC;
   float centrality;   
   
   //Loop over events
   for(int i=0;i < nV0;i++){
     
     t->GetEvent(i);
     
     if(t->GetLeaf("fTreeVariablePID")->GetValue() != pdg) continue;
     
     if (t->GetLeaf("fTreeVariableRapMC")->GetValue() > 0.5) continue;
     if (t->GetLeaf("fTreeVariableRapMC")->GetValue() < -0.5) continue;
     
     ptMC = t->GetLeaf("fTreeVariablePtMC")->GetValue();
     yMC = t->GetLeaf("fTreeVariableRapMC")->GetValue();
     centrality = t->GetLeaf("fTreeVariableCentrality")->GetValue();
       
     if (centrality > cmin[j] && centrality <= cmax[j] ){
       hReco[j]->Fill(ptMC,yMC,centrality);
     }
     
   }  
   
   hGenptAll[j] = hGen[j]->ProjectionX(Form("hGenptAll%i",j),6,15,cmin[j],cmax[j]);
   hRecoptAll[j] = hReco[j]->ProjectionX();
   
   hRecoptAll[j]->SetName(Form("hRecoptAll%i",j));
   hGenptAll[j]->Sumw2();
   
   for (int i = 1; i <= hPt[j]->GetNbinsX();i++) {
     
     double min = hPt[j]->GetBinCenter(i)-(hPt[j]->GetBinWidth(i)/2);
     double max = min + hPt[j]->GetBinWidth(i);
     
     for (int k = 1; k <= hGenptAll[j]->GetNbinsX();k++) {
       
       double entries = hGenptAll[j]->GetBinContent(k);
       if (hGenptAll[j]->GetBinCenter(k) < max && hGenptAll[j]->GetBinCenter(k) >= min) 
	 hPt[j]->AddBinContent(i,entries);
     }
   } 
   
   hEffptAll[j] = hReco[j]->ProjectionX();
   SetBinX(hEffptAll[j]);
   hEffptAll[j]->SetName(Form("hEffptAll%i",j));
  
   new TCanvas;   
   hEffptAll[j]->Divide(hEffptAll[j],hPt[j],1,1,"B");
   hEffptAll[j]->SetTitle(Form("#varepsilon for %i % < centr < %i %",cmin[j],cmax[j]));
   hEffptAll[j]->SetLineColor(1);
   hEffptAll[j]->SetMarkerStyle(21);
   hEffptAll[j]->SetMarkerColor(1);
   hEffptAll[j]->SetMarkerSize(0.8);
   //hEffptAll[j]->Draw("SAME");
   
   new TCanvas;
   hPt[j]->SetLineColor(4);
   //hPt[j]->Draw();
   
   hGenptAll[j]->SetLineColor(1);
   //hGenptAll[j]->Draw("SAME");
   
   new TCanvas;
   hRecoptAll[j]->SetLineColor(1);
   //hRecoptAll[j]->Draw();
   
   }
   
   //output file
   TFile* Write = new TFile ("efficienzeLambda.root", "recreate");
    for (int j = 0 ; j < ctr ; j++) {
      hEffptAll[j]->Write();
    }
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
