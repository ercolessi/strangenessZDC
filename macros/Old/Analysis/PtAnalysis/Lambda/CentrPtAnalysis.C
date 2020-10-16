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


void CentrPtAnalysis() {

  const char *outputname = "yielsPt.root";

  const int ctr = 11;
  int cmin [ctr]= {0,10,20,30,40,50,60,70,80,90,0};
  int cmax [ctr]= {10,20,30,40,50,60,70,80,90,100,100};
  double mult[ctr],n[ctr];
  
  //Efficiencies 
  TFile *e = new TFile("efficienzeLambda.root");
  
  TH1D* hEffpt[ctr];
  TH2F *hMinvPt[ctr]; 
  TH2F *hMinvPt2[ctr];
  TH1D* hPT[ctr];
  TH1D* hBinWidth[ctr];

  for (int j = 0 ; j < ctr ; j++) {
    
    hEffpt[j] = (TH1D*) e->Get(Form("hEffptAll%i",j)); 
    hEffpt[j]->Sumw2();
    hMinvPt[j]= new TH2F(Form("MinvPt%i",j)," ",55,0.1,10.,50,1.09,1.14);
    SetBinX(hMinvPt[j]);
    hMinvPt2[j]= new TH2F(Form("MinvPt2%i",j)," ",55,0.1,10.,50,1.09,1.14);
    SetBinX(hMinvPt2[j]);
    
  }
   
  //Fit function for invariant mass
  TF1 *ffit = new TF1("ffit","gaus(0)+pol2(3)",1.09,1.14);
  ffit->SetParameter(0,0.0002);
  ffit->SetParameter(1,1.1159);
  ffit->SetParLimits(1,1.1155,1.1165);
  ffit->SetParameter(2,1.393e-03);
  ffit->SetParLimits(2,1.39E-3,1.4E-3);
  ffit->SetNpx(1000);

  //Open data files
  TFile *f = new TFile("/alistorage/fercolessi/TreeCreation/StrangenessAnalysis13TeV/Strangeness15h.root");
  f->ls();
  TTree *te = (TTree *) f->Get("PWGLF_StrVsMult/fTreeEvent");
  te->ls();
  TTree *tc = (TTree *) f->Get("PWGLF_StrVsMult/fTreeV0");
  tc->ls();
  
  for (int j = 0 ; j < ctr ; j++) {

    //Multiplicity and number of events
    for(int i=0;i < te->GetEntries();i++){
      te->GetEvent(i);
          
      if(te->GetLeaf("fCentrality")->GetValue() < cmin[j]) continue;
      if(te->GetLeaf("fCentrality")->GetValue() > cmax[j]) continue;
      
      n[j]++;
      mult[j] += te->GetLeaf("fNTracksGlobal2015")->GetValue(); 
    }
    

    //Loop over V0 Tree entries 
    for(int i=0;i < tc->GetEntries();i++) {
      
      tc->GetEvent(i);  

      if(tc->GetLeaf("fTreeVariableCentrality")->GetValue() < cmin[j]) continue;
      if(tc->GetLeaf("fTreeVariableCentrality")->GetValue() > cmax[j]) continue;
      
      if (tc->GetLeaf("fTreeVariableRapLambda")->GetValue() > 0.5) continue;
      if (tc->GetLeaf("fTreeVariableRapLambda")->GetValue() < -0.5) continue;
      
      double pt = tc->GetLeaf("fTreeVariablePt")->GetValue();
      double Lambdas = (tc->GetLeaf("fTreeVariableInvMassLambda")->GetValue());
      
      hMinvPt[j]->Fill(pt,Lambdas);
      hMinvPt2[j]->Fill(pt,Lambdas);
    }
    
    
    //TCanvas* c = new TCanvas ();
    //c->SetLogz();
    //hMinvPt[j]->SetXTitle("p_{T} (GeV/c)");
    //hMinvPt[j]->SetYTitle("M_{inv} (GeV/c^2)");  
    // hMinvPt[j]->Draw("colz");
    
    hPT[j] = hMinvPt[j]->ProjectionX(); //su x hai pt così
    hPT[j]->Reset();
    
    hBinWidth[j] = hMinvPt2[j]->ProjectionX(); //su x hai pt così
    hBinWidth[j]->Reset();
    
    const int nbins = hPT[j]->GetNbinsX();
    float yield[nbins], eyield[nbins];
    
    new TCanvas;
    //loop over bins in PT
    for (int i=0; i<= nbins; i++) {
      
      TH1D* hProjMinv = hMinvPt[j]->ProjectionY("hProjMinv",i,i);
      double binwidthMinv = hProjMinv->GetBinWidth(i);

      hProjMinv->Fit(ffit,"","",1.1,1.13);
      ffit->SetParameter(3,0);
      ffit->SetParameter(4,0);
      ffit->SetParameter(5,0);
     
      yield[i] = (ffit->Integral(1.09,1.14))/binwidthMinv/(hPT[j]->GetBinWidth(i));
      eyield[i] = TMath::Sqrt(TMath::Abs(yield[i]));

      //Fill PT Spectra
      hPT[j]->SetBinContent(i,yield[i]);
      hPT[j]->SetBinError(i,eyield[i]);
      
      hBinWidth[j]->SetBinContent(i,hBinWidth[j]->GetBinWidth(i));

  }

  //new TCanvas;
  hPT[j]->SetLineColor(j+1);
  hPT[j]->SetMarkerStyle(20);
  hPT[j]->SetMarkerSize(1.);
  if (j!=9) hPT[j]->SetMarkerColor(j+1);
  if (j==9) hPT[j]->SetMarkerColor(46);
  hPT[j]->SetYTitle("yields");
 
  new TCanvas;
  hBinWidth[j]->SetYTitle("Bin Width");
  //hBinWidth[j]->Draw();

  hPT[j]->Divide(hPT[j],hEffpt[j],1,1);
  hPT[j]->SetXTitle("p_{T} (GeV/c)");
  hPT[j]->SetName(Form("hPT_%i-%i",cmin[j],cmax[j]));
  hPT[j]->SetTitle(Form("%i % < V0 percentile < %i %",cmin[j],cmax[j]));
  
  }
  
  new TCanvas;
  hPT[0]->Draw();
  for(int j=1;j < (ctr-1);j++){
    hPT[j]->Draw("SAME");
  }

  TFile* Write = new TFile (outputname, "recreate");
    for (int j = 0 ; j < ctr ; j++) {
      hPT[j]->Write();
    }
   Write->Close();
   
  
  for(int j=0;j < ctr;j++){
    printf("%i-%i) n = %i -- mult = %f\n",cmin[j],cmax[j],n[j],mult[j]);
  }
  
  
}


//---------------------------------------------------------------------------------------------
void BinLogX(TH2 *h) {

  // Method for logarithmic binning of histograms

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
