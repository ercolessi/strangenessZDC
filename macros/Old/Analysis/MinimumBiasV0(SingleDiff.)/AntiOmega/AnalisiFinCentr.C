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

void SetBinX(TH1 *h);

void AnalisiFinCentr() {

  int pdg = -3334;
  const char *outputname = "PtSpectraAntiOmega.root";
  const int ctr = 5;
  const int eff = 4;
  Double_t cmin [ctr]= {0.0, 5.0, 15.0, 50.0 ,0.0};
  Double_t cmax [ctr]= {5.0, 15.0, 50.0, 100.0, 100.0};
  double mult[ctr];
  int n[ctr];  
  TH1D* hEffpt;
  TH2F *hMinvPt[ctr]; 
  TH2F *hMinvPt2[ctr];
  TH1D* hPT[ctr];
  TH1D* hBinWidth[ctr];
  TH1D* MuPT[ctr];
  TH1D* SigPT[ctr];
  TH1D* ChiPT[ctr];
  const int nbins = 9;
  float yield[nbins], eyield[nbins];
  TH1D* hmult = new TH1D ("hmult","",ctr,0.,ctr);
  TH1D* hevt = new TH1D ("hevt","",ctr,0.,ctr);
  
  //Efficiencies 
  TFile *e = new TFile("EfficienciesAntiOmega.root");
  hEffpt = (TH1D*) e->Get("hEffptReco3");


  for (int j = 0 ; j < ctr ; j++) {    
    hMinvPt[j]= new TH2F(Form("MinvPt%i",j)," ",nbins,0.1,10.,100,1.62,1.74);
    SetBinX(hMinvPt[j]);
    hMinvPt2[j]= new TH2F(Form("MinvPt2%i",j)," ",nbins,0.1,10.,100,1.62,1.74);
    SetBinX(hMinvPt2[j]);
    
    n[j]=0;
    mult[j]=0;
  }
  
  //Fit function for invariant mass
  TF1 *ffit = new TF1("ffit","gaus(0)/sqrt(2*TMath::Pi())/[2]+pol2(3)",1.62,1.74);
  ffit->SetParameter(0,0.0002);
  ffit->SetParameter(1,1.67251);
  ffit->SetParLimits(1,1.672,1.673);
  ffit->SetParameter(2,2.34086e-03);
  ffit->SetParLimits(2,2.4E-3,2.3E-3);
  ffit->SetNpx(1000);


  //Open data files
  TFile *f = new TFile("/alistorage/fercolessi/TreeCreation/StrangenessAnalysis13TeV/StrangenessTOT15f15h18i.root");
  f->ls();
  TTree *te = (TTree *) f->Get("PWGLF_StrVsMult/fTreeEvent");
  te->ls();
  TTree *tc = (TTree *) f->Get("PWGLF_StrVsMult/fTreeCascade");
  tc->ls();
  
  //loop over 12 centralities
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

      if(tc->GetLeaf("fTreeCascVarCentrality")->GetValue() < cmin[j]) continue;
      if(tc->GetLeaf("fTreeCascVarCentrality")->GetValue() > cmax[j]) continue;
      
      if (tc->GetLeaf("fTreeCascVarRapOmega")->GetValue() > 0.5) continue;
      if (tc->GetLeaf("fTreeCascVarRapOmega")->GetValue() < -0.5) continue;

      //Cuts
      double CascDecayX = tc->GetLeaf("fTreeCascVarCascadeDecayX")->GetValue();
      double CascDecayY = tc->GetLeaf("fTreeCascVarCascadeDecayY")->GetValue();
      double CascDecay2D = TMath::Sqrt(CascDecayX*CascDecayX + CascDecayY*CascDecayY);
      if (CascDecay2D < 0.5) continue;

      double V0DecayX = tc->GetLeaf("fTreeCascVarV0DecayX")->GetValue();
      double V0DecayY = tc->GetLeaf("fTreeCascVarV0DecayY")->GetValue();
      double V0Decay2D = TMath::Sqrt(V0DecayX*V0DecayX + V0DecayY*V0DecayY);
      if (V0Decay2D < 1.1) continue;

      if (tc->GetLeaf("fTreeCascVarDCABachToPrimVtx")->GetValue() < 0.04) continue;
      if (tc->GetLeaf("fTreeCascVarDCAV0ToPrimVtx")->GetValue() < 0.06) continue;
      if (tc->GetLeaf("fTreeCascVarDCAV0Daughters")->GetValue() >1.5 ) continue;
      if (tc->GetLeaf("fTreeCascVarCascCosPointingAngle")->GetValue() <0.97 ) continue;
      if (tc->GetLeaf("fTreeCascVarV0CosPointingAngle")->GetValue() <0.97 ) continue;
      if (tc->GetLeaf("fTreeCascVarDCABachToBaryon")->GetValue() >1.3 ) continue;
      if (pdg == -3334){
	if (tc->GetLeaf("fTreeCascVarDCAPosToPrimVtx")->GetValue() <0.04 ) continue;
	if (tc->GetLeaf("fTreeCascVarDCANegToPrimVtx")->GetValue() <0.03 ) continue;
      }
      if (pdg==3334){
	if (tc->GetLeaf("fTreeCascVarDCAPosToPrimVtx")->GetValue() <0.03 ) continue;
	if (tc->GetLeaf("fTreeCascVarDCANegToPrimVtx")->GetValue() <0.04 ) continue;
      }
      
      if (tc->GetLeaf("fTreeCascVarLeastNbrClusters")->GetValue() <70 ) continue; 
      
      

      if(pdg==3334){
	if (tc->GetLeaf("fTreeCascVarPosNSigmaProton")->GetValue() > 4) continue;
	if (tc->GetLeaf("fTreeCascVarPosNSigmaProton")->GetValue() < -4) continue
									     ;
	if (tc->GetLeaf("fTreeCascVarNegNSigmaPion")->GetValue() > 4) continue;
	if (tc->GetLeaf("fTreeCascVarNegNSigmaPion")->GetValue() < -4) continue;   }
      
      if(pdg==-3334){
	if (tc->GetLeaf("fTreeCascVarNegNSigmaProton")->GetValue() > 4) continue;
	if (tc->GetLeaf("fTreeCascVarNegNSigmaProton")->GetValue() < -4) continue
									     ;
	if (tc->GetLeaf("fTreeCascVarPosNSigmaPion")->GetValue() > 4) continue;
	if (tc->GetLeaf("fTreeCascVarPosNSigmaPion")->GetValue() < -4) continue;   }
      


      double pt = tc->GetLeaf("fTreeCascVarPt")->GetValue();
      double Xi = (tc->GetLeaf("fTreeCascVarMassAsOmega")->GetValue());
      
      hMinvPt[j]->Fill(pt,Xi);
      hMinvPt2[j]->Fill(pt,Xi);
    }
        
    //TCanvas* c = new TCanvas ();
    //c->SetLogz();
    //hMinvPt[j]->SetXTitle("p_{T} (GeV/c)");
    hMinvPt[j]->SetYTitle("M_{inv} (GeV/c^2)");  
    hMinvPt[j]->Draw("colz");
    
    //Prepare histos
    hPT[j] = hMinvPt[j]->ProjectionX(); //su x hai pt così
    hPT[j]->Reset();
    MuPT[j]= (TH1D*)hPT[j]->Clone("MuPT");
    MuPT[j]->Reset();
    SigPT[j] = (TH1D*)hPT[j]->Clone("SigPT");
    SigPT[j]->Reset();
    ChiPT[j] = (TH1D*)hPT[j]->Clone("ChiPT");
    ChiPT[j]->Reset();
    
    //histo for binning
    hBinWidth[j] = hMinvPt2[j]->ProjectionX(); //su x hai pt così
    hBinWidth[j]->Reset();  
           
    //loop over bins in PT -> Get yields and fit parameters 
    for (int i=0; i<= nbins; i++) {
      
      TH1D* hProjMinv = hMinvPt[j]->ProjectionY("hProjMinv",i,i);
      double binwidthMinv = hProjMinv->GetBinWidth(i);
      
      hProjMinv->Fit(ffit,"","",1.62,1.74);
      ffit->SetParameter(3,0);
      ffit->SetParameter(4,0);
      ffit->SetParameter(5,0);
      double mu=ffit->GetParameter(1);
      double emu=ffit->GetParError(1);
      double sigma=ffit->GetParameter(2);
      double esigma=ffit->GetParError(2);
      double Chi=ffit->GetChisquare();
      double NDF=ffit->GetNDF();
      double ChiR;
      if (NDF!=0) ChiR = Chi/NDF;
      else ChiR=0;

      yield[i] =( (ffit->Integral(1.62,1.74))/binwidthMinv/(hPT[j]->GetBinWidth(i)));
      //eyield[i] = TMath::Sqrt(TMath::Abs(yield[i]));
      eyield[i] = (ffit->GetParError(0)/ffit->GetParameter(0))*yield[i];
      
      hPT[j]->SetBinContent(i,TMath::Abs(yield[i]));
      hPT[j]->SetBinError(i,eyield[i]);

      MuPT[j]->SetBinContent(i,mu);
      MuPT[j]->SetBinError(i,emu);
      SigPT[j]->SetBinContent(i,sigma);
      SigPT[j]->SetBinError(i,esigma);
      ChiPT[j]->SetBinContent(i,ChiR);      
      
      hBinWidth[j]->SetBinContent(i,hBinWidth[j]->GetBinWidth(i));
    }
  }

  for(int j=0;j < ctr;j++){
    mult[j]/=n[j];
    printf("%f-%f) n = %i -- mult = %f\n",cmin[j],cmax[j],n[j],mult[j]);
    hmult->SetBinContent(j+1,mult[j]);
    hevt->SetBinContent(j+1,n[j]);
  }
    
  
  //Fill PT Spectra
  for(int j=0;j < ctr;j++){
           
    hPT[j]->SetMarkerStyle(20);
    hPT[j]->SetMarkerSize(1.);
    hPT[j]->SetYTitle("yields");
    
        
    hPT[j]->SetXTitle("p_{T} (GeV/c)");
    hPT[j]->SetYTitle("1/N_{evt} #frac{dN}{dp_{T} dy} ");
    hPT[j]->SetName(Form("hPT_%f-%f",cmin[j],cmax[j]));
    hPT[j]->SetTitle(Form("%f % < V0 percentile < %f %",cmin[j],cmax[j]));
   
  }
  
 hPT[0]->Divide(hPT[0],hEffpt,1,1);
  hPT[0]->Scale(1./n[0]);
  hPT[1]->Divide(hPT[1],hEffpt,1,1);
  hPT[1]->Scale(1./n[1]);

  hPT[2]->Divide(hPT[2],hEffpt,1,1);
  hPT[2]->Scale(1./n[2]);
  hPT[3]->Divide(hPT[3],hEffpt,1,1);
  hPT[3]->Scale(1./n[3]);
 
  hPT[4]->Divide(hPT[4],hEffpt,1,1);
  hPT[4]->Scale(1./n[4]);
   

  hPT[0]->SetLineColor(kRed+1);
  hPT[0]->SetMarkerColor(kRed+1);
  hPT[1]->SetLineColor(kSpring+8);
  hPT[1]->SetMarkerColor(kSpring+8);
  hPT[2]->SetLineColor(kAzure-4);
  hPT[2]->SetMarkerColor(kAzure-4);
  hPT[3]->SetLineColor(kBlue+3);
  hPT[3]->SetMarkerColor(kBlue+3);
  hPT[4]->SetLineColor(kBlack);
  hPT[4]->SetMarkerColor(kBlack);
  
  
  TCanvas* c2 = new  TCanvas();
  c2->SetLogy();
  c2->SetGridx();
  c2->SetGridy();
  c2->BuildLegend(0.5,0.70,0.9,0.9,"","LEP");
  hPT[0]->GetXaxis()->SetRangeUser(0.5,10.);
  hPT[0]->Draw();
  for(int j=1;j < (ctr);j++){
    hPT[j]->Draw("SAME");
  }
  
  TFile* Write = new TFile (outputname, "recreate");
  for (int j = 0 ; j < ctr ; j++) {    
    hPT[j]->Write();
    hevt->Write();
    hmult->Write();
  }

  Write->Close();


/*hPT[0]->Divide(hPT[0],hEffpt[0],1,1);
  hPT[0]->Scale(1./n[0]/mult[0]);
  hPT[1]->Divide(hPT[1],hEffpt[0],1,1);
  hPT[1]->Scale(1./n[1]/mult[1]);

  hPT[2]->Divide(hPT[2],hEffpt[0],1,1);
  hPT[2]->Scale(1./n[2]/mult[2]);
  hPT[3]->Divide(hPT[3],hEffpt[0],1,1);
  hPT[3]->Scale(1./n[3]/mult[3]);
 
  hPT[4]->Divide(hPT[4],hEffpt[1],1,1);
  hPT[4]->Scale(1./n[4]/mult[4]);
  hPT[5]->Divide(hPT[5],hEffpt[1],1,1);
  hPT[5]->Scale(1./n[5]/mult[5]);
  hPT[6]->Divide(hPT[6],hEffpt[2],1,1);
  hPT[6]->Scale(1./n[6]/mult[6]);
  hPT[7]->Divide(hPT[7],hEffpt[2],1,1);
  hPT[7]->Scale(1./n[7]/mult[7]);
  hPT[8]->Divide(hPT[8],hEffpt[3],1,1);
  hPT[8]->Scale(1./n[8]/mult[8]);*/ 
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

  for (int i = 1; i < 3; i++)  {new_bins[i] = 0.5 + new_bins[i-1];}
  for (int i = 3; i < 5; i++) {new_bins[i] = 0.6 + new_bins[i-1];}
  for (int i = 5; i < 6; i++) {new_bins[i] = 0.8 + new_bins[i-1];}
  for (int i = 6; i < 8; i++) {new_bins[i] = 1.0 + new_bins[i-1];}
  for (int i = 8; i < 9; i++) {new_bins[i] = 2.0 + new_bins[i-1];}
  for (int i = 9; i < 10; i++) {new_bins[i] = 3.0 + new_bins[i-1];}
  
  axis->Set(bins, new_bins);
  delete [] new_bins;
}
