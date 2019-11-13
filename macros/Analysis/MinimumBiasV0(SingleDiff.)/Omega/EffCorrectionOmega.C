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

void EffCorrectionOmega() {
  
  //Open file
  TFile *f = new TFile("/alistorage/fercolessi/MCAnalysis/MonteCarloMerged15f.root");

  int pdg = 3334;
  const int ctr = 4;
  const int nbins = 9;
  Double_t cmin [ctr]= {0.0, 15.0 , 50.0 ,0.0};
  Double_t cmax [ctr]= {15.0, 50.0 , 100.0 , 100.0};

  // const int ctr=1;
  TH2F *hMinvPt[ctr]; 
  TH3D* hGen[ctr];
  TH3D *hReco[ctr];
  TH1D* hEffReco[ctr];
  TH1D* hEffMinv[ctr];
  TH1D* hPtMinv[ctr];
  TH1D* hPtReco[ctr];
  TH1F* hPtGen[ctr];
  TH1D* hGenProjX[ctr];
  float yield[nbins], eyield[nbins];

  //TH1D* hEffptAll[ctr];
 
  //Fit function for invariant mass
  TF1 *ffit = new TF1("ffit","gaus(0)/sqrt(2*TMath::Pi())/[2]+pol2(3)",1.62,1.74);
  ffit->SetParameter(0,0.0002);
  ffit->SetParameter(1,1.67251);
  ffit->SetParLimits(1,1.672,1.673);
  ffit->SetParameter(2,2.34086e-03);
  ffit->SetParLimits(2,2.4E-3,2.3E-3);
  ffit->SetNpx(1000);


  //Get Branches
  TTree *te = (TTree *) f->Get("PWGLF_StrVsMult/fTreeEvent");
  TTree *t = (TTree *) f->Get("PWGLF_StrVsMult_MC/fTreeCascade");
  TList *l = (TList *) f->Get("PWGLF_StrVsMult_MC/cList");
  //V0 entries
  int nV0=t->GetEntries();
 
  for (int j = 0 ; j < ctr ; j++) {

    hMinvPt[j]= new TH2F(Form("MinvPt%i",j)," ",nbins,0.001,10.,50,1.62,1.74);
    SetBinX(hMinvPt[j]);

    hGen[j] = (TH3D *) l->FindObject("fHistGeneratedPtVsYVsCentralityOmegaMinus");

    hReco[j] = new TH3D(Form("hReco%i",j)," ", nbins, 0., 10.,20,-0.5,0.5,100,cmin[j],cmax[j]);
    SetBinX(hReco[j]);

    hPtGen[j] = new TH1F(Form("hPtGen%i",j)," ",nbins,0.,10.);
    SetBinX(hPtGen[j]); 

  }

  for (int j = 0 ; j < ctr ; j++) {
  
    double ptMC,yMC;
    double centrality;       
    printf(Form("cmin = %f --- cmax = %f ",cmin[j],cmax[j]));
    
    //Loop over events
    for(int i=0;i < nV0;i++){
      
      t->GetEvent(i);
        
      if(t->GetLeaf("fTreeCascVarPID")->GetValue() != pdg) continue;
     
      if (t->GetLeaf("fTreeCascVarRapMC")->GetValue() > 0.5) continue;
      if (t->GetLeaf("fTreeCascVarRapMC")->GetValue() < -0.5) continue;

      //Cuts
      double CascDecayX = t->GetLeaf("fTreeCascVarCascadeDecayXMC")->GetValue();
      double CascDecayY = t->GetLeaf("fTreeCascVarCascadeDecayYMC")->GetValue();
      double CascDecay2D = TMath::Sqrt(CascDecayX*CascDecayX + CascDecayY*CascDecayY);
      if (CascDecay2D < 0.5) continue;

      double V0DecayX = t->GetLeaf("fTreeCascVarV0DecayXMC")->GetValue();
      double V0DecayY = t->GetLeaf("fTreeCascVarV0DecayYMC")->GetValue();
      double V0Decay2D = TMath::Sqrt(V0DecayX*V0DecayX + V0DecayY*V0DecayY);
      if (V0Decay2D < 1.1) continue;

      if (t->GetLeaf("fTreeCascVarDCABachToPrimVtx")->GetValue() < 0.04) continue;
      if (t->GetLeaf("fTreeCascVarDCAV0ToPrimVtx")->GetValue() < 0.06) continue;
      if (t->GetLeaf("fTreeCascVarDCAV0Daughters")->GetValue() >1.5 ) continue;
      if (t->GetLeaf("fTreeCascVarCascCosPointingAngle")->GetValue() <0.97 ) continue;
      if (t->GetLeaf("fTreeCascVarV0CosPointingAngle")->GetValue() <0.97 ) continue;
      if (t->GetLeaf("fTreeCascVarDCABachToBaryon")->GetValue() >1.3 ) continue;
      if (pdg == -3334){
	if (t->GetLeaf("fTreeCascVarDCAPosToPrimVtx")->GetValue() <0.04 ) continue;
	if (t->GetLeaf("fTreeCascVarDCANegToPrimVtx")->GetValue() <0.03 ) continue;
      }
      if (pdg==3334){
	if (t->GetLeaf("fTreeCascVarDCAPosToPrimVtx")->GetValue() <0.03 ) continue;
	if (t->GetLeaf("fTreeCascVarDCANegToPrimVtx")->GetValue() <0.04 ) continue;
      }
      
      if (t->GetLeaf("fTreeCascVarLeastNbrClusters")->GetValue() <70 ) continue;     

     
      ptMC = t->GetLeaf("fTreeCascVarPtMC")->GetValue();
      yMC = t->GetLeaf("fTreeCascVarRapMC")->GetValue();
      centrality = t->GetLeaf("fTreeCascVarCentrality")->GetValue();
      double XiMC = (t->GetLeaf("fTreeCascVarMassAsOmega")->GetValue());
       
       
      if (centrality > cmin[j] && centrality <= cmax[j] ){      
	hMinvPt[j]->Fill(ptMC,XiMC);    
	hReco[j]->Fill(ptMC,yMC,centrality);
      }
    }  
   
    hGenProjX[j] = hGen[j]->ProjectionX(Form("hGenptProjX%i",j),6,15,cmin[j],cmax[j]);
    hGenProjX[j]->Sumw2();
    hPtReco[j] = hReco[j]->ProjectionX();   
    hPtReco[j]->SetName(Form("hRecoProjX%i",j));
    hPtReco[j]->Sumw2();
   
    //Rebin dell'istogramma generato così da poterlo dividere poi 
    for (int i = 1; i <= hPtGen[j]->GetNbinsX();i++) {
      
      double min = hPtGen[j]->GetBinCenter(i)-(hPtGen[j]->GetBinWidth(i)/2);
      double max = min + hPtGen[j]->GetBinWidth(i);
      
      for (int k = 1; k <= hGenProjX[j]->GetNbinsX();k++) {
	
	double entries = (hGenProjX[j]->GetBinContent(k));
	if (hGenProjX[j]->GetBinCenter(k) < max && hGenProjX[j]->GetBinCenter(k) >= min) 
	  hPtGen[j]->AddBinContent(i,entries);
      }
    } 
   
    hPtMinv[j] = hMinvPt[j]->ProjectionX(); //su x hai pt così
    hPtMinv[j]->Reset();

    //loop over bins in PT -> Get yields and fit parameters 
    for (int i=1; i<= nbins; i++) {      
      TH1D* hProjMinv = hMinvPt[j]->ProjectionY("hProjMinv",i,i);
      double binwidthMinv = hProjMinv->GetBinWidth(i);
            
      hProjMinv->Fit(ffit,"","",1.62,1.74);
      ffit->SetParameter(3,0);
      ffit->SetParameter(4,0);
      ffit->SetParameter(5,0);
      
      yield[i] =(ffit->Integral(1.62,1.74))/binwidthMinv;
      eyield[i] = (ffit->GetParError(0)/ffit->GetParameter(0))*yield[i];
      //eyield[i] = TMath::Sqrt(TMath::Abs(yield[i]));
      
      hPtMinv[j]->SetBinContent(i,TMath::Abs(yield[i]));
      hPtMinv[j]->SetBinError(i,eyield[i]);
    }



    TCanvas* c = new TCanvas ();
    c->SetLogz();
    hMinvPt[j]->SetXTitle("p_{T} (GeV/c)");
    hMinvPt[j]->SetYTitle("M_{inv} (GeV/c^2)");  
    hMinvPt[j]->Draw("colz");
    
    new TCanvas;    
    hEffMinv[j]= (TH1D*)hPtMinv[j]->Clone("hPtMinv");
    hEffReco[j]= (TH1D*)hPtReco[j]->Clone("hPtReco");
    //hEptAll[j]->SetName(Form("hEffptAll%i",j));
    //hEffReco[j]->Draw();
  }

   for (int j = 0 ; j < ctr ; j++) {

     //new TCanvas;   
    hPtMinv[j]->Divide(hPtMinv[j],hPtReco[j],1,1,"B");
    hPtReco[j]->Divide(hPtReco[j],hPtGen[j],1,1,"B");
    
    new TCanvas;
    hPtGen[j]->SetLineColor(4);
   }
   
      
   //output file
   TFile* Write = new TFile ("EfficienciesOmega.root", "recreate");
  for (int j = 0 ; j < ctr ; j++) {
    
    hPtMinv[j]->SetName(Form("hEffptMinvReco%i",j));
    hPtMinv[j]->SetTitle(Form("#varepsilon for %f % < centr < %f %",cmin[j],cmax[j]));
    hPtMinv[j]->SetLineColor(1);
    hPtMinv[j]->SetMarkerStyle(21);
    hPtMinv[j]->SetMarkerColor(1);
    hPtMinv[j]->SetMarkerSize(0.8);
    hPtMinv[j]->GetXaxis()->SetRangeUser(0.5,10.);
    hPtMinv[j]->GetYaxis()->SetRangeUser(0.0,1.1);

    hPtReco[j]->SetName(Form("hEffptReco%i",j));
    hPtReco[j]->SetTitle(Form("#varepsilon for %f % < centr < %f %",cmin[j],cmax[j]));
    hPtReco[j]->SetLineColor(2);
    hPtReco[j]->SetMarkerStyle(21);
    hPtReco[j]->SetMarkerColor(2);
    hPtReco[j]->SetMarkerSize(0.8);
    
    
    hPtMinv[j]->Write();
    hPtReco[j]->Write();
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
  
  for (int i = 1; i < 3; i++)  {new_bins[i] = 0.5 + new_bins[i-1];}
  for (int i = 3; i < 5; i++) {new_bins[i] = 0.6 + new_bins[i-1];}
  for (int i = 5; i < 6; i++) {new_bins[i] = 0.8 + new_bins[i-1];}
  for (int i = 6; i < 8; i++) {new_bins[i] = 1.0 + new_bins[i-1];}
  for (int i = 8; i < 9; i++) {new_bins[i] = 2.0 + new_bins[i-1];}
  for (int i = 9; i < 10; i++) {new_bins[i] = 3.0 + new_bins[i-1];}
  
  
  axis->Set(bins, new_bins);
  delete [] new_bins;
}
