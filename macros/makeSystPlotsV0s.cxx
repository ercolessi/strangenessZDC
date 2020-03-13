#include <iostream>
#include <fstream>
//#define max(a,b) (a>b ? a : b)
#include <TString.h>
#include <TFile.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TClonesArray.h>
#include <TLatex.h>
#include <TMath.h>
#include <TF1.h>
#include <TLeaf.h>
#include <TH3.h>

const int nfiles = 37;
double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr );

void makeCutAnalysis(TString lV0Type = "Lambda", TString lWhichEstimator = "V0M", Double_t lMultBoundLo = 0.0, Double_t lMultBoundHi = 100.0, TString lWhichSystVar = "V0Radius");

//Main Function
void makeSystPlotsV0s(TString Type = "Lambda", Double_t Low = 0.0, Double_t High = 100.0){
  
  makeCutAnalysis(Type,"V0M",Low,High,"V0Radius");
  //makeCutAnalysis(Type,"V0M",Low,High,"V0CosPA");
  /*makeSystPlotsV0s(Type,"V0M",Low,High,"DCANegToPV");
  makeSystPlotsV0s(Type,"V0M",Low,High,"DCAPosToPV");
  makeSystPlotsV0s(Type,"V0M",Low,High,"DCAV0Daughters");*/

  return;
}

//---------------------------------------------------------------
void makeCutAnalysis(TString lV0Type = "Lambda",
		      TString lWhichEstimator = "V0M",
		      Double_t lMultBoundLo = 0.0,
          Double_t lMultBoundHi = 100.0,
		      TString lWhichSystVar = "V0Radius"){

  Double_t loosestcut = 0.0;
  Double_t veryloosecut = 0.0;
  Double_t loosecut = 0.0;
  Double_t tightcut = 0.0;
  Double_t verytightcut = 0.0;
  Double_t minrange = 0.0;
  Double_t maxrange = 0.0;
  if(lWhichSystVar == "V0Radius") {
    loosestcut = 0.300; 
    veryloosecut = 0.300;
    loosecut = 0.4;
    tightcut = 0.6;
    verytightcut = 0.7;
    minrange = 0.;
    maxrange = 1.;
  }
  if(lWhichSystVar == "V0CosPA") {
    loosestcut = 0.95; 
    veryloosecut = 0.95;
    loosecut = 0.96;
    tightcut = 0.98;
    verytightcut = 0.99;
  }
  if(lWhichSystVar == "DCAPosToPV" || lWhichSystVar == "DCANegToPV") {
    loosestcut = 0.050; 
    veryloosecut = 0.050;
    loosecut = 0.055;
    tightcut = 0.070;
    verytightcut = 0.080;
  }
  if(lWhichSystVar == "DCAV0Daughters") {
    loosestcut = 1.5; 
    veryloosecut = 1.5;
    loosecut = 1.25;
    tightcut = 0.75;
    verytightcut = 0.50;
  }
  if(lWhichSystVar == "ProperLifetime") {
    loosestcut = 40.; 
    veryloosecut = 40.;
    loosecut = 30.;
    tightcut = 12.;
  }
   if(lWhichSystVar == "SigmaForSignalExtraction") {
    loosestcut = 7.; 
    loosecut = 7.;
    tightcut = 5.;
    verytightcut = 4.;
  }
   if(lWhichSystVar == "NumberOfCrossedRows") {
    loosestcut = 75.; 
    tightcut = 75.;
    verytightcut = 80.;
  }
  
  cout<<"  "<<endl;
  cout<<"-------------------------------------------------------------"<<endl;
  cout<<"        Macro to make cut study of systematics for V0s       "<<endl;
  cout<<"-------------------------------------------------------------"<<endl;
  cout<<endl;

  //Set data files
  TString lSystFile = "Results-Systematics-V0M_000_100_ZDC_000_100_UseMCRatioFD"; 
  if(lV0Type == "Lambda"    ) lSystFile.Append("-Lambda-");
  if(lV0Type == "AntiLambda") lSystFile.Append("-AntiLambda-");

  TString lDataFilename[nfiles];
  for (int i = 0; i < nfiles; i++){ //[0] is loosest cut
    lDataFilename[i] = lSystFile + lWhichSystVar + Form("-%i.root",i+1);
  }

  TFile* InputFile[nfiles]; 
  for (int i = 0; i < nfiles; i++){
    InputFile[i] = new TFile(lDataFilename[i].Data());
  }
 
   
  //Get real raw histograms from files
  TH1D* lHistPtRaw[nfiles];
  for (int i = 0; i < nfiles; i++){
    lHistPtRaw[i] = (TH1D *) InputFile[i]->Get("lInvMassReal/lInvMassRealRawData/fHistPtRaw");
  }

  //Get MC raw histograms from files
  TH1D* lHistPtMC[nfiles];
  for (int i = 0; i < nfiles; i++){
    lHistPtMC[i] = (TH1D *) InputFile[i]->Get("lInvMassMC/lInvMassMCRawData/fHistEffNumerator");
  }

  //Fill Sgn Loss Fraction
  Double_t binlimits[nfiles];
  for (int i = 0; i < nfiles; i++){
  	binlimits[i] = (minrange + (maxrange-minrange)/nfiles*i);
    cout <<binlimits[i] << endl;
  }
  Long_t binnumb = sizeof(binlimits)/sizeof(Double_t) - 1;

  TH1D* SgnLossFraction = new TH1D("SgnLossFraction", Form(" ; Fraction of signal loss; %s",lWhichSystVar.Data()), binnumb, binlimits);
  TH1D* SgnLossFractionMC = new TH1D("SgnLossFractionMC", Form(" ; Fraction of signal loss; %s",lWhichSystVar.Data()), binnumb, binlimits);

  const int nbins = lHistPtRaw[0]->GetNbinsX();
  Double_t maxy = 0.0;
  Double_t errori = 0.0;
  Double_t error0 = 0.0;
  Double_t erroriMC = 0.0;
  Double_t error0MC = 0.0;

  Double_t entryi = 0.0;
  Double_t entry0 = 0.0;
  for (int i = 1; i < binnumb; i++){
     	
    Double_t entry = 1. - (double)(lHistPtRaw[i]->Integral(1,nbins))/(double)(lHistPtRaw[0]->Integral(1,nbins));
    //Double_t entryMC = 1. - (double)(lHistPtMC[i]->Integral(1,nbins))/(double)(lHistPtMC[0]->Integral(1,nbins));
   
    for (int j=0; j<lHistPtMC[0]->GetNbinsX();j++){
      entryi += lHistPtMC[i]->GetBinContent(j);
      entry0 += lHistPtMC[11]->GetBinContent(j);
    }
    cout << entryi << "  " << entry0 << endl;
    Double_t entryMC = 1. - entryi/entry0;

    Double_t Inti = (double)lHistPtRaw[i]->IntegralAndError(1,nbins,errori);
    Double_t Int0 = (double)lHistPtRaw[11]->IntegralAndError(1,nbins,error0); 
    Double_t IntiMC = (double)lHistPtMC[i]->IntegralAndError(1,nbins,erroriMC);
    Double_t Int0MC = (double)lHistPtMC[11]->IntegralAndError(1,nbins,error0MC);    
    
    Double_t error = (double) (TMath::Sqrt(errori*errori/(Inti*Inti) + error0*error0/(Int0*Int0))*entry);
    Double_t errorMC = (double) (TMath::Sqrt(erroriMC*erroriMC/(IntiMC*IntiMC) + error0MC*error0MC/(Int0MC*Int0MC))*entryMC);

    maxy = max(entry,maxy);

    if (entry!=0) SgnLossFraction->SetBinContent(i,entry);
    else SgnLossFraction->SetBinContent(i,0.00000001);
    //cout<< i <<"  "<<entry << endl;
    SgnLossFraction->SetBinError(i,error);
    SgnLossFractionMC->SetBinContent(i,entryMC);
    SgnLossFractionMC->SetBinError(i,errorMC);
  }
  

  //Cut Markers
  TLine* fvl = new TLine(veryloosecut,0.,veryloosecut,maxy);
  fvl->SetLineColor(kRed);
  fvl->SetLineWidth(2);
  fvl->SetLineStyle(7);
  TLine* fl = new TLine(loosecut,0.,loosecut,maxy);
  fl->SetLineColor(kBlue+1);
  fl->SetLineWidth(2);
  fl->SetLineStyle(7);
  TLine* ft = new TLine(tightcut,0.,tightcut,maxy);
  ft->SetLineColor(kGreen+3);
  ft->SetLineWidth(2);
  ft->SetLineStyle(7);
  TLine* fvt = new TLine(verytightcut,0.,verytightcut,maxy);
  fvt->SetLineColor(kAzure+9);
  fvt->SetLineWidth(2);
  fvt->SetLineStyle(7);

  //Draw Options
  SgnLossFraction->SetMarkerStyle(20);
  SgnLossFraction->SetMarkerSize(0.8);
  SgnLossFraction->SetMarkerColor(kBlack);
  SgnLossFractionMC->SetMarkerStyle(20);
  SgnLossFractionMC->SetMarkerSize(0.8);
  SgnLossFractionMC->SetMarkerColor(kRed);

  new TCanvas;
  SgnLossFraction->Draw();
  //SgnLossFractionMC->Draw();
  if(veryloosecut!=0) fvl->Draw("SAME");
  if(loosecut!=0) fl->Draw("SAME");
  if(tightcut!=0) ft->Draw("SAME");
  if(verytightcut!=0) fvt->Draw("SAME");
}
  