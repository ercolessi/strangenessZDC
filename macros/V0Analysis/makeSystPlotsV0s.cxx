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

const int nfiles = 41;
double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr );

void makeCutAnalysis(TString lV0Type = "Lambda", TString lWhichEstimator = "V0M", Double_t lMultBoundLo = 0.0, Double_t lMultBoundHi = 100.0, TString lWhichSystVar = "V0Radius");

//Main Function
void makeSystPlotsV0s(TString Type = "Lambda", Double_t Low = 0.0, Double_t High = 100.0){
  
  makeCutAnalysis(Type,"V0M",Low,High,"V0Radius");
 // makeCutAnalysis(Type,"V0M",Low,High,"V0CosPA");
  //makeCutAnalysis(Type,"V0M",Low,High,"DCANegToPV");
  /*makeCutAnalysis(Type,"V0M",Low,High,"DCAPosToPV");
  makeCutAnalysis(Type,"V0M",Low,High,"DCAV0Daughters");
  makeCutAnalysis(Type,"V0M",Low,High,"ProperLifetime");*/
  //makeCutAnalysis(Type,"V0M",Low,High,"NumberOfCrossedRows");
  //makeCutAnalysis(Type,"V0M",Low,High,"NumberOfCrossedRowsOverFindable"); 
  /*makeCutAnalysis(Type,"V0M",Low,High,"TPCPIDNSigmas");
  makeCutAnalysis(Type,"V0M",Low,High,"SigmaForSignalExtraction");*/

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
  Double_t truecut = 0.0;
  if(lWhichSystVar == "V0Radius") {
    truecut = 0.5;
    loosestcut = 0.300; 
    veryloosecut = 0.300;
    loosecut = 0.4;
    tightcut = 0.6;
    verytightcut = 0.7;
    minrange = 0.00001;
    maxrange = 1.00;
  }
  if(lWhichSystVar == "V0CosPA") {
    truecut = 0.995;
    loosestcut = 0.993; 
    veryloosecut = 0.993;
    loosecut = 0.994;
    tightcut = 0.996;
    verytightcut = 0.997;
    minrange = 0.99;
    maxrange = 1.0;
  }
  if(lWhichSystVar == "DCAPosToPV" || lWhichSystVar == "DCANegToPV") {
    truecut = 0.060;
    loosestcut = 0.050; 
    veryloosecut = 0.050;
    loosecut = 0.055;
    tightcut = 0.070;
    verytightcut = 0.080;
    minrange = 0.04;
    maxrange = 0.1;
  }
  if(lWhichSystVar == "DCAV0Daughters") {
    truecut = 1.;
    loosestcut = 1.5; 
    veryloosecut = 1.5;
    loosecut = 1.25;
    tightcut = 0.75;
    verytightcut = 0.50;
    minrange = 0.0;
    maxrange = 2.0;
  }
  if(lWhichSystVar == "ProperLifetime") {
    truecut = 30;
    loosestcut = 40.; 
    loosecut = 40.;
    tightcut = 20.;
    minrange = 0.;
    maxrange = 60.;
  }
   if(lWhichSystVar == "SigmaForSignalExtraction") {
    truecut = 6;
    loosestcut = 7.; 
    loosecut = 7.;
    tightcut = 5.;
    verytightcut = 4.;
    minrange = 1.;
    maxrange = 8.;
  }
   if(lWhichSystVar == "NumberOfCrossedRows") {
    truecut = 70;
    loosestcut = 70.; 
    tightcut = 75.;
    verytightcut = 80.;
    minrange = 60.;
    maxrange = 90.;
  }
   if(lWhichSystVar == "TPCPIDNSigmas") {
    truecut = 5;
    loosestcut = 7.0; 
    loosecut = 7.0;
    tightcut = 6.0;
    verytightcut = 4.0;
    minrange = 1.0;
    maxrange = 8.0;
  }
   if(lWhichSystVar == "NumberOfCrossedRowsOverFindable") {
    truecut = 0.8;
    loosestcut = 0.95; 
    verytightcut = 0.95;
    minrange = 0.75;
    maxrange = 1.0;
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
  lDataFilename[0] = lSystFile + lWhichSystVar + "-loosest.root"; //[0] is loosest cut
  for (int i = 1; i < nfiles; i++){ 
    lDataFilename[i] = lSystFile + lWhichSystVar + Form("-%i.root",i);
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
    lHistPtMC[i] = (TH1D *) InputFile[i]->Get("lInvMassMC/lInvMassMCRawData/fHistPtRawMC");
  }

  //Fill Sgn Loss Fraction
  Double_t binlimits[nfiles];
  binlimits[0] = minrange;
  for (int i = 1; i < nfiles; i++){
  	binlimits[i] = binlimits[i-1] + (maxrange-minrange)/(nfiles-1);
  }
  Long_t binnumb = sizeof(binlimits)/sizeof(Double_t) - 1;

  TH1D* SgnLossFraction = new TH1D("SgnLossFraction", Form(" ; %s; Fraction of signal loss",lWhichSystVar.Data()), binnumb, binlimits);
  TH1D* SgnLossFractionMC = new TH1D("SgnLossFractionMC", Form(" ; %s; Fraction of signal loss",lWhichSystVar.Data()), binnumb, binlimits);

  const int nbins = lHistPtRaw[0]->GetNbinsX();
  Double_t maxy = 0.0;
  Double_t errori = 0.0;
  Double_t error0 = 0.0;
  Double_t erroriMC = 0.0;
  Double_t error0MC = 0.0;

  for (int i = 1; i < binnumb; i++){
     	
    Double_t entry = 1. - TMath::Abs((double)(lHistPtRaw[i]->Integral(1,nbins))/(double)(lHistPtRaw[0]->Integral(1,nbins)));
    Double_t entryMC = 1. - TMath::Abs((double)(lHistPtMC[i]->Integral(1,nbins))/(double)(lHistPtMC[0]->Integral(1,nbins)));
    cout << "Real " << lHistPtRaw[i]->Integral(1,nbins) << "loose " << lHistPtRaw[0]->Integral(1,nbins) << endl;

    Double_t Inti = (double)lHistPtRaw[i]->IntegralAndError(1,nbins,errori);
    Double_t Int0 = (double)lHistPtRaw[0]->IntegralAndError(1,nbins,error0); 
    Double_t IntiMC = (double)lHistPtMC[i]->IntegralAndError(1,nbins,erroriMC);
    Double_t Int0MC = (double)lHistPtMC[0]->IntegralAndError(1,nbins,error0MC);    
    
    Double_t error = (double) (TMath::Sqrt(errori*errori/(Inti*Inti) + error0*error0/(Int0*Int0))*entry);
    Double_t errorMC = (double) (TMath::Sqrt(erroriMC*erroriMC/(IntiMC*IntiMC) + error0MC*error0MC/(Int0MC*Int0MC))*entryMC);

    maxy = max(entry,maxy);

    if (entry>0) SgnLossFraction->SetBinContent(i,entry);
    else SgnLossFraction->SetBinContent(i,0.000000001);
    if (entryMC>0) SgnLossFractionMC->SetBinContent(i,entryMC);
    else SgnLossFractionMC->SetBinContent(i,0.00000001);
   
    SgnLossFraction->SetBinError(i,error);
    SgnLossFractionMC->SetBinError(i,errorMC);
  }
  
  //Draw Options
  SgnLossFraction->SetMarkerStyle(20);
  SgnLossFraction->SetMarkerSize(0.8);
  SgnLossFraction->SetMarkerColor(kBlack);
  SgnLossFractionMC->SetMarkerStyle(20);
  SgnLossFractionMC->SetMarkerSize(0.8);
  SgnLossFractionMC->SetMarkerColor(kRed);
  SgnLossFractionMC->SetLineColor(kRed);

  TCanvas* c = new TCanvas;
  SgnLossFraction->SetStats(kFALSE);  
  SgnLossFraction->GetYaxis()->SetTitleSize(0.05);
  SgnLossFraction->GetYaxis()->SetTitleOffset(0.9);
  SgnLossFraction->GetXaxis()->SetTitleSize(0.05);
  SgnLossFraction->GetXaxis()->SetTitleOffset(0.8);
  SgnLossFraction->Draw();
  SgnLossFractionMC->Draw("SAME");
  //Cut Markers
  TLine* fvl = new TLine(veryloosecut,-0.001,veryloosecut,maxy);
  fvl->SetLineColor(kRed);
  fvl->SetLineWidth(2);
  fvl->SetLineStyle(7);
  TLine* fl = new TLine(loosecut,-0.001,loosecut,maxy);
  fl->SetLineColor(kBlue+1);
  fl->SetLineWidth(2);
  fl->SetLineStyle(7);
  TLine* ft = new TLine(tightcut,-0.001,tightcut,maxy);
  ft->SetLineColor(kGreen+1);
  ft->SetLineWidth(2);
  ft->SetLineStyle(7);
  TLine* fvt = new TLine(verytightcut,-0.001,verytightcut,maxy);
  fvt->SetLineColor(kAzure+8);
  fvt->SetLineWidth(2);
  fvt->SetLineStyle(7);
  TLine* ftrue = new TLine(truecut,-0.001,truecut,maxy);
  ftrue->SetLineColor(kBlack);
  ftrue->SetLineWidth(2);
  ftrue->SetLineStyle(7);

  if(veryloosecut!=0) fvl->Draw("SAME");
  if(loosecut!=0) fl->Draw("SAME");
  if(tightcut!=0) ft->Draw("SAME");
  if(verytightcut!=0) fvt->Draw("SAME");
  ftrue->Draw("SAME");
  c->SaveAs(Form("CS%s.png",lWhichSystVar.Data()));
  return; 
}
  