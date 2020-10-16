#include <iostream>
#include <fstream>
//#define max(a,b) (a>b ? a : b)

double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr );

void makeCutAnalysis(TString lV0Type = "Lambda", TString lWhichEstimator = "V0M", Double_t lMultBoundLo = 0.0, Double_t lMultBoundHi = 100.0, TString lWhichSystVar = "V0Radius");

//Main Function
void PlotCutAnalysisV0s(TString Type = "Lambda", Double_t Low = 0.0, Double_t High = 100.0){
  
  makeCutAnalysis(Type,"V0M",Low,High,"V0Radius");

  return
}

//---------------------------------------------------------------
void makeSystPlotsV0s(TString lV0Type = "Lambda",
		      TString lWhichEstimator = "V0M",
		      Double_t lMultBoundLo = 0.0,
              Double_t lMultBoundHi = 100.0,
		      TString lWhichSystVar = "V0Radius"){

  Double_t loosestCut = 0.0;
  if (lWhichSystVar == "V0Radius") loosestCut = 0.300;
  
  cout<<"  "<<endl;
  cout<<"-------------------------------------------------------------"<<endl;
  cout<<"        Macro to make cut study of systematics for V0s       "<<endl;
  cout<<"-------------------------------------------------------------"<<endl;
  cout<<endl;

  //Set data files
  TString lSystFile = "Results-Systematics-V0M_000_100_ZDC_000_100_UseMCRatioFD"; 
  if(lV0Type == "Lambda"    ) lSystFile.Append("-Lambda-");
  if(lV0Type == "AntiLambda") lSystFile.Append("-AntiLambda-");

  const int nfiles = 20;
  TString lDataFilename[nfiles];
  for (int i = 0; i < nfiles; i++){ //[0] is loosest cut
    lDataFilename[i] = lSystFile + lWhichSystVar + Form("-%i.root",i);
  }
  
  TFile* InputFile[nfiles]; 
  for (int i = 0; i < nfiles; i++){
    InputFile[i] = new TFile(lDataFilename[i].Data());
  }
   
  //Get raw histograms from files
  TH1D* lHistPtRaw[nfiles];
  for (int i = 0; i < nfiles; i++){
    lHistPtRaw[i] = (TH1D *) InputFile[i]->Get("lInvMassReal/lInvMassRealRawData/fHistPtRaw");
  }


  //Fill Sgn Loss Fraction
  Double_t binlimits[20];
  for (int i = 0; i < 20; i++){
  	binlimits[i] = (loosestCut + 0.2*i);
  }
  //Long_t binnumb = sizeof(binlimits)/sizeof(Double_t) - 1;


  TH1D* SgnLossFraction = new TH1D("SgnLossFraction", Form(" ; Fraction of signal loss; %s",lWhichSystVar), 20, binlimits);

  const int nbins = lHistPtRaw[0]->GetNbinsX();
  for (int i = 1; i <= 20; i++){
  	Double_t entry = 1 - lHistPtRaw[i]->Integral(1,nbins)/lHistPtRaw[i]->Integral(1,nbins);
  	SgnLossFraction->SetBinContent(i,entry);
  }

  new TCanvas;
  SgnLossFraction->Draw();



}
  