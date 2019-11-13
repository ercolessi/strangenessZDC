
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

void ShowMinv() {

  const int nbins = 23;

  TFile *f = new TFile("Omega/PtSpectraOmega.root");
  TH2D* hMinvPt = (TH2D*) f->Get("MinvPt4");
 
  //Fit function for invariant mass
   TF1 *ffit = new TF1("ffit","gaus(0)/sqrt(2*TMath::Pi())/[2]+pol2(3)",1.62,1.74);
  ffit->SetParameter(0,0.0002);
  ffit->SetParameter(1,1.67251);
  ffit->SetParLimits(1,1.672,1.673);
  ffit->SetParameter(2,2.34086e-03);
  ffit->SetParLimits(2,2.4E-3,2.3E-3);
  ffit->SetNpx(1000);
  
  TH1D* hProjMinv = hMinvPt->ProjectionY("hProjMinv",5,5);
  cout<< "PTcenter = " << hMinvPt->GetXaxis()->GetBinCenter(5)<< endl;
  
  //Fit
  hProjMinv->Fit(ffit,"","",1.62,1.74);
  ffit->SetParameter(3,0);
  ffit->SetParameter(4,0);
  ffit->SetParameter(5,0);
}


