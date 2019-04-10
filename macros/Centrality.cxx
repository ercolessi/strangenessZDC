#include <TFile.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TClonesArray.h>
#include <TLatex.h>

// STL includes
#include <iostream>
#include <fstream>
using std::cout;
using std::endl;

void BinLogX(TH1 *h);

void Centrality()
{
  
  Float_t ZDCP1[5], ZDCP1Sum, ZDCN1[5], ZDCN1Sum, ZDCP2[5], ZDCP2Sum, ZDCN2[5], ZDCN2Sum, EZDC;
  Double_t cP1[5];
  Double_t cP2[5];
  Double_t cN1[5];
  Double_t cN2[5];
  const Int_t bin = 100000;
  
  TFile * Calib = new TFile ("Calib12b.root");
 
  TH1D * hcP1_1 = (TH1D *)Calib->Get("hcP1_1");
  TH1D * hcP1_2 = (TH1D *)Calib->Get("hcP1_2");
  TH1D * hcP1_3 = (TH1D *)Calib->Get("hcP1_3");
  TH1D * hcP1_4 = (TH1D *)Calib->Get("hcP1_4");
  TH1D * hcP2_1 = (TH1D *)Calib->Get("hcP2_1");
  TH1D * hcP2_2 = (TH1D *)Calib->Get("hcP2_2");
  TH1D * hcP2_3 = (TH1D *)Calib->Get("hcP2_3");
  TH1D * hcP2_4 = (TH1D *)Calib->Get("hcP2_4");
  TH1D * hcN1_1 = (TH1D *)Calib->Get("hcN1_1");
  TH1D * hcN1_2 = (TH1D *)Calib->Get("hcN1_2");
  TH1D * hcN1_3 = (TH1D *)Calib->Get("hcN1_3");
  TH1D * hcN1_4 = (TH1D *)Calib->Get("hcN1_4");
  TH1D * hcN2_1 = (TH1D *)Calib->Get("hcN2_1");
  TH1D * hcN2_2 = (TH1D *)Calib->Get("hcN2_2");
  TH1D * hcN2_3 = (TH1D *)Calib->Get("hcN2_3");
  TH1D * hcN2_4 = (TH1D *)Calib->Get("hcN2_4");

  cP1[1] = hcP1_1->GetMean();
  cP1[2] = hcP1_2->GetMean();
  cP1[3] = hcP1_3->GetMean();
  cP1[4] = hcP1_4->GetMean();
  cP2[1] = hcP2_1->GetMean();
  cP2[2] = hcP2_2->GetMean();
  cP2[3] = hcP2_3->GetMean();
  cP2[4] = hcP2_4->GetMean();
  cN1[1] = hcN1_1->GetMean();
  cN1[2] = hcN1_2->GetMean();
  cN1[3] = hcN1_3->GetMean();
  cN1[4] = hcN1_4->GetMean();
  cN2[1] = hcN2_1->GetMean();
  cN2[2] = hcN2_2->GetMean();
  cN2[3] = hcN2_3->GetMean();
  cN2[4] = hcN2_4->GetMean();

  Calib->Close();//-----------------------------------------------
  
  TFile* Read = new TFile ("LeadingMerged12b.root");
  TTree * T = (TTree *)Read->Get("tEvents");
  
  T->SetBranchAddress("adcZDCP1",ZDCP1);
  T->SetBranchAddress("adcZDCN1",ZDCN1);
  T->SetBranchAddress("adcZDCP2",ZDCP2);
  T->SetBranchAddress("adcZDCN2",ZDCN2);

  TH1D * hEZDC = new TH1D("hEZDC","EZDC", bin ,-1,4000);
  //BinLogX(hEZDC);
   
  cout << "Looping on Tree entries: " <<T->GetEntries() << " entries..."<< endl;

  for(Int_t i=0; i<T->GetEntries();i++)//loop over Tree entries
    {
      T->GetEvent(i);//get Tree entries

      if (i%1000000 == 0)
	cout<<"Analizing event = "<<i+1<<" / "<< T->GetEntries() << endl;
      
      for (Int_t l=1; l<5; l++)
	{
	  ZDCP1[l]/=cP1[l];
	  ZDCP2[l]/=cP2[l];
	  ZDCN1[l]/=cN1[l];
	  ZDCN2[l]/=cN2[l];
	}
      
      //Sums
      ZDCP1Sum=ZDCP1[1]+ZDCP1[2]+ZDCP1[3]+ZDCP1[4];
      ZDCP2Sum=ZDCP2[1]+ZDCP2[2]+ZDCP2[3]+ZDCP2[4];
      ZDCN1Sum=ZDCN1[1]+ZDCN1[2]+ZDCN1[3]+ZDCN1[4];
      ZDCN2Sum=ZDCN2[1]+ZDCN2[2]+ZDCN2[3]+ZDCN2[4];

      EZDC = ZDCP1[0]+ZDCP2[0]+ZDCN1[0]+ZDCN2[0];
      hEZDC->Fill(EZDC);
    }

  TFile* Write = new TFile ("Centrality12b.root", "recreate");
  
  hEZDC->Write();

  Double_t IntTOT = hEZDC->Integral(0,bin); // == total number of events
  const Int_t nclass = 10; // number of percentile classes to divide in 
  Double_t  E, x[nclass], kbin[nclass], y[nclass], Int[bin], ratio[bin],D[bin]; //some variables
  
  for (Int_t k=1; k<bin; k++) //loop through the bins
    {
      E = hEZDC->GetBinCenter(k);
      Int[k]= hEZDC->Integral(0,k); //cumulative sum of bins untill bin k
      ratio[k] = Int[k]/IntTOT; //fraction of the events in the rangee [0,k]
      D[k] = hEZDC->GetBinContent(k);//value of histo in bin k
      Double_t perc;//=n+1/nclass    
      for (Int_t n = 0; n<nclass; n++)
	{
	  perc = ((Double_t)(n+1))/((Double_t)nclass);
	  if(ratio[k] <= (perc + (0.01/nclass)) && ratio[k] >= (perc - (0.01/nclass)))    {
	    x[n]=hEZDC->GetBinCenter(k);
	    kbin[n]=k;
	    y[n]=(ratio[k-1]+D[k]/IntTOT)*100;
	  }	    
	}
      
      x[nclass-1] = hEZDC->GetBinCenter(bin);     
    }
  

  TGraph * Centr = new TGraph (nclass,x,y);
  // BinLogX(Centr);
  Centr->SetMarkerStyle(21);
  Centr->SetMarkerColor(2);
  
  Centr->Write();

  
  for (Int_t l =0; l<nclass;l++)
    {
      cout << "y = " << y[l] << endl;
      cout << "x = " << x[l] << endl;
      cout << "kbin = " << kbin[l] << endl;
    }
  
  Write->Close();
  
}

void BinLogX(TH1 *h) {

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
