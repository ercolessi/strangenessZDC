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

void Test_Calibration()
{
  Float_t ZDCP1[5], ZDCP1Sum, ZDCN1[5], ZDCN1Sum, ZDCP2[5], ZDCP2Sum, ZDCN2[5], ZDCN2Sum, EZDC, LogE;
  Double_t cP1[5];
  Double_t cP2[5];
  Double_t cN1[5];
  Double_t cN2[5];
  const Int_t bin = 1000;
  const Int_t range = 4000;
      
  TFile* Read = new TFile ("LeadingMerged12b.root");
  TTree * T = (TTree *)Read->Get("tEvents");
  
  T->SetBranchAddress("adcZDCP1",ZDCP1);
  T->SetBranchAddress("adcZDCN1",ZDCN1);
  T->SetBranchAddress("adcZDCP2",ZDCP2);
  T->SetBranchAddress("adcZDCN2",ZDCN2);

  TH1D * hEZDC = new TH1D("hEZDC","EZDC", bin ,-1,range);
  TH1D * hEZDClog = new TH1D("hEZDClog","EZDC", bin, 0.,4.);
  TH1D * hEZDC8 = new TH1D("hEZDC8","EZDC", bin ,-1,range);
  TH1D * hEZDC6 = new TH1D("hEZDC6","EZDC", bin ,-1, range);
  TH1D * hEZDC4 = new TH1D("hEZDC4","EZDC", bin ,-1,range);
  TH1D * hEZDC2 = new TH1D("hEZDC2","EZDC", bin ,-1,range);
  TH1D * hEZDC1 = new TH1D("hEZDC1","EZDC", bin ,-1,range);
  TH1D * hCentr = new TH1D("hCentr","EZDC", bin ,-1,range);
  TH1D * hEZDC8log = new TH1D("hEZDC8log","EZDC", bin ,0.,4.);
  TH1D * hEZDC6log = new TH1D("hEZDC6log","EZDC", bin ,0.,4.);
  TH1D * hEZDC4log = new TH1D("hEZDC4log","EZDC", bin ,0.,4.);
  TH1D * hEZDC2log = new TH1D("hEZDC2log","EZDC", bin ,0.,4.);
  TH1D * hEZDC1log = new TH1D("hEZDC1log","EZDC", bin ,0.,4.);
  TH1D * hCentrlog = new TH1D("hCentrlog","EZDC", bin ,-1,range);
  TH1D * hEZDC_perc = new TH1D("hEZDC_perc","EZDC_perc", 100 ,-1.,101.);
  TH1D * hEZDC_perclog = new TH1D("hEZDC_perclog","EZDC_perc", 100 ,-1.,101.);
  
  
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
      
      LogE = TMath::Log10(TMath::Abs(EZDC)+1);
      hEZDClog->Fill(LogE);
    }

  TFile* Write = new TFile("Centrality12b.root", "recreate");

  hEZDC->Write();
  hEZDClog->Write();
  
  TH1D * hCumu = new TH1D(*hEZDC);
  hCumu->SetName("Cumulative distribution");
  hCumu->Reset();  
  TH1D * hCumulog = new TH1D(*hEZDClog);
  hCumulog->SetName("Cumulative distribution Elog");
  hCumulog->Reset();
  
  Double_t IntTOT = hEZDC->Integral(1,bin);// == total number of events
  Double_t IntTOTlog = hEZDClog->Integral(1,bin); 
  Double_t val = 0;
  Double_t vallog = 0;
  
  for (Int_t k=1; k<bin; k++) //loop through the bins
    {
      val+= (hEZDC->GetBinContent(k))/IntTOT;
      vallog+= (hEZDClog->GetBinContent(k))/IntTOTlog;
      hCumu->SetBinContent(k,val);
      hCumulog->SetBinContent(k,vallog);   
    }
  
  hCumu->Write();
  hCumulog->Write();

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
      LogE = TMath::Log10(TMath::Abs(EZDC)+1);
    
      Double_t E = hCumu->Interpolate(EZDC);
      hEZDC_perc->Fill(E*100);
      if(E>0.9) hEZDC8->Fill(EZDC);
      if(E>0.7 && E<0.8) hEZDC6->Fill(EZDC);
      if(E>0.5 && E<0.6) hEZDC4->Fill(EZDC);
      if(E>0.3 && E<0.4) hEZDC2->Fill(EZDC);
      if(E>0.1 && E<0.2) hEZDC1->Fill(EZDC);

      Double_t Elog = hCumulog->Interpolate(LogE);
      hEZDC_perclog->Fill(Elog*100);
      if(Elog>0.9) hEZDC8log->Fill(LogE);
      if(Elog>0.7 && Elog<0.8) hEZDC6log->Fill(LogE);
      if(Elog>0.5 && Elog<0.6) hEZDC4log->Fill(LogE);
      if(Elog>0.3 && Elog<0.4) hEZDC2log->Fill(LogE);
      if(Elog>0.1 && Elog<0.2) hEZDC1log->Fill(LogE);
    }

 
 hEZDC8->SetFillColor(16);
 hEZDC8->Write();
 hEZDC6->SetFillColor(16);
 hEZDC6->Write();
 hEZDC4->SetFillColor(16);
 hEZDC4->Write();
 hEZDC2->SetFillColor(16);
 hEZDC2->Write();
 hEZDC1->SetFillColor(16);
 hEZDC1->Write();
 hEZDC_perc->Write();
 hEZDC8log->SetFillColor(16);
 hEZDC8log->Write();
 hEZDC6log->SetFillColor(16);
 hEZDC6log->Write();
 hEZDC4log->SetFillColor(16);
 hEZDC4log->Write();
 hEZDC2log->SetFillColor(16);
 hEZDC2log->Write();
 hEZDC1log->SetFillColor(16);
 hEZDC1log->Write();
 hEZDC_perclog->Write();
  // Write->Close();

}
