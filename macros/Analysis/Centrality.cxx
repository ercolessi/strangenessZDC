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

void Centrality()
{
  
  Float_t ZDCP1[5], ZDCP1Sum, ZDCN1[5], ZDCN1Sum, ZDCP2[5], ZDCP2Sum, ZDCN2[5], ZDCN2Sum, EZDC;
  Double_t cP1[5];
  Double_t cP2[5];
  Double_t cN1[5];
  Double_t cN2[5];
  const Int_t bin = 1000;
  const Int_t range = 4000;
  
  TFile * Calib = new TFile ("Calib18i.root");
 
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
  
  TFile* Read = new TFile ("18ipass1.root");
  TTree * T = (TTree *)Read->Get("tEvents");
  
  T->SetBranchAddress("adcZDCP1",ZDCP1);
  T->SetBranchAddress("adcZDCN1",ZDCN1);
  T->SetBranchAddress("adcZDCP2",ZDCP2);
  T->SetBranchAddress("adcZDCN2",ZDCN2);

  TH1D * hEZDC = new TH1D("hEZDC","EZDC", bin ,-1,range);
  TH1D * hEZDC8 = new TH1D("hEZDC8","EZDC", bin ,-1,range);
  TH1D * hEZDC6 = new TH1D("hEZDC6","EZDC", bin ,-1, range);
  TH1D * hEZDC4 = new TH1D("hEZDC4","EZDC", bin ,-1,range);
  TH1D * hEZDC2 = new TH1D("hEZDC2","EZDC", bin ,-1,range);
  TH1D * hCentr = new TH1D("hCentr","EZDC", bin ,-1,range);
  TH1D * hEZDC_perc = new TH1D("hEZDC_perc","EZDC_perc", bin ,-1.,100.);
   
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

  TFile* Write = new TFile ("Centrality18i.root", "recreate");
  
  hEZDC->Write();

  Double_t IntTOT = hEZDC->Integral(0,bin); // == total number of events
  const Int_t nclass = 10; // number of percentile classes to divide in 
  Double_t  E, x[nclass], kbin[nclass], y[nclass], Int[bin], ratio[bin],D[bin]; //some variables
  
  for (Int_t k=1; k<bin; k++) //loop through the bins
    {
      E = hEZDC->GetBinCenter(k);
      Int[k]= hEZDC->Integral(0,k);//cumulative sum of bins untill bin k
      ratio[k] = Int[k]/IntTOT; //fraction of the events in the rangee [0,k]
      D[k] = hEZDC->GetBinContent(k);//value of histo in bin k
      Double_t perc;//=n+1/nclass
      
      for (Int_t n = 0; n<nclass; n++)
	{
	  perc = ((Double_t)(n+1))/((Double_t)nclass);
	  Int_t const slices = 100;
	  Double_t IntRatio;
	  for (Int_t l=0; l<slices; l++)
	    {
	      IntRatio = (ratio[k-1]+((D[k]*range*l)/(IntTOT*bin*slices)));
	      if ((IntRatio <= (perc + 5*(0.001/nclass))) && (IntRatio >= (perc - 5*(0.001/nclass)))) {
		x[n]=hEZDC->GetBinCenter(k);
		kbin[n]=k;
		y[n]=IntRatio;
		break;
	      }
	      
	      if (IntRatio >=perc && IntRatio < ( ((Double_t)(n+2))/((Double_t)nclass))) hCentr->Fill(x[n]);
	    }
	}
      x[nclass-1] = hEZDC->GetBinCenter(bin);     
    }
 
  for (Int_t j =1; j<nclass;j++)
    {
      cout << "y = " << y[j] << endl;
      cout << "x = " << x[j] << endl;
      Int_t r = 0;
      for(r=0; r<(j*10); r++) hCentr->Fill(x[j-1]);
    }

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
      if(EZDC>x[8]) hEZDC8->Fill(EZDC);
      if(EZDC>x[6] && EZDC<x[7]) hEZDC6->Fill(EZDC);
      if(EZDC>x[4] && EZDC<x[5]) hEZDC4->Fill(EZDC);
      if(EZDC>x[2] && EZDC<x[3]) hEZDC2->Fill(EZDC);

      for(Int_t k = 0; k<9; k++)
	{
	  Double_t p = hCentr->Interpolate(EZDC);
	  hEZDC_perc->Fill(p);
	}
    }

    hEZDC8->SetFillColor(16);
    hEZDC8->Write();
    hEZDC6->SetFillColor(16);
    hEZDC6->Write();
    hEZDC4->SetFillColor(16);
    hEZDC4->Write();
    hEZDC2->SetFillColor(16);
    hEZDC2->Write();
    hCentr->Write();
    hEZDC_perc->Write();
    
  Write->Close();
















   /*  if(ratio[k] <= (perc + (0.1/nclass)) && ratio[k] >= (perc - (0.1/nclass)))    {
	    cout<< hEZDC->Interpolate(perc) << endl;
	    z[n]=hEZDC->GetBinCenter(k);
	    zkbin[n]=k;
	    zy[n]=(ratio[k-1]+D[k]/IntTOT)*100;
	     }
	  */
}

