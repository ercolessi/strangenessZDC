#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

void ZDCCentrality()
{  
  //Output file name
   const char *outputname = "EZDCprova.root";
  
  //Declare some variables
  Float_t ZDCP1[5], ZDCP1Sum, ZDCN1[5], ZDCN1Sum, ZDCP2[5], ZDCP2Sum, ZDCN2[5], ZDCN2Sum, EZDC, LogE;
  const Int_t bin = 1000;
  const Double_t range = 4000;
  const Double_t Logrange = 4.; 
  
  //Open Leading File and Set Variables 
  TFile* Read = new TFile ("LeadingMerged12f.root");
  TTree * T = (TTree *)Read->Get("tEvents");
  
  T->SetBranchAddress("adcZDCP1",ZDCP1);
  T->SetBranchAddress("adcZDCN1",ZDCN1);
  T->SetBranchAddress("adcZDCP2",ZDCP2);
  T->SetBranchAddress("adcZDCN2",ZDCN2);
  
  TH1F * hESumZDC = new TH1F("hESumZDC","hSumEZDC", bin ,0.,Logrange);
  
  cout << "Looping on Tree entries: " <<T->GetEntries() << " entries..."<< endl;
  for(Int_t i=0; i<T->GetEntries();i++)//loop over Tree entries
    {
      T->GetEvent(i);//get Tree entries

      if (i%1000000 == 0)
	cout<<"Analizing event = "<<i+1<<" / "<< T->GetEntries() << endl;
      
      EZDC = ZDCP1[0]+ZDCP2[0]+ZDCN1[0]+ZDCN2[0];    
      LogE = TMath::Log10(TMath::Abs(EZDC)+1);
      hESumZDC->Fill(LogE);        
    }
  
  TFile* Write = new TFile (outputname, "recreate");
  hESumZDC->Write();
  
  //Cumulative distribution gives centrality
  TH1F * hCumulative = new TH1F(*hESumZDC);
  hCumulative->SetName("hCumulative");
  hCumulative->Reset();  
  
  Double_t IntTOT = hESumZDC->Integral(1,bin);// == total number of events
  Double_t val = 0;
  
  for (Int_t k=1; k<bin; k++) //loop through the bins
    {
      val+= (hESumZDC->GetBinContent(k))/IntTOT;
      hCumulative->SetBinContent(k,val);
    }
  
  hCumulative->Write(); 
  Write->Close();
}

