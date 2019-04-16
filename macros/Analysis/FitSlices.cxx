{
  Double_t sigmaP1[100];//==bin
  Double_t sigmaP2[100];
  Double_t sigmaN1[100];
  Double_t sigmaN2[100];
  Int_t y=1;
  Int_t bin = 100;
  Float_t ZDCP1[5], ZDCP1Sum, ZDCN1[5], ZDCN1Sum, ZDCP2[5], ZDCP2Sum, ZDCN2[5], ZDCN2Sum,dP1,meanP1, dN1, meanN1,dP2,meanP2, dN2, meanN2;
  Double_t cP1[5];
  Double_t cP2[5];
  Double_t cN1[5];
  Double_t cN2[5];
  
  /* TFile * File = new TFile ("LHC12b.root");
  
  TH2D * resP1 = (TH2D *)File->Get("resP1");
  TH2D * resN1 = (TH2D *)File->Get("resN1");
  TH2D * resP2 = (TH2D *)File->Get("resP2");
  TH2D * resN2 = (TH2D *)File->Get("resN2");
  TH1D * resP1_1 = (TH1D *)File->Get("resP1_1");
  TH1D * resN1_1 = (TH1D *)File->Get("resN1_1");
  TH1D * resP2_1 = (TH1D *)File->Get("resP2_1");
  TH1D * resN2_1 = (TH1D *)File->Get("resN2_1");
  TH1D * resP1_2 = (TH1D *)File->Get("resP1_2");
  TH1D * resN1_2 = (TH1D *)File->Get("resN1_2");
  TH1D * resP2_2 = (TH1D *)File->Get("resP2_2");
  TH1D * resN2_2 = (TH1D *)File->Get("resN2_2");

  for (Int_t i =0; i<100; i++)//loop sui bin x
    {
      sigmaP1[i] = resP1_2->GetBinContent(i);//prendo la sigma dei dati y su quel bin i di x
    }
  
  File->Close();//------------------------------------------
  */
   
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

  TH2D * resP1= new TH2D("resP1","Resolution P1;sum-common;(sum+common)/2;",bin,-1,4000,bin,-2000,2000);
  TH2D * resN1= new TH2D("resN1","Resolution N1;sum-common;sum+comm /2;",bin,-1,4000,bin,-2000,2000);
  TH2D * resP2= new TH2D("resP2","Resolution P2;sum-common;(sum+common)/2;",bin,-1,4000,bin,-2000,2000);
  TH2D * resN2= new TH2D("resN2","Resolution N2;sum-common;sum+comm /2 ;",bin,-1,4000,bin,-2000,2000);
  TH2D * hresP1= new TH2D("hresP1","Resolution P1;sum-common;(sum+common)/2;",bin,-1,4000,bin,-200,200);
  TH2D * hresN1= new TH2D("hresN1","Resolution N1;sum-common;sum+comm /2;",bin,-1,4000,bin,-200,200);
  TH2D * hresP2= new TH2D("hresP2","Resolution P2;sum-common;(sum+common)/2;",bin,-1,4000,bin,-200,200);
  TH2D * hresN2= new TH2D("hresN2","Resolution N2;sum-common;sum+comm /2 ;",bin,-1,4000,bin,-200,200);
  
  cout << "Looping on Tree entries: " <<T->GetEntries() << " entries..."<< endl;

  for(Int_t p=0; p<T->GetEntries();p++)//loop over Tree entries
    {
      T->GetEvent(p);//get Tree entries
      
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
      
      dP1=ZDCP1Sum-ZDCP1[0];
      meanP1=(ZDCP1Sum+ZDCP1[0])/2;
      dN1=ZDCN1Sum-ZDCN1[0];
      meanN1=(ZDCN1Sum+ZDCN1[0])/2;
      dP2=ZDCP2Sum-ZDCP2[0];
      meanP2=(ZDCP2Sum+ZDCP2[0])/2;
      dN2=ZDCN2Sum-ZDCN2[0];
      meanN2=(ZDCN2Sum+ZDCN2[0])/2;
      
      resP1->Fill(meanP1,dP1);
      resN1->Fill(meanN1,dN1);
      resP2->Fill(meanP2,dP2);
      resN2->Fill(meanN2,dN2);
      
    } // Tree loop

  TF1* gaus = new TF1 ("gaus", "gaus", -100,100);

  TFile* Write = new TFile ("Resolution12b.root", "recreate");
  resP1->Write();
  resP1->FitSlicesY(gaus);
  TH1D * resP1_1 = (TH1D*)gDirectory->Get("resP1_1");
  TH1D * resP1_2 = (TH1D*)gDirectory->Get("resP1_2");
  resP2->Write();
  resP2->FitSlicesY(gaus);
  TH1D * resP2_1 = (TH1D*)gDirectory->Get("resP2_1");
  TH1D * resP2_2 = (TH1D*)gDirectory->Get("resP2_2");
  resN1->Write();
  resN1->FitSlicesY(gaus);
  TH1D * resN1_1 = (TH1D*)gDirectory->Get("resN1_1");
  TH1D * resN1_2 = (TH1D*)gDirectory->Get("resN1_2");
  resN2->Write();
  resN2->FitSlicesY(gaus);
  TH1D * resN2_1 = (TH1D*)gDirectory->Get("resN2_1");
  TH1D * resN2_2 = (TH1D*)gDirectory->Get("resN2_2"); 

  for (Int_t i =0; i<100; i++)//loop sui bin x
    {
      sigmaP1[i] = resP1_2->GetBinContent(i);//prendo la sigma dei dati y su quel bin i di x
      sigmaP2[i] = resP2_2->GetBinContent(i);
      sigmaN1[i] = resN1_2->GetBinContent(i);
      sigmaN2[i] = resN2_2->GetBinContent(i);
    }

  for(Int_t p=0; p<T->GetEntries();p++)//loop over Tree entries
    {
      T->GetEvent(p);//get Tree entries

      if (p%1000000 == 0)
	cout<<"Analizing event = "<<p+1<<" / "<< T->GetEntries() << endl;
      
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
      
      dP1=ZDCP1Sum-ZDCP1[0];
      meanP1=(ZDCP1Sum+ZDCP1[0])/2;
      dN1=ZDCN1Sum-ZDCN1[0];
      meanN1=(ZDCN1Sum+ZDCN1[0])/2;
      dP2=ZDCP2Sum-ZDCP2[0];
      meanP2=(ZDCP2Sum+ZDCP2[0])/2;
      dN2=ZDCN2Sum-ZDCN2[0];
      meanN2=(ZDCN2Sum+ZDCN2[0])/2;
   
      for (Int_t k=-1; k<3999; k+=40)
	{
	  if (meanP1<(k+40) && meanP1>=(k)) {hresP1->Fill(meanP1,dP1/sigmaP1[y]);};
	  if (meanN1<(k+40) && meanN1>=(k)) {hresN1->Fill(meanN1,dN1/sigmaN1[y]);};
	  if (meanP2<(k+40) && meanP2>=(k)) {hresP2->Fill(meanP2,dP2/sigmaP2[y]);};
	  if (meanN2<(k+40) && meanN2>=(k)) {hresN2->Fill(meanN2,dN2/sigmaN2[y]);};
	  y=+1;
	}
    } // Tree loop
  
  hresP1->Write();
  hresP1->FitSlicesY(gaus);
  TH1D * hresP1_1 = (TH1D*)gDirectory->Get("hresP1_1");
  hresP1_1->Write();
  TH1D * hresP1_2 = (TH1D*)gDirectory->Get("hresP1_2"); 
  hresP1_2->Scale(3);
  hresP1_2->SetLineColor(5);
  hresP1_2->Write();
  hresP2->Write();
  hresP2->FitSlicesY(gaus);
  TH1D * hresP2_1 = (TH1D*)gDirectory->Get("hresP2_1");
  hresP2_1->Write();
  TH1D * hresP2_2 = (TH1D*)gDirectory->Get("hresP2_2"); 
  hresP2_2->Scale(3);
  hresP2_2->Write();
  hresN1->Write();
  hresN1->FitSlicesY(gaus);
  TH1D * hresN1_1 = (TH1D*)gDirectory->Get("hresN1_1");
  hresN1_1->Write();
  TH1D * hresN1_2 = (TH1D*)gDirectory->Get("hresN1_2"); 
  hresN1_2->Scale(3);
  hresN1_2->Write();
  hresN2->Write();
  hresN2->FitSlicesY(gaus);
  TH1D * hresN2_1 = (TH1D*)gDirectory->Get("hresN2_1");
  hresN2_1->Write();
  TH1D * hresN2_2 = (TH1D*)gDirectory->Get("hresN2_2"); 
  hresN2_2->Scale(3);
  hresN2_2->Write();

  Write->Close();
}
