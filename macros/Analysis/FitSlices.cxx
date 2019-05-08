{
  
  //Definition of variables
  const Int_t bin = 1000;
  Double_t sigmaP1[bin];
  Double_t sigmaP2[bin];
  Double_t sigmaN1[bin];
  Double_t sigmaN2[bin];
  Float_t ZDCP1[5], ZDCP1Sum, ZDCN1[5], ZDCN1Sum, ZDCP2[5], ZDCP2Sum, ZDCN2[5],
    ZDCN2Sum,dP1,meanP1, dN1, meanN1,dP2,meanP2, dN2, meanN2,EZDC, LogE;
  Float_t V0M, V0A, V0C, CL0, CL1, SPDClusters, SPDTracklets, ZNA, ZNC, ZNApp,
    ZNCpp, RefMult08_p, RefMult05_p, V0MEstimator_p, V0CEstimator_p, V0AEstimator_p,
    RefMult08_abs, RefMult05_abs, V0MEstimator_abs, V0CEstimator_abs, V0AEstimator_abs;
  Double_t cP1[5];
  Double_t cP2[5];
  Double_t cN1[5];
  Double_t cN2[5];

  //Get calibration coefficients for ZDC channels 
  TFile * Calib = new TFile ("Calib12f.root");
 
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

  Calib->Close();//-------------------------------------

  //Open the Tree and get the information
  TFile* Read = new TFile ("LeadingMerged12f.root");
  TTree * T = (TTree *)Read->Get("tEvents");
  
  T->SetBranchAddress("adcZDCP1",ZDCP1);
  T->SetBranchAddress("adcZDCN1",ZDCN1);
  T->SetBranchAddress("adcZDCP2",ZDCP2);
  T->SetBranchAddress("adcZDCN2",ZDCN2);
  T->SetBranchAddress("V0M",&V0M);
  T->SetBranchAddress("V0A",&V0A);
  T->SetBranchAddress("V0C",&V0C);
  T->SetBranchAddress("CL0",&CL0);
  T->SetBranchAddress("CL1",&CL1);
  T->SetBranchAddress("SPDClusters",&SPDClusters);
  T->SetBranchAddress("SPDTracklets",&SPDTracklets);
  T->SetBranchAddress("ZNA",&ZNA);
  T->SetBranchAddress("ZNC",&ZNC);
  T->SetBranchAddress("ZNApp",&ZNApp);
  T->SetBranchAddress("ZNCpp",&ZNCpp);
  T->SetBranchAddress("RefMult08_p",&RefMult08_p);
  T->SetBranchAddress("RefMult05_p",&RefMult05_p);
  T->SetBranchAddress("V0MEstimator_p",&V0MEstimator_p);
  T->SetBranchAddress("V0AEstimator_p",&V0AEstimator_p);
  T->SetBranchAddress("V0CEstimator_p",&V0CEstimator_p);
  T->SetBranchAddress("RefMult08_abs",&RefMult08_abs);
  T->SetBranchAddress("RefMult05_abs",&RefMult05_abs);
  T->SetBranchAddress("V0MEstimator_abs",&V0MEstimator_abs);
  T->SetBranchAddress("V0AEstimator_abs",&V0AEstimator_abs);
  T->SetBranchAddress("V0CEstimator_abs",&V0CEstimator_abs);

  //Histos
  //ZDCs histos on C hemisphere
  TH2D *hNCV0A = new TH2D("hNCV0A",";ZDCN-C[0];V0A;",100,-1,4000,100,-1,700);
  TH2D *hNCV0C = new TH2D("hNCV0C",";ZDCN-C[0];V0C;",100,-1,4000,100,-1,700);
  TH2D *hNCV0M = new TH2D("hNCV0M",";ZDCN-C[0];V0M;",100,-1,4000,100,-1,700);
  TH2D *hNCeta8 = new TH2D("hNCeta8",";ZDCN-C[0];RefMult08;",100,-1,4000,100,-1,700);
  TH2D *hPCV0A = new TH2D("hPCV0A",";ZDCP-C[0];V0A;",100,-1,4000,100,-1,700);
  TH2D *hPCV0C = new TH2D("hPCV0C",";ZDCP-C[0];V0C;",100,-1,4000,100,-1,700);
  TH2D *hPCV0M = new TH2D("hPCV0M",";ZDCP-C[0];V0M;",100,-1,4000,100,-1,700);
  TH2D *hPCeta8 = new TH2D("hPCeta8",";ZDCP-C[0];RefMult08;",100,-1,4000,100,-1,700);
  TH2D *totCV0M = new TH2D("totCV0M",";ZDCP-C[0]+ZDCN-C[0];V0M;",100,-1,4000,100,-1,700);
  TH2D *totCeta8 = new TH2D("totCeta8",";ZDCP-C[0]+ZDCN-C[0];RefMult08;",100,-1,4000,100,-1,700);
  TH2D *totCV0A = new TH2D("totCV0A",";ZDCP-C[0]+ZDCN-C[0];V0A;",100,-1,4000,100,-1,700);
  TH2D *totCV0C = new TH2D("totCV0C",";ZDCP-C[0]+ZDCN-C[0];V0C;",100,-1,4000,100,-1,700);
  //A hemisphere
  TH2D *hNAV0A = new TH2D("hNAV0A",";ZDCN-C[0];V0A;",100,-1,4000,100,-1,700);
  TH2D *hNAV0C = new TH2D("hNAV0C",";ZDCN-C[0];V0C;",100,-1,4000,100,-1,700);
  TH2D *hNAV0M = new TH2D("hNAV0M",";ZDCN-C[0];V0M;",100,-1,4000,100,-1,700);
  TH2D *hNAeta8 = new TH2D("hNAeta8",";ZDCN-C[0];RefMult08;",100,-1,4000,100,-1,700);
  TH2D *hPAV0A = new TH2D("hPAV0A",";ZDCP-A[0];V0A;",100,-1,4000,100,-1,700);
  TH2D *hPAV0C = new TH2D("hPAV0C",";ZDCP-A[0];V0C;",100,-1,4000,100,-1,700);
  TH2D *hPAV0M = new TH2D("hPAV0M",";ZDCP-A[0];V0M;",100,-1,4000,100,-1,700);
  TH2D *hPAeta8 = new TH2D("hPAeta8",";ZDCP-A[0];RefMult08;",100,-1,4000,100,-1,700);
  TH2D *totAV0M = new TH2D("totAV0M",";ZDCP-A[0]+ZDCN-A[0];V0M;",100,-1,4000,100,-1,700);
  TH2D *totAeta8 = new TH2D("totAeta8",";ZDCP-A[0]+ZDCN-A[0];RefMult08;",100,-1,4000,100,-1,700);
  TH2D *totAV0A = new TH2D("totAV0A",";ZDCP-A[0]+ZDCN-A[0];V0A;",100,-1,4000,100,-1,700);
  TH2D *totAV0C = new TH2D("totAV0C",";ZDCP-A[0]+ZDCN-A[0];V0C;",100,-1,4000,100,-1,700);
  //EZDC histos (Sum all)
  TH2D *ACV0M = new TH2D("ACV0M",";ZDCP+ZDCN;V0M;",100,-1,5000,100,-1,700);
  TH2D *ACeta8 = new TH2D("ACeta8",";ZDCP+ZDC;RefMult08;",100,-1,5000,100,-1,700);
  TH2D *ACV0A = new TH2D("ACV0A",";ZDCP+ZDCN;V0A;",100,-1,5000,100,-1,700);
  TH2D *ACV0C = new TH2D("ACV0C",";ZDCP+ZDCN;V0C;",100,-1,5000,100,-1,700);
  //Resolution histos
  TH2D * resP1 = new TH2D("resP1","Resolution P1;(sum+common)/2;sum-common;",bin,-1,3999,bin,-2000,2000);
  TH2D * resN1 = new TH2D("resN1","Resolution N1;(sum+common)/2;sum-common;",bin,-1,3999,bin,-2000,2000);
  TH2D * resP2 = new TH2D("resP2","Resolution P2;(sum+common)/2;sum-common;",bin,-1,3999,bin,-2000,2000);
  TH2D * resN2 = new TH2D("resN2","Resolution N2;(sum+common)/2;sum-common ;",bin,-1,3999,bin,-2000,2000);
  TH2D * hresP1 = new TH2D("hresP1","Resolution P1;(sum+common)/2;sum-common;",bin,-1,3999,bin,-100,100);
  TH2D * hresN1 = new TH2D("hresN1","Resolution N1;(sum+common)/2;sum-common;",bin,-1,3999,bin,-100,100);
  TH2D * hresP2 = new TH2D("hresP2","Resolution P2;(sum+common)/2;sum-common;",bin,-1,3999,bin,-100,100);
  TH2D * hresN2 = new TH2D("hresN2","Resolution N2;(sum+common)/2;sum-common ;",bin,-1,3999,bin,-100,100);
  //ZDCSum vs common histos
  TH2D * hZDCP1 = new TH2D("hZDCP1",";ZDCP1[0];ZDCP1Sum;",bin,-1,4000,bin,-1,4000);
  TH2D * hZDCP2 = new TH2D("hZDCP2",";ZDCP2[0];ZDCP2Sum;",bin,-1,4000,bin,-1,4000);
  TH2D * hZDCN1 = new TH2D("hZDCN1",";ZDCN1[0];ZDCN1Sum;",bin,-1,4000,bin,-1,4000);
  TH2D * hZDCN2 = new TH2D("hZDCN2",";ZDCN2[0];ZDCN2Sum;",bin,-1,4000,bin,-1,4000);

  TH1D * hEZDC = new TH1D("hEZDC","EZDC", bin ,-1,4000);
  TH1D * hEZDC_perc = new TH1D("hEZDC_perc","EZDC_perc", bin ,-1.,100.);
  TProfile* profileV0M =  new TProfile ("profileV0M","ZDC percentile;V0M;",100,-1,110);
  TProfile* profileRefMult08= new TProfile ("profileRefMult08","ZDC percentile;V0M;",100,-1,110);
  TProfile* profileV0M_E =  new TProfile ("profileV0M_E","ZDC percentile;V0M;",100,-1,701);
  TProfile* profileRefMult08_E= new TProfile ("profileRefMult08_E","ZDC percentile;V0M;",100,-1,201);
  TH2D *percentileV0M = new TH2D("percentileV0M",";ZDC percentile (%);V0M;",100,-1,110,100,-1,700);
  TH1D * hEZDClog = new TH1D("hEZDClog","EZDC", bin, 0.,4.);
  TH1D * hEZDC_perclog = new TH1D("hEZDC_perclog","EZDC_perc", 100 ,-1.,101.);
  //--------------------------------------------------------------

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

      EZDC = ZDCP1[0]+ZDCP2[0]+ZDCN1[0]+ZDCN2[0];
      hEZDC->Fill(EZDC);

      LogE = TMath::Log10(TMath::Abs(EZDC)+1);
      hEZDClog->Fill(LogE);
      
    } // Tree loop

  //Now get nsigma and normalize histos
  TF1* gaus = new TF1 ("gaus", "gaus", -100,100);

  TFile* Write = new TFile ("Res12fpass2.root", "recreate");
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
  
  hEZDC->Write();
  hEZDClog->Write();

  for (Int_t i =0; i<bin; i++)//loop sui bin x
    {
      sigmaP1[i] = resP1_2->GetBinContent(i+1);//prendo la sigma dei dati y su quel bin i di x
      sigmaP2[i] = resP2_2->GetBinContent(i+1);
      sigmaN1[i] = resN1_2->GetBinContent(i+1);
      sigmaN2[i] = resN2_2->GetBinContent(i+1);
    }

  TH1D * hCumu = new TH1D(*hEZDC);
  hCumu->SetName("Cumulative distribution");
  hCumu->Reset();
  TH1D * hCumulog = new TH1D(*hEZDClog);
  hCumulog->SetName("Cumulative distribution Elog");
  hCumulog->Reset();
  
  Double_t IntTOT = hEZDC->Integral(1,bin);// == total number of events
  Double_t val = 0;
  Double_t IntTOTlog = hEZDClog->Integral(1,bin); 
  Double_t vallog = 0;
 
  
  for (Int_t k=1; k<bin; k++) //loop through the bins
    {
      val+= (hEZDC->GetBinContent(k))/IntTOT;     
      hCumu->SetBinContent(k,val);
      vallog+= (hEZDClog->GetBinContent(k))/IntTOTlog;
      hCumulog->SetBinContent(k,vallog); 
    }
  
  hCumu->Write();
  hCumulog->Write();
  
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

      Int_t iN1 = hresN1->GetXaxis()->FindBin(meanN1);
      Int_t iP1 = hresP1->GetXaxis()->FindBin(meanP1);
      Int_t iN2 = hresN2->GetXaxis()->FindBin(meanN2);
      Int_t iP2 = hresP2->GetXaxis()->FindBin(meanP2);
      Double_t fN1 = dN1/sigmaN1[iN1-1];
      Double_t fP1 = dP1/sigmaP1[iP1-1];
      Double_t fN2 = dN2/sigmaN2[iN2-1];
      Double_t fP2 = dP2/sigmaP2[iP2-1];
      Float_t cutZDCP1,  cutZDCN1,  cutZDCP2,  cutZDCN2;

      if(sigmaN1[iN1-1]!=0) {	
	hresN1->Fill(meanN1,fN1);
	if (fN1<3 && fN1>(-3)) {
	  hZDCN1->Fill(ZDCN1[0],ZDCN1Sum);
	  if(ZDCP1[0]<100 ) {
	    hNCV0A->Fill(ZDCN1[0],V0AEstimator_abs);
	    hNCV0C->Fill(ZDCN1[0],V0CEstimator_abs);
	    hNCV0M->Fill(ZDCN1[0],V0MEstimator_abs);
	    hNCeta8->Fill(ZDCN1[0],RefMult08_abs);
	  }
	}
      }
      
      if(sigmaP1[iP1-1]!=0) {
	hresP1->Fill(meanP1,fP1);
	if (fP1<3 && fP1>(-3)) {
	  hZDCP1->Fill(ZDCP1[0],ZDCP1Sum);
	  if(ZDCN1[0]<100) {
	    hPCV0A->Fill(ZDCP1[0],V0AEstimator_abs);
	    hPCV0C->Fill(ZDCP1[0],V0CEstimator_abs);
	    hPCV0M->Fill(ZDCP1[0],V0MEstimator_abs);
	    hPCeta8->Fill(ZDCP1[0],RefMult08_abs);
	  }
	}
      }
      
      if(sigmaP2[iP2-1]!=0) {
	hresP2->Fill(meanP2,fP2);
	if (fP2<3 && fP2>(-3)) {
	  hZDCP2->Fill(ZDCP2[0],ZDCP2Sum);
	  if(ZDCN2[0]<100) {
	    hPAV0A->Fill(ZDCP2[0],V0AEstimator_abs);
	    hPAV0C->Fill(ZDCP2[0],V0CEstimator_abs);
	    hPAV0M->Fill(ZDCP2[0],V0MEstimator_abs);
	    hPAeta8->Fill(ZDCP2[0],RefMult08_abs);
	  }	  
	}
      }
      
      if(sigmaN2[iN2-1]!=0) {
	hresN2->Fill(meanN2,fN2);
	if (fN2<3 && fN2>(-3)) {
	  hZDCN2->Fill(ZDCN2[0],ZDCN2Sum);
	  if(ZDCP2[0]<100) {
	    hNAV0A->Fill(ZDCN2[0],V0AEstimator_abs);
	    hNAV0C->Fill(ZDCN2[0],V0CEstimator_abs);
	    hNAV0M->Fill(ZDCN2[0],V0MEstimator_abs);
	    hNAeta8->Fill(ZDCN2[0],RefMult08_abs);
	  }
	}
      }      
      
      /*if (sigmaP1[iP1]==0) {hresP1->Fill(meanP1,0.);};
	if (sigmaN1[iN1]==0) {hresN1->Fill(meanN1,0.);};
	if (sigmaP2[iP2]==0) {hresP2->Fill(meanP2,0.);};
	if (sigmaN2[iN2]==0) {hresN2->Fill(meanN2,0.);};*/
      /*
      if (fN2<3 && fN2>(-3) && sigmaN2[iN2-1]!=0) cutZDCN2 = ZDCN2[0];
      else cutZDCN2 = 0;
      
      if (fN1<3 && fN1>(-3) && sigmaN1[iN1-1]!=0) cutZDCN1 = ZDCN1[0];
      else cutZDCN1 = 0;
      
      if (fP2<3 && fP2>(-3) && sigmaP2[iP2-1]!=0) cutZDCP2 = ZDCP2[0];
      else cutZDCP2 = 0;
      
      if (fP1<3 && fP1>(-3) && sigmaP1[iP1-1]!=0) cutZDCP1 = ZDCP1[0];
      else cutZDCP1 = 0;
      */
      cutZDCN2=ZDCN2[0];
      cutZDCN1=ZDCN1[0];
      cutZDCP1=ZDCP1[0];
      cutZDCP2=ZDCP2[0];

      EZDC = ZDCP1[0]+ZDCP2[0]+ZDCN1[0]+ZDCN2[0];
      Double_t E = hCumu->Interpolate(EZDC);
      Double_t E100=E*100;

      LogE = TMath::Log10(TMath::Abs(EZDC)+1);
      Double_t ELog = hCumulog->Interpolate(LogE);
      Double_t ELog100=ELog*100;
      hEZDC_perclog->Fill(ELog100);
      
      hEZDC_perc->Fill(ELog100);
      profileV0M->Fill(ELog100,V0MEstimator_abs);
      profileRefMult08->Fill(ELog100,RefMult08_abs);
      profileV0M_E->Fill(V0MEstimator_abs,ELog100);
      profileRefMult08_E->Fill(RefMult08_abs,ELog100);
      percentileV0M->Fill(ELog100,V0MEstimator_abs);
      
      totCV0M->Fill(cutZDCN1+cutZDCP1,V0MEstimator_abs);
      totCeta8->Fill(cutZDCN1+cutZDCP1,RefMult08_abs);
      totCV0A->Fill(cutZDCN1+cutZDCP1,V0AEstimator_abs);
      totCV0C->Fill(cutZDCN1+cutZDCP1,V0CEstimator_abs);
      totAV0M->Fill(cutZDCN2+cutZDCP2,V0MEstimator_abs);
      totAeta8->Fill(cutZDCN2+cutZDCP2,RefMult08_abs);
      totAV0A->Fill(cutZDCN2+cutZDCP2,V0AEstimator_abs);
      totAV0C->Fill(cutZDCN2+cutZDCP2,V0CEstimator_abs);
      ACV0M->Fill(cutZDCN1+cutZDCP1+cutZDCN2+cutZDCP2,V0MEstimator_abs);
      ACeta8->Fill(cutZDCN1+cutZDCP1+cutZDCN2+cutZDCP2,RefMult08_abs);
      ACV0A->Fill(cutZDCN1+cutZDCP1+cutZDCN2+cutZDCP2,V0AEstimator_abs);
      ACV0C->Fill(cutZDCN1+cutZDCP1+cutZDCN2+cutZDCP2,V0CEstimator_abs);
	 
    } // Tree loop

  hEZDC_perc->Write();
  hEZDC_perclog->Write();
  percentileV0M->Write();
  percentileV0M->ProfileY()->Write();
  percentileV0M->ProfileX()->Write();
  
  //Write normalized histos  
  hresP1->Write();
  hresP1->FitSlicesY(gaus);
  TH1D * hresP1_1 = (TH1D*)gDirectory->Get("hresP1_1");
  hresP1_1->Write();
  TH1D * hresP1_2 = (TH1D*)gDirectory->Get("hresP1_2");
  TH1D * hminusresP1_2 = (TH1D*)gDirectory->Get("hresP1_2"); 
  hresP1_2->Scale(3);
  hresP1_2->SetLineColor(2);
  hresP1_2->Write();
  hresP2->Write();
  hresP2->FitSlicesY(gaus);
  TH1D * hresP2_1 = (TH1D*)gDirectory->Get("hresP2_1");
  hresP2_1->Write();
  TH1D * hresP2_2 = (TH1D*)gDirectory->Get("hresP2_2"); 
  hresP2_2->Scale(3);
  hresP2_2->SetLineColor(2);
  hresP2_2->Write();
  hresN1->Write();
  hresN1->FitSlicesY(gaus);
  TH1D * hresN1_1 = (TH1D*)gDirectory->Get("hresN1_1");
  hresN1_1->Write();
  TH1D * hresN1_2 = (TH1D*)gDirectory->Get("hresN1_2"); 
  hresN1_2->Scale(3);
  hresN1_2->SetLineColor(2);
  hresN1_2->Write();
  hresN2->Write();
  hresN2->FitSlicesY(gaus);
  TH1D * hresN2_1 = (TH1D*)gDirectory->Get("hresN2_1");
  hresN2_1->Write();
  TH1D * hresN2_2 = (TH1D*)gDirectory->Get("hresN2_2"); 
  hresN2_2->Scale(3);
  hresN2_2->SetLineColor(2);
  hresN2_2->Write();
  hZDCP1->Write();
  hZDCP2->Write();
  hZDCN1->Write();
  hZDCN2->Write();
  hNCV0A->Write();
  hNCV0A->ProfileY()->Write();
  hNCV0A->ProfileX()->Write();
  hNCV0C->Write();
  hNCV0C->ProfileY()->Write();
  hNCV0C->ProfileX()->Write();
  hNCV0M->Write();
  hNCV0M->ProfileY()->Write();
  hNCV0M->ProfileX()->Write();
  hNCeta8->Write();
  hNCeta8->ProfileY()->Write();
  hNCeta8->ProfileX()->Write();
  hPCV0A->Write();
  hPCV0A->ProfileY()->Write();
  hPCV0A->ProfileX()->Write();
  hPCV0C->Write();
  hPCV0C->ProfileY()->Write();
  hPCV0C->ProfileX()->Write();
  hPCV0M->Write();
  hPCV0M->ProfileY()->Write();
  hPCV0M->ProfileX()->Write();
  hPCeta8->Write();
  hPCeta8->ProfileY()->Write();
  hPCeta8->ProfileX()->Write();
  hNAV0A->Write();
  hNAV0A->ProfileY()->Write();
  hNAV0A->ProfileX()->Write();
  hNAV0C->Write();
  hNAV0C->ProfileY()->Write();
  hNAV0C->ProfileX()->Write();
  hNAV0M->Write();
  hNAV0M->ProfileY()->Write();
  hNAV0M->ProfileX()->Write();
  hNAeta8->Write();
  hNAeta8->ProfileY()->Write();
  hNAeta8->ProfileX()->Write();
  hPAV0A->Write();
  hPAV0A->ProfileY()->Write();
  hPAV0A->ProfileX()->Write();
  hPAV0C->Write();
  hPAV0C->ProfileY()->Write();
  hPAV0C->ProfileX()->Write();
  hPAV0M->Write();
  hPAV0M->ProfileY()->Write();
  hPAV0M->ProfileX()->Write();
  hPAeta8->Write();
  hPAeta8->ProfileY()->Write();
  hPAeta8->ProfileX()->Write();
  totCV0M->Write();
  totCV0M ->ProfileY()->Write();
  totCV0M->ProfileX()->Write();
  totCeta8->Write();
  totCeta8->ProfileY()->Write();
  totCeta8->ProfileX()->Write();
  totCV0A->Write();
  totCV0A ->ProfileY()->Write();
  totCV0A->ProfileX()->Write();
  totCV0C->Write();
  totCV0C->ProfileY()->Write();
  totCV0C->ProfileX()->Write();
  totAV0M->Write();
  totAV0M->ProfileY()->Write();
  totAV0M->ProfileX()->Write();
  totAeta8->Write();
  totAeta8->ProfileY()->Write();
  totAeta8->ProfileX()->Write();
  totAV0A->Write();
  totAV0A ->ProfileY()->Write();
  totAV0A->ProfileX()->Write();
  totAV0C->Write();
  totAV0C->ProfileY()->Write();
  totAV0C->ProfileX()->Write();
  ACV0M->Write();
  ACV0M->ProfileY()->Write();
  ACV0M->ProfileX()->Write();
  ACeta8->Write();
  ACeta8->ProfileY()->Write();
  ACeta8->ProfileX()->Write();
  ACV0A->Write();
  ACV0A->ProfileY()->Write();
  ACV0A->ProfileX()->Write();
  ACV0C->Write();
  ACV0C->ProfileY()->Write();
  ACV0C->ProfileX()->Write();

  profileV0M->Write();
  profileRefMult08->Write();
  profileV0M_E->Write();
  profileRefMult08_E->Write();
  
  Write->Close();
}
