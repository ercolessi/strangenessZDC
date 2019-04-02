
void AnalysisMult(Int_t Period=122)
{ 
  if (Period == 0) {
    cout<< "No Period specified...exiting" << "\n\n\n"
	<< " ex: LHC12b --> 122 (b second letter of the alphabet" << "\n\n\n" ;
    return;}
  
  Double_t cP1[5];
  Double_t cP2[5];
  Double_t cN1[5];
  Double_t cN2[5];

  TFile * Calib;

  if (Period==122)  Calib = new TFile ("PROVACalib12b_kINT7.root");
  if (Period==123)  Calib = new TFile ("ZDCCalib_12c.root");
  if (Period==126)  Calib = new TFile ("LHC12f.root");
  if (Period==189)  Calib = new TFile ("ZDCCalib_18i.root");

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

  Calib->Close();

  //-------------------------------
  TFile* Read;
  if (Period==122)  Read = new TFile ("LeadingMerged12b_kINT7.root");
  if (Period==123)  Read = new TFile ("LeadingMerged12c.root");
  if (Period==126)  Read = new TFile ("LeadingMerged12f.root");
  if (Period==189)  Read = new TFile ("LeadingMerged18i.root");

  TTree * T = (TTree *)Read->Get("tEvents");
  
  Int_t bin = 1000;
  Float_t ZDCP1[5], ZDCP1Sum, ZDCN1[5], ZDCN1Sum, ZDCP2[5], ZDCP2Sum, ZDCN2[5],
    ZDCN2Sum,dP1,meanP1, dN1, meanN1,dP2,meanP2, dN2, meanN2;
  Float_t V0M, V0A, V0C, CL0, CL1, SPDClusters, SPDTracklets, ZNA, ZNC, ZNApp,
    ZNCpp, RefMult08_p, RefMult05_p, V0MEstimator_p, V0CEstimator_p, V0AEstimator_p,
    RefMult08_abs, RefMult05_abs, V0MEstimator_abs, V0CEstimator_abs, V0AEstimator_abs;
  
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
  //C ZDCs
  TH2D *hNCV0A = new TH2D("hNCV0A",";ZDCN-C[0];V0A;",100,-1,4000,100,-1,700);
  TH2D *hNCV0C = new TH2D("hNCV0C",";ZDCN-C[0];V0C;",100,-1,4000,100,-1,700);
  TH2D *hNCV0M = new TH2D("hNCV0M",";ZDCN-C[0];V0M;",100,-1,4000,100,-1,700);
  TH2D *hNCeta8 = new TH2D("hNCeta8",";ZDCN-C[0];RefMult08;",100,-1,4000,100,-1,700);
  TH2D *hPCV0A = new TH2D("hPCV0A",";ZDCP-C[0];V0A;",100,-1,4000,100,-1,700);
  TH2D *hPCV0C = new TH2D("hPCV0C",";ZDCP-C[0];V0C;",100,-1,4000,100,-1,700);
  TH2D *hPCV0M = new TH2D("hPCV0M",";ZDCP-C[0];V0M;",100,-1,4000,100,-1,700);
  TH2D *hPCeta8 = new TH2D("hPCeta8",";ZDCP-C[0];RefMult08;",100,-1,4000,100,-1,700);
  //A ZDCs
  TH2D *hNAV0A = new TH2D("hNAV0A",";ZDCN-C[0];V0A;",100,-1,4000,100,-1,700);
  TH2D *hNAV0C = new TH2D("hNAV0C",";ZDCN-C[0];V0C;",100,-1,4000,100,-1,700);
  TH2D *hNAV0M = new TH2D("hNAV0M",";ZDCN-C[0];V0M;",100,-1,4000,100,-1,700);
  TH2D *hNAeta8 = new TH2D("hNAeta8",";ZDCN-C[0];RefMult08;",100,-1,4000,100,-1,700);
  TH2D *hPAV0A = new TH2D("hPAV0A",";ZDCP-C[0];V0A;",100,-1,4000,100,-1,700);
  TH2D *hPAV0C = new TH2D("hPAV0C",";ZDCP-C[0];V0C;",100,-1,4000,100,-1,700);
  TH2D *hPAV0M = new TH2D("hPAV0M",";ZDCP-C[0];V0M;",100,-1,4000,100,-1,700);
  TH2D *hPAeta8 = new TH2D("hPAeta8",";ZDCP-C[0];RefMult08;",100,-1,4000,100,-1,700);
  //A+C SUM
  TH2D *totCV0M = new TH2D("totCV0M",";ZDCP-C[0]+ZDCN-C[0];V0M;",100,-1,4000,100,-1,700);
  TH2D *totCeta8 = new TH2D("totCeta8",";ZDCP-C[0]+ZDCN-C[0];RefMult08;",100,-1,4000,100,-1,700);
  TH2D *totCV0A = new TH2D("totCV0A",";ZDCP-C[0]+ZDCN-C[0];V0A;",100,-1,4000,100,-1,700);
  TH2D *totCV0C = new TH2D("totCV0C",";ZDCP-C[0]+ZDCN-C[0];V0C;",100,-1,4000,100,-1,700);
  TH2D *totAV0M = new TH2D("totAV0M",";ZDCP-C[0]+ZDCN-C[0];V0M;",100,-1,4000,100,-1,700);
  TH2D *totAeta8 = new TH2D("totAeta8",";ZDCP-C[0]+ZDCN-C[0];RefMult08;",100,-1,4000,100,-1,700);
  TH2D *totAV0A = new TH2D("totAV0A",";ZDCP-C[0]+ZDCN-C[0];V0A;",100,-1,4000,100,-1,700);
  TH2D *totAV0C = new TH2D("totAV0C",";ZDCP-C[0]+ZDCN-C[0];V0C;",100,-1,4000,100,-1,700);
  TH2D *ACV0M = new TH2D("ACV0M",";ZDCP+ZDCN;V0M;",100,-1,5000,100,-1,700);
  TH2D *ACeta8 = new TH2D("ACeta8",";ZDCP+ZDC;RefMult08;",100,-1,5000,100,-1,700);
  TH2D *ACV0A = new TH2D("ACV0A",";ZDCP+ZDCN;V0A;",100,-1,5000,100,-1,700);
  TH2D *ACV0C = new TH2D("ACV0C",";ZDCP+ZDCN;V0C;",100,-1,5000,100,-1,700);
  
  TH2D * hZDCP1= new TH2D("hZDCP1",";ZDCP1[0];ZDCP1Sum;",bin,-1,4000,bin,-1,4000);
  TH2D * hZDCP2= new TH2D("hZDCP2",";ZDCP2[0];ZDCP2Sum;",bin,-1,4000,bin,-1,4000);
  TH2D * hZDCN1= new TH2D("hZDCN1",";ZDCN1[0];ZDCN1Sum;",bin,-1,4000,bin,-1,4000);
  TH2D * hZDCN2= new TH2D("hZDCN2",";ZDCN2[0];ZDCN2Sum;",bin,-1,4000,bin,-1,4000);
  TH2D * resP1= new TH2D("resP1","Resolution P1;sum-common;(sum+common)/2;",bin,-1,4000,bin,-2000,2000);
  TH2D * resN1= new TH2D("resN1","Resolution N1;sum-common;sum+comm /2;",bin,-1,4000,bin,-2000,2000);
  TH2D * resP2= new TH2D("resP2","Resolution P2;sum-common;(sum+common)/2;",bin,-1,4000,bin,-2000,2000);
  TH2D * resN2= new TH2D("resN2","Resolution N2;sum-common;sum+comm /2 ;",bin,-1,4000,bin,-2000,2000);
  TH2D * hZDCNP1= new TH2D("hZDCNP1","ZDCN1[0] vs ZDCP1[0];ZDCN1[0];ZDCP1[0];",bin,-1,2000,bin,-1,2000);
  TH2D * hZDCNP2= new TH2D("hZDCNP2","ZDCN2[0] vs ZDCP2[0];ZDCN2[0];ZDCP2[0];",bin,-1,2000,bin,-1,2000);
 
  //----------------------
   for(Int_t i=0; i<T->GetEntries();i++)//loop over Tree entries
     {
       T->GetEvent(i);//get Tree entries
       
       for (Int_t l=1; l<5; l++)
	 {
	   ZDCP1[l]/=cP1[l];
	   ZDCP2[l]/=cP2[l];
	   ZDCN1[l]/=cN1[l];
	   ZDCN2[l]/=cN2[l];
	 }
       ZDCP1Sum=ZDCP1[1]+ZDCP1[2]+ZDCP1[3]+ZDCP1[4];
       ZDCP2Sum=ZDCP2[1]+ZDCP2[2]+ZDCP2[3]+ZDCP2[4];
       ZDCN1Sum=ZDCN1[1]+ZDCN1[2]+ZDCN1[3]+ZDCN1[4];
       ZDCN2Sum=ZDCN2[1]+ZDCN2[2]+ZDCN2[3]+ZDCN2[4];
       
       hZDCNP1->Fill(ZDCN1[0],ZDCP1[0]);
       hZDCNP2->Fill(ZDCN2[0],ZDCP2[0]);
       
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

       hZDCP1->Fill(ZDCP1[0],ZDCP1Sum);
       hZDCP2->Fill(ZDCP2[0],ZDCP2Sum);
       hZDCN1->Fill(ZDCN1[0],ZDCN1Sum);
       hZDCN2->Fill(ZDCN2[0],ZDCN2Sum);
      
       //Pcut
       if(ZDCP1[0]<100 )// && ZDCN1[1]>150 && ZDCN1[2]>150 && ZDCN1[0]>150)
	 {
	   hNCV0A->Fill(ZDCN1[0],V0AEstimator_abs);
	   hNCV0C->Fill(ZDCN1[0],V0CEstimator_abs);
	   hNCV0M->Fill(ZDCN1[0],V0MEstimator_abs);
	   hNCeta8->Fill(ZDCN1[0],RefMult08_abs);
	 }
       //Ncut
       if(ZDCN1[0]<100 )//&& ZDCP1[1]>150 && ZDCP1[2]>150 && ZDCP1[0]>150)
	 {
	   hPCV0A->Fill(ZDCP1[0],V0AEstimator_abs);
	   hPCV0C->Fill(ZDCP1[0],V0CEstimator_abs);
	   hPCV0M->Fill(ZDCP1[0],V0MEstimator_abs);
	   hPCeta8->Fill(ZDCP1[0],RefMult08_abs);
	 }
       
       //A ZDCs
       //Pcut
       if(ZDCP2[0]<100)//  && ZDCN2[1]>100 && ZDCN2[2]>100 && ZDCN2[0]>100)
	 {
	   hNAV0A->Fill(ZDCN2[0],V0AEstimator_abs);
	   hNAV0C->Fill(ZDCN2[0],V0CEstimator_abs);
	   hNAV0M->Fill(ZDCN2[0],V0MEstimator_abs);
	   hNAeta8->Fill(ZDCN2[0],RefMult08_abs);
	 }

       //Ncut
       if(ZDCN2[0]<100)// && ZDCP2[1]>100 && ZDCP2[2]>100 && ZDCP2[0]>100)
	 {
	   hPAV0A->Fill(ZDCP2[0],V0AEstimator_abs);
	   hPAV0C->Fill(ZDCP2[0],V0CEstimator_abs);
	   hPAV0M->Fill(ZDCP2[0],V0MEstimator_abs);
	   hPAeta8->Fill(ZDCP2[0],RefMult08_abs);
	 }

       totCV0M->Fill(ZDCN1[0]+ZDCP1[0],V0MEstimator_abs);
       totCeta8->Fill(ZDCN1[0]+ZDCP1[0],RefMult08_abs);
       totCV0A->Fill(ZDCN1[0]+ZDCP1[0],V0AEstimator_abs);
       totCV0C->Fill(ZDCN1[0]+ZDCP1[0],V0CEstimator_abs);

       totAV0M->Fill(ZDCN2[0]+ZDCP2[0],V0MEstimator_abs);
       totAeta8->Fill(ZDCN2[0]+ZDCP2[0],RefMult08_abs);
       totAV0A->Fill(ZDCN2[0]+ZDCP2[0],V0AEstimator_abs);
       totAV0C->Fill(ZDCN2[0]+ZDCP2[0],V0CEstimator_abs);
       
       ACV0M->Fill(ZDCN1[0]+ZDCP1[0]+ZDCN2[0]+ZDCP2[0],V0MEstimator_abs);
       ACeta8->Fill(ZDCN1[0]+ZDCP1[0]+ZDCN2[0]+ZDCP2[0],RefMult08_abs);
       ACV0A->Fill(ZDCN1[0]+ZDCP1[0]+ZDCN2[0]+ZDCP2[0],V0AEstimator_abs);
       ACV0C->Fill(ZDCN1[0]+ZDCP1[0]+ZDCN2[0]+ZDCP2[0],V0CEstimator_abs);
     }

   //Write on root file
   TFile* Write;
   if (Period==122)  Write = new TFile ("PROVA12b_kINT7.root", "recreate");   
   if (Period==123)  Write = new TFile ("LHC12c.root", "recreate");   
   if (Period==126)  Write = new TFile ("CUTonZDCMLHC12f.root", "recreate");   
   if (Period==189)  Write = new TFile ("LHC18i.root", "recreate");
   gStyle->SetOptFit();

   hZDCNP1->Write();
   hZDCNP2->Write();
   resP1->Write();
   resN1->Write();
   resP2->Write();
   resN2->Write();
   hZDCP1->Write();
   hZDCP2->Write();
   hZDCN1->Write();
   hZDCN2->Write();
   
   hNCV0A->Write();
   hNCV0A->ProfileY()->Write();
   hNCV0A->ProfileX()->Write();
   //hNCV0A_pfy->Write();
   //hNCV0A_pfx->Write();
   
   hNCV0C->Write();
   hNCV0C->ProfileY()->Write();
   hNCV0C->ProfileX()->Write();
   //hNCV0C_pfy->Write();
   //hNCV0C_pfx->Write();
   
   hNCV0M->Write();
   hNCV0M->ProfileY()->Write();
   hNCV0M->ProfileX()->Write();
   //hNCV0M_pfy->Write();
   //hNCV0M_pfx->Write();
   
   hNCeta8->Write();
   hNCeta8->ProfileY()->Write();
   hNCeta8->ProfileX()->Write();
   //hNCeta8_pfy->Write();
   //hNCeta8_pfx->Write();
   
   hPCV0A->Write();
   hPCV0A->ProfileY()->Write();
   hPCV0A->ProfileX()->Write();
   //hPCV0A_pfy->Write();
   //hPCV0A_pfx->Write();
   
   hPCV0C->Write();
   hPCV0C->ProfileY()->Write();
   hPCV0C->ProfileX()->Write();
   //hPCV0C_pfy->Write();
   //hPCV0C_pfx->Write();
   
   hPCV0M->Write();
   hPCV0M->ProfileY()->Write();
   hPCV0M->ProfileX()->Write();
   //hPCV0M_pfy->Write();
   //hPCV0M_pfx->Write();
   
   hPCeta8->Write();
   hPCeta8->ProfileY()->Write();
   hPCeta8->ProfileX()->Write();
   //hPCeta8_pfy->Write();
   //hPCeta8_pfx->Write();

   hNAV0A->Write();
   hNAV0A->ProfileY()->Write();
   hNAV0A->ProfileX()->Write();
   //hNAV0A_pfy->Write();
   //hNAV0A_pfx->Write();
   
   hNAV0C->Write();
   hNAV0C->ProfileY()->Write();
   hNAV0C->ProfileX()->Write();
   //hNAV0C_pfy->Write();
   //hNAV0C_pfx->Write();
   
   hNAV0M->Write();
   hNAV0M->ProfileY()->Write();
   hNAV0M->ProfileX()->Write();
   //hNAV0M_pfy->Write();
   //hNAV0M_pfx->Write();
   
   hNAeta8->Write();
   hNAeta8->ProfileY()->Write();
   hNAeta8->ProfileX()->Write();
   //hNAeta8_pfy->Write();
   //hNAeta8_pfx->Write();
   
   hPAV0A->Write();
   hPAV0A->ProfileY()->Write();
   hPAV0A->ProfileX()->Write();
   //hPAV0A_pfy->Write();
   //hPAV0A_pfx->Write();
   
   hPAV0C->Write();
   hPAV0C->ProfileY()->Write();
   hPAV0C->ProfileX()->Write();
   //hPAV0C_pfy->Write();
   //hPAV0C_pfx->Write();
   
   hPAV0M->Write();
   hPAV0M->ProfileY()->Write();
   hPAV0M->ProfileX()->Write();
   //hPAV0M_pfy->Write();
   //hPAV0M_pfx->Write();
   
   hPAeta8->Write();
   hPAeta8->ProfileY()->Write();
   hPAeta8->ProfileX()->Write();
   //hPAeta8_pfy->Write();
   //hPAeta8_pfx->Write();
   
   totCV0M->Write();
   totCV0M ->ProfileY()->Write();
   totCV0M->ProfileX()->Write();
   //totCV0M_pfy->Write();
   //totCV0M_pfx->Write();
   
   totCeta8->Write();
   totCeta8->ProfileY()->Write();
   totCeta8->ProfileX()->Write();
   //totCeta8_pfy->Write();
   //totCeta8_pfx->Write();
   
   totCV0A->Write();
   totCV0A ->ProfileY()->Write();
   totCV0A->ProfileX()->Write();
   //totCV0A_pfy->Write();
   // totCV0A_pfx->Write();
   
   totCV0C->Write();
   totCV0C->ProfileY()->Write();
   totCV0C->ProfileX()->Write();
   //totCV0C_pfy->Write();
   //totCV0C_pfx->Write();

   totAV0M->Write();
   totAV0M ->ProfileY()->Write();
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
   //ACV0M_pfy->Write();
   //ACV0M_pfx->Write();
   
   ACeta8->Write();
   ACeta8->ProfileY()->Write();
   ACeta8->ProfileX()->Write();
   //ACeta8_pfy->Write();
   //ACeta8_pfx->Write();
   
   ACV0A->Write();
   ACV0A->ProfileY()->Write();
   ACV0A->ProfileX()->Write();
   //ACV0A_pfy->Write();
   //ACV0A_pfx->Write();
   
   ACV0C->Write();
   ACV0C->ProfileY()->Write();
   ACV0C->ProfileX()->Write();
   //ACV0C_pfy->Write();
   //ACV0C_pfx->Write();

   Write->Close();
}




 
