#include <iostream>
#include <fstream>
//#define max(a,b) (a>b ? a : b)

double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr );
TH1D* makeSystPlotsV0s(TString lV0Type = "Lambda", TString lWhichEstimator = "V0M", Double_t lMultBoundLo = 0.0, Double_t lMultBoundHi = 100.0, TString lWhichSystVar = "V0Radius");

void ComputeSystematics(TString Type = "Lambda", Double_t Low = 0.0, Double_t High = 100.0){
  
  TH1D* hRadius = makeSystPlotsV0s(Type,"V0M",Low,High,"V0Radius");
  TH1D* hNegToPV = makeSystPlotsV0s(Type,"V0M",Low,High,"DCANegToPV");
  TH1D* hPosToPV = makeSystPlotsV0s(Type,"V0M",Low,High,"DCAPosToPV");
  TH1D* hDaught = makeSystPlotsV0s(Type,"V0M",Low,High,"DCAV0Daughters");
  TH1D* hCosPA = makeSystPlotsV0s(Type,"V0M",Low,High,"V0CosPA");
  TH1D* hPLT = makeSystPlotsV0s(Type,"V0M",Low,High,"ProperLifetime");
  TH1D* hNOfCRows = makeSystPlotsV0s(Type,"V0M",Low,High,"NumberOfCrossedRows");
  TH1D* hNOfCRowsOverFind = makeSystPlotsV0s(Type,"V0M",Low,High,"NumberOfCrossedRowsOverFindable"); 
  TH1D* hTPCNSigmas = makeSystPlotsV0s(Type,"V0M",Low,High,"TPCPIDNSigmas");
  TH1D* hV0Rej = makeSystPlotsV0s(Type,"V0M",Low,High,"CompetingV0Rejection");
  TH1D* hSigmasForExtr = makeSystPlotsV0s(Type,"V0M",Low,High,"SigmaForSignalExtraction");

  //Other contributions:
  Double_t MaterialBudget = 0.04;
  Double_t MultIndipEffic = 0.02;
  Double_t MultIndipFeedDown = 0.02;

  TH1D* hSystTopological = (TH1D*)hRadius->Clone("hSystTopological");
  hSystTopological->Reset();
  Double_t Radius = 0;
  Double_t NegToPV= 0;
  Double_t PosToPV= 0;
  Double_t Daught= 0;
  Double_t CosPA= 0;

  for (int i = 1; i <= hRadius->GetNbinsX(); i++){
    if (TMath::Abs(hRadius->GetBinContent(i)/hRadius->GetBinError(i)) > 2) Radius = hRadius->GetBinContent(i);
    else Radius = 0;
    if (TMath::Abs(hNegToPV->GetBinContent(i)/hNegToPV->GetBinError(i)) > 2) NegToPV= hNegToPV->GetBinContent(i);
    else  NegToPV = 0;
    if (TMath::Abs(hPosToPV->GetBinContent(i)/hPosToPV->GetBinError(i)) > 2) PosToPV = hPosToPV->GetBinContent(i);
    else PosToPV = 0;
    if (TMath::Abs(hDaught->GetBinContent(i)/hDaught->GetBinError(i)) > 2) Daught = hDaught->GetBinContent(i);
    else Daught = 0;
    if (TMath::Abs(hCosPA->GetBinContent(i)/hCosPA->GetBinError(i)) > 2) CosPA = hCosPA->GetBinContent(i);
    else CosPA = 0;
    
    hSystTopological->SetBinContent(i, 
      TMath::Sqrt(Radius*Radius+NegToPV*NegToPV+PosToPV*PosToPV+Daught*Daught+CosPA*CosPA)    
    );
  }

  TH1D* hSystOthers = (TH1D*)hRadius->Clone("hSystOthers");
  hSystOthers->Reset();
  Double_t V0Rej = 0;
  Double_t PLT= 0;
  Double_t NOfCRows = 0;
  Double_t NOfCRowsOverFind = 0;
  Double_t TPCNSigmas = 0;
  Double_t SigmasForExtr = 0;

  for (int i = 1; i<= hRadius->GetNbinsX(); i++){
      if (TMath::Abs(hV0Rej->GetBinContent(i)/hV0Rej->GetBinError(i)) > 2) V0Rej = hV0Rej->GetBinContent(i);
      else V0Rej = 0;
      if (TMath::Abs(hPLT->GetBinContent(i)/hPLT->GetBinError(i)) > 2) PLT = hPLT->GetBinContent(i);
      else  PLT= 0;
      if (TMath::Abs(hNOfCRows->GetBinContent(i)/hNOfCRows->GetBinError(i)) > 2) NOfCRows = hNOfCRows->GetBinContent(i);
      else  NOfCRows= 0;
      if (TMath::Abs(hNOfCRowsOverFind->GetBinContent(i)/hNOfCRowsOverFind->GetBinError(i)) > 2) NOfCRowsOverFind = hNOfCRowsOverFind->GetBinContent(i);
      else NOfCRowsOverFind = 0;
      if (TMath::Abs(hTPCNSigmas->GetBinContent(i)/hTPCNSigmas->GetBinError(i)) > 2) TPCNSigmas = hTPCNSigmas->GetBinContent(i);
      else TPCNSigmas = 0;
      if (TMath::Abs(hSigmasForExtr->GetBinContent(i)/hSigmasForExtr->GetBinError(i)) > 2) SigmasForExtr  = hSigmasForExtr->GetBinContent(i);
      else SigmasForExtr = 0;
      cout << SigmasForExtr << "     " << TMath::Abs(hSigmasForExtr->GetBinContent(i)/hSigmasForExtr->GetBinError(i)) << endl;
   
      hSystOthers->SetBinContent(i, 
      TMath::Sqrt(PLT*PLT+NOfCRows*NOfCRows+NOfCRowsOverFind*NOfCRowsOverFind+TPCNSigmas*TPCNSigmas+SigmasForExtr*SigmasForExtr+V0Rej*V0Rej)
    );
  }

  TH1D* hSystTot = (TH1D*)hRadius->Clone("hSystTot");  
  hSystTot->Reset();
  for (int i = 1; i<= hRadius->GetNbinsX(); i++){
      hSystTot->SetBinContent(i, 
        TMath::Sqrt(hSystOthers->GetBinContent(i)*hSystOthers->GetBinContent(i) + 
          hSystTopological->GetBinContent(i)*hSystTopological->GetBinContent(i)) 
      );
  }

  TCanvas* cn = new TCanvas("cn","",950,900);
  TLegend* legend = new TLegend (0.14,0.7,0.6,0.9);
  legend->AddEntry(hSystTopological,"Topological systematics","L");
  legend->AddEntry(hSystOthers,"Selection cuts systematics","L");
  legend->AddEntry(hSystTot,"Total systematics","L");
  legend->SetTextSize(0.03);
  hSystTot->SetLineWidth(2);
  hSystTot->SetLineColor(kBlack);
  hSystOthers->SetLineWidth(2);
  hSystOthers->SetLineColor(kRed);
  hSystTopological->SetLineWidth(2);
  hSystTopological->SetLineColor(kBlue);
  hSystTot->Draw();
  hSystTot->SetYTitle("Systematics");
  hSystTot->GetYaxis()->SetRangeUser(-0.005,0.25);
  hSystTot->GetYaxis()->SetTitleOffset(1.);
  hSystTot->SetTitle("Systematics contributions for #Lambda");
  hSystOthers->Draw("SAME");
  hSystTopological->Draw("SAME");
  legend->Draw("SAME");
  cn->SaveAs(Form("%sSystematicsDisplayed_%03.0f_%03.0f.png",Type.Data(),Low, High));


  TCanvas* contr = new TCanvas("contr","",1300,900);
  TLegend* l = new TLegend (0.14,0.7,0.6,0.9);
  l->AddEntry(hRadius,"V0 Radius","L");
  l->AddEntry(hNegToPV,"DCA Neg To PV","L");
  l->AddEntry(hPosToPV,"DCA Pos To PV","L");
  l->AddEntry(hDaught,"DCA V0 Daughters","L");
  l->AddEntry(hCosPA,"V0 Cosine of PA","L");
  l->SetTextSize(0.03);
  hRadius->SetLineColor(kRed+1);
  hNegToPV->SetLineColor(kOrange+1);
  hPosToPV->SetLineColor(kSpring-1);
  hDaught->SetLineColor(kAzure-4);
  hCosPA->SetLineColor(kBlue+2);
  hRadius->SetMarkerStyle(1);
  hNegToPV->SetMarkerStyle(1);
  hPosToPV->SetMarkerStyle(1);
  hDaught->SetMarkerStyle(1);
  hCosPA->SetMarkerStyle(1);
  hRadius->SetLineWidth(2);
  hNegToPV->SetLineWidth(2);
  hPosToPV->SetLineWidth(2);
  hDaught->SetLineWidth(2);
  hCosPA->SetLineWidth(2);
  hRadius->SetTitle("Contribution of Topological Variables");
  hRadius->Draw("L");
  hNegToPV->Draw("SAME L");
  hPosToPV->Draw("SAME L");
  hDaught->Draw("SAME L");
  hCosPA->Draw("SAME L");
  l->Draw("SAME");
  contr->SaveAs(Form("%sSystematicsTopological_%03.0f_%03.0f.png",Type.Data(),Low,High));


  TCanvas* contr1 = new TCanvas("contr1","",1300,900);
  TLegend* l1 = new TLegend (0.14,0.7,0.7,0.9);
  l1->AddEntry(hPLT,"Proper Life Time","L");
  l1->AddEntry(hNOfCRows,"N. Of Crossed Rows","L");
  l1->AddEntry(hNOfCRowsOverFind,"N. Of Crossed Rows Over Findable","L");
  l1->AddEntry(hTPCNSigmas,"PID N. Sigmas","L");
  l1->AddEntry(hSigmasForExtr,"Sigmas For Signal Extraction","L");
  l1->AddEntry(hV0Rej,"Competing V0 Rejection","L");
  l1->SetTextSize(0.03);
  hPLT->SetLineColor(kRed+1);
  hNOfCRows->SetLineColor(kOrange+1);
  hNOfCRowsOverFind->SetLineColor(kSpring-1);
  hTPCNSigmas->SetLineColor(kAzure-4);
  hSigmasForExtr->SetLineColor(kBlue+2);
  hV0Rej->SetLineColor(kBlack);
  hPLT->SetMarkerStyle(1);
  hNOfCRows->SetMarkerStyle(1);
  hNOfCRowsOverFind->SetMarkerStyle(1);
  hTPCNSigmas->SetMarkerStyle(1);
  hSigmasForExtr->SetMarkerStyle(1);
  hV0Rej->SetMarkerStyle(1);
  hPLT->SetLineWidth(2);
  hNOfCRows->SetLineWidth(2);
  hNOfCRowsOverFind->SetLineWidth(2);
  hTPCNSigmas->SetLineWidth(2);
  hSigmasForExtr->SetLineWidth(2);
  hV0Rej->SetLineWidth(2);
  hPLT->SetTitle("Contribution of Selection Variables");
  hPLT->Draw("L");
  hNOfCRows->Draw("SAME L");
  hNOfCRowsOverFind->Draw("SAME L");
  hTPCNSigmas->Draw("SAME L");
  hSigmasForExtr->Draw("SAME L");
  hV0Rej->Draw("SAME L");
  l1->Draw("SAME");
  contr1->SaveAs(Form("%sSystematicsSelection%03.0f_%03.0f.png",Type.Data(),Low,High));

  //Spettri
  TString lDataFilename;
  lDataFilename = Form("Files_%03.0f-%03.0f/Results-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_000_100_UseMCRatioFD.root", Low, High, Type.Data(), Low, High);
  TFile* InputFile; 
  InputFile = new TFile(lDataFilename.Data());
  
  //Get histograms from files
  TH1D* lHistPt;
  lHistPt = (TH1D *) InputFile->Get(Form("fHistPt%s", Type.Data()));
  const int nbins = lHistPt->GetNbinsX();
  TH1D* hSystPt = (TH1D*)lHistPt->Clone("hSystPt");
  for (int i = 1; i<= nbins; i++){
    hSystPt->SetBinContent(i,lHistPt->GetBinContent(i));
    hSystPt->SetBinError(i, TMath::Sqrt(
          MaterialBudget*MaterialBudget +
          MultIndipEffic*MultIndipEffic +
          MultIndipFeedDown*MultIndipFeedDown +
          hSystTot->GetBinContent(i)*hSystTot->GetBinContent(i))
      *lHistPt->GetBinContent(i));
    hSystPt->SetMarkerStyle(1);
  }

  TCanvas* s = new TCanvas("s","",900,1300);
  s->SetLogy();
  hSystPt->SetStats(kFALSE);
  hSystPt->GetYaxis()->SetTitleOffset(1.);
  if (Low == 0. && High == 100.){  
    hSystPt->SetLineColor(kBlack); 
    hSystPt->SetFillStyle(0); 
    lHistPt->SetLineColor(kBlack);   
    lHistPt->SetMarkerColor(kBlack);
  } 
  if (Low == 0. && High == 10.){  
    hSystPt->SetLineColor(kRed+1); 
    hSystPt->SetFillStyle(0);
    lHistPt->SetLineColor(kRed+1);  
    lHistPt->SetMarkerColor(kRed+1); 
  }
  if (Low == 10. && High == 20.){  
    hSystPt->SetLineColor(kRed-4); 
    hSystPt->SetFillStyle(0);
    lHistPt->SetLineColor(kRed-4);  
    lHistPt->SetMarkerColor(kRed-4); 
  }
  if (Low == 20. && High == 30.){  
    hSystPt->SetLineColor(kOrange+7); 
    hSystPt->SetFillStyle(0);
    lHistPt->SetMarkerColor(kOrange+7);
    lHistPt->SetLineColor(kOrange+7);   
  }
  if (Low == 30. && High == 40.){  
    hSystPt->SetLineColor(kYellow); 
    hSystPt->SetFillStyle(0);
    lHistPt->SetMarkerColor(kYellow);
    lHistPt->SetLineColor(kYellow);   
  }
  if (Low == 40. && High == 50.){  
     hSystPt->SetLineColor(kGreen+1); 
    hSystPt->SetFillStyle(0);
    lHistPt->SetMarkerColor(kGreen+1);
     lHistPt->SetLineColor(kGreen+1);   
  }
   if (Low == 50. && High == 70.){  
     hSystPt->SetLineColor(kAzure+8); 
    hSystPt->SetFillStyle(0);
    lHistPt->SetMarkerColor(kAzure+8);
    lHistPt->SetLineColor(kAzure+8);   
  }
   if (Low == 70. && High == 100.){  
    hSystPt->SetLineColor(kBlue+3); 
    hSystPt->SetFillStyle(0);
    lHistPt->SetMarkerColor(kBlue+3); 
    lHistPt->SetLineColor(kBlue+3);   
  }


  lHistPt->SetLineWidth(1);
  lHistPt->SetMarkerStyle(20);
  lHistPt->SetMarkerSize(0.8); 
  hSystPt->Draw("E2");
  lHistPt->Draw("SAME");
  s->SaveAs(Form("%sSpettriSystematics_%03.0f_%03.0f.png",Type.Data(),Low,High));

  //Rapporti con Fiorella
  if ((Low == 0. && High == 100.)||(Low==40.&&High==50.)||(Low==70.&&High==100)){
    TFile* InputFileFiorella; 
    InputFileFiorella = new TFile(Form("~/Scaricati/resultsSpectra/Results-%s-13TeV-V0M_%.1f_%.1f_UseMCRatioFD_WithV0refitAndImprovedDCA.root",Type.Data(),Low,High),"READ");
    TH1D* hSpettroFior = (TH1D *) InputFileFiorella->Get(Form("fHistPt%s",Type.Data()));
    TH1D* hRapStat = (TH1D*) hSpettroFior->Clone("hRapStat");
    hRapStat->Reset();
    for (int i = 1; i < hSpettroFior->GetNbinsX(); i++){
      hRapStat->SetBinContent(i,lHistPt->GetBinContent(i)/hSpettroFior->GetBinContent(i));
      hRapStat->SetBinError(i,ErrorInRatio(lHistPt->GetBinContent(i),lHistPt->GetBinError(i),hSpettroFior->GetBinContent(i),hSpettroFior->GetBinError(i)));
    }
    
    TH1D* hRapSyst = (TH1D*) hSpettroFior->Clone("hRapSyst");
    hRapSyst->Reset();
    for (int i = 1; i < hSpettroFior->GetNbinsX(); i++){
      hRapSyst->SetBinContent(i,hSystPt->GetBinContent(i)/hSpettroFior->GetBinContent(i));
      hRapSyst->SetBinError(i,hSystPt->GetBinError(i)/hSpettroFior->GetBinContent(i));
        //ErrorInRatio(hSystPt->GetBinContent(i),hSystPt->GetBinError(i),hSpettroFior->GetBinContent(i),hSpettroFior->GetBinError(i)));
    }
    TCanvas* n = new TCanvas();
    n->SetGridy();
    TLegend* l2 = new TLegend (0.7,0.7,0.85,0.85);
    l2->SetTextSize(0.03);
    l2->AddEntry(hRapSyst,"syst","f");
    l2->AddEntry(hRapStat,"stat","LE");
    TLine* line = new TLine(0.8,1.,8.,1.);
    line->SetLineColor(kRed);
    line->SetLineWidth(2);
    line->SetLineStyle(7);
    hRapSyst->SetStats(kFALSE);
    hRapSyst->GetXaxis()->SetRangeUser(0.8,8.);
    hRapSyst->GetYaxis()->SetRangeUser(0.8,1.2);
    hRapSyst->SetYTitle("CorrSpectra_{THIS}/CorrSpectra_{FIORELLA}");
    hRapSyst->SetTitle(Form("V0M %.1f_%.1f",Low,High));
    hRapSyst->GetYaxis()->SetTitleOffset(1.3);
    hRapSyst->SetFillColor(16);
    hRapSyst->Draw("E2");
    hRapStat->SetLineWidth(2);
    hRapStat->SetLineColor(kBlack);
    hRapStat->SetMarkerStyle(20);
    hRapStat->SetMarkerSize(0.8);
    hRapStat->Draw("SAME");
    line->Draw("SAME");
    l2->Draw("SAME");
    n->SaveAs(Form("%sRapportiFiorella-V0M_%.1f_%.1f.png",Type.Data(),Low,High));
  }

  TFile* Write = new TFile (Form("SystematicsFinalResults-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_000_100.root",Type.Data(),Low, High), "recreate");
  lHistPt->Write();
  hSystPt->Write();
  hSystTopological->Write();
  hSystOthers->Write();
  hSystTot->Write(); 
  
  return;
}

//---------------------------------------------------------------
TH1D* makeSystPlotsV0s(
          TString lV0Type = "Lambda",
          TString lWhichEstimator = "V0M",
          Double_t lMultBoundLo = 0.0,
          Double_t lMultBoundHi = 100.0,
          TString lWhichSystVar = "V0Radius")
{
  cout<<"  "<<endl;
  cout<<"-------------------------------------------------------------"<<endl;
  cout<<"                  Make systematics for V0s                   "<<endl;
  cout<<"-------------------------------------------------------------"<<endl;
  cout<<endl;
  

  //Definition of variables
  //Number of variation cuts
  int n = 0;
  if (lWhichSystVar == "NumberOfCrossedRowsOverFindable" || lWhichSystVar == "CompetingV0Rejection")  n = 2;
  else if (lWhichSystVar == "NumberOfCrossedRows" || lWhichSystVar == "ProperLifetime")  n = 3;
  else if (lWhichSystVar == "TPCPIDNSigmas" || lWhichSystVar == "SigmaForSignalExtraction") n = 4;
  else n = 5;
  const int nfiles = n;  
  TH1D* hMaxDev;
  TH1D* lHistPt[nfiles];
  TH1D* lHistPtRaw[nfiles];
  TH1D* hCut[nfiles-1];
  TString lDataFilename[nfiles];
  TFile* InputFile[nfiles]; 
  

  //Set data files
  TString lSystFile = Form("Files_%03.0f-%03.0f/Results-Systematics-V0M_%03.0f_%03.0f_ZDC_000_100_UseMCRatioFD",lMultBoundLo,lMultBoundHi,lMultBoundLo,lMultBoundHi); 
  if(lV0Type == "Lambda"    ) lSystFile.Append("-Lambda-");
  if(lV0Type == "AntiLambda") lSystFile.Append("-AntiLambda-");
  
  //Files
  lDataFilename[0] = Form("Files_%03.0f-%03.0f/Results-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_000_100_UseMCRatioFD.root",lMultBoundLo,lMultBoundHi,lV0Type.Data(),lMultBoundLo,lMultBoundHi);
  InputFile[0] = new TFile(lDataFilename[0].Data(),"READ");
  for (int i = 1; i <= nfiles-1; i++){
    lDataFilename[i] = lSystFile + lWhichSystVar + Form("-%i.root",i);
    InputFile[i] = new TFile(lDataFilename[i].Data(),"READ");
  }  
  cout<<endl<<"\n\n---> Do Systematics for: "<< lWhichSystVar << endl;
  cout<<endl<<"\n\n---> Number of syst files: "<< nfiles-1 << endl;
  cout << " \n\n--- Set Minimum Bias file   :  " << lDataFilename[0].Data() << endl;
  for (int i = 1; i <= nfiles-1; i++){
    cout << " --- Set Systematics file #"<< i << " :  " << lDataFilename[i].Data() << endl;
  }  


  //Inizialize Histos
  lHistPt[0] = (TH1D *) InputFile[0]->Get(Form("fHistPt%s", lV0Type.Data()));
  lHistPtRaw[0] = (TH1D *) InputFile[0]->Get("lInvMassReal/lInvMassRealRawData/fHistPtRaw");
  hMaxDev = (TH1D*)lHistPt[0]->Clone("hMaxDev");
  hMaxDev->Reset();
  for (int i = 1; i < nfiles; i++){
    lHistPt[i] = (TH1D *) InputFile[i]->Get(Form("fHistPt%s", lV0Type.Data()));
    lHistPtRaw[i] = (TH1D *) InputFile[i]->Get("lInvMassReal/lInvMassRealRawData/fHistPtRaw");      
    hCut[i-1] = (TH1D*)lHistPt[0]->Clone(Form("hCut%i",i));
    hCut[i-1]->Reset();
  }
    

  //Fill histos for variation cuts
  for (int k = 1; k < nfiles; k++){
    for (int i = 1; i<= lHistPt[0]->GetNbinsX(); i++){
        hCut[k-1]->SetBinContent(i,lHistPt[k]->GetBinContent(i)/lHistPt[0]->GetBinContent(i));
        hCut[k-1]->SetBinError(i,ErrorInRatio(lHistPt[k]->GetBinContent(i),lHistPt[k]->GetBinError(i),lHistPt[0]->GetBinContent(i),lHistPt[0]->GetBinError(i)));
     }
  }
    

  double binvalue[nfiles-1];
  double maxvalue = 0.;
  for (int i = 1; i<= lHistPt[0]->GetNbinsX(); i++){
    for (int k = 1; k < nfiles; k++){ 
      binvalue[k-1] = TMath::Abs(hCut[k-1]->GetBinContent(i)-1);
    }
    maxvalue = binvalue[0];
    int counter = 0;
    for (int k = 1; k < nfiles; k++){ 
      maxvalue = max(maxvalue,binvalue[k]);
      if (maxvalue == binvalue[k]) counter = k;
    }       
    hMaxDev->SetBinError(i,hCut[counter]->GetBinError(i));
    hMaxDev->SetBinContent(i,TMath::Abs(maxvalue));
    if (TMath::Abs(maxvalue) == 0)hMaxDev->SetBinContent(i,0.00000001);
  }
    

  //Prepare Canvas
  //Max Deviation
  hMaxDev->GetXaxis()->SetRangeUser(0.7,8.);
  hMaxDev->GetYaxis()->SetRangeUser(-0.005,0.1);
  hMaxDev->SetYTitle("max rel. dev.");
  hMaxDev->SetTitle(Form("%s",lWhichSystVar.Data()));
  hMaxDev->GetYaxis()->SetTitleSize(0.05);
  hMaxDev->GetYaxis()->SetTitleOffset(1.0);
  hMaxDev->GetXaxis()->SetTitleSize(0.05);
  hMaxDev->GetXaxis()->SetTitleOffset(0.8);
  hMaxDev->SetMarkerStyle(20);
  hMaxDev->SetMarkerSize(1.1);
  hMaxDev->SetMarkerColor(kBlack);
  hMaxDev->SetLineColor(kBlack);
  TCanvas* maxdev = new TCanvas("maxdev"," ",1000,800);
  maxdev->SetGridy();
  maxdev->SetGridx();
  hMaxDev->SetStats(kFALSE);  
  hMaxDev->Draw();
  maxdev->SaveAs(Form("MaxRelDev%s_%03.0f_%03.0f.png",lWhichSystVar.Data(),lMultBoundLo,lMultBoundHi));

  //Cut Variation
  TCanvas* cutvar = new TCanvas("cutvar"," ",1000,800);
  cutvar->SetGridy();
  TLegend* legend = new TLegend (0.6,0.75,0.8,0.9);  
  if (lWhichSystVar == "NumberOfCrossedRowsOverFindable"){
    hCut[0]->SetMarkerColor(kGreen+3);
    hCut[0]->SetLineColor(kGreen+3);
    legend->AddEntry(hCut[0],"tight","LEP");
  }
  else if (lWhichSystVar == "CompetingV0Rejection"){
    hCut[0]->SetMarkerColor(kBlue+1);
    hCut[0]->SetLineColor(kBlue+1);
    legend->AddEntry(hCut[0],"loose","LEP");
  }

  else if (lWhichSystVar == "NumberOfCrossedRows"){
      hCut[0]->SetMarkerColor(kGreen+3);
      hCut[1]->SetMarkerColor(kAzure+9);
      hCut[0]->SetLineColor(kGreen+3);
      hCut[1]->SetLineColor(kAzure+9);
      legend->AddEntry(hCut[0],"tight","LEP");
      legend->AddEntry(hCut[1],"very tight","LEP");
    }
  else if (lWhichSystVar == "ProperLifetime"){
      hCut[0]->SetMarkerColor(kBlue+1);
      hCut[1]->SetMarkerColor(kGreen+3);
      hCut[0]->SetLineColor(kBlue+1);
      hCut[1]->SetLineColor(kGreen+3);
      legend->AddEntry(hCut[0],"loose","LEP");
      legend->AddEntry(hCut[1],"tight","LEP");
  }
  else if (lWhichSystVar == "TPCPIDNSigmas"){
      hCut[0]->SetMarkerColor(kRed);
      hCut[1]->SetMarkerColor(kBlue+1);
      hCut[2]->SetMarkerColor(kGreen+3);
      hCut[0]->SetLineColor(kRed);
      hCut[1]->SetLineColor(kBlue+1);
      hCut[2]->SetLineColor(kGreen+3);
      legend->AddEntry(hCut[0],"very loose","LEP");
      legend->AddEntry(hCut[1],"loose","LEP");
      legend->AddEntry(hCut[2],"tight","LEP");
  }
  else if (lWhichSystVar == "SigmaForSignalExtraction"){
      hCut[2]->SetMarkerColor(kAzure+9);
      hCut[0]->SetMarkerColor(kBlue+1);
      hCut[1]->SetMarkerColor(kGreen+3);
      hCut[2]->SetLineColor(kAzure+9);
      hCut[0]->SetLineColor(kBlue+1);
      hCut[1]->SetLineColor(kGreen+3);
      legend->AddEntry(hCut[0],"loose","LEP");
      legend->AddEntry(hCut[1],"tight","LEP");
      legend->AddEntry(hCut[2],"very tight","LEP");
  }  
  else {
    hCut[0]->SetMarkerColor(kRed);
    hCut[1]->SetMarkerColor(kBlue+1);
    hCut[2]->SetMarkerColor(kGreen+3);
    hCut[3]->SetMarkerColor(kAzure+9);
    hCut[0]->SetLineColor(kRed);
    hCut[1]->SetLineColor(kBlue+1);
    hCut[2]->SetLineColor(kGreen+3);
    hCut[3]->SetLineColor(kAzure+9);    
    legend->AddEntry(hCut[0],"very loose","LEP");
    legend->AddEntry(hCut[1],"loose","LEP");
    legend->AddEntry(hCut[2],"tight","LEP");
    legend->AddEntry(hCut[3],"very tight","LEP");
  }

  for (int k = 0; k < nfiles-1; k++){
    hCut[k]->GetXaxis()->SetRangeUser(0.7,8.);
    hCut[k]->GetYaxis()->SetRangeUser(0.9,1.1);
    hCut[k]->SetYTitle("Yield^{syst-cut} / Yield^{def-cut}");
    hCut[k]->SetTitle(Form("%s",lWhichSystVar.Data()));
    hCut[k]->GetYaxis()->SetTitleSize(0.05);
    hCut[k]->GetYaxis()->SetTitleOffset(0.9);
    hCut[k]->GetXaxis()->SetTitleSize(0.05);
    hCut[k]->GetXaxis()->SetTitleOffset(0.8);
    hCut[k]->SetMarkerStyle(20);
    hCut[k]->SetMarkerSize(1.1);
  }
  hCut[0]->Draw();
  hCut[0]->SetStats(kFALSE);
  for (int k = 1; k < nfiles-1; k++){
  hCut[k]->Draw("SAME");
  }
  legend->SetTextSize(0.035);
  legend->Draw("SAME");    
  cutvar->SaveAs(Form("%s_%03.0f_%03.0f.png",lWhichSystVar.Data(),lMultBoundLo,lMultBoundHi));


  //Return Max Dev Histo
  return hMaxDev;
}


//---------------------------------------------------------------------------------------------------
double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr ) {
    //Error in a Ratio
    if(B!=0) {
        Double_t errorfromtop = Aerr*Aerr / (B*B) ;
        Double_t errorfrombottom = ((A*A)/(B*B*B*B)) * Berr * Berr;
        return TMath::Sqrt( TMath::Abs(errorfromtop - errorfrombottom) );
    }
    return 1.;
}