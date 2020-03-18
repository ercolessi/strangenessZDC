#include <iostream>
#include <fstream>
//#define max(a,b) (a>b ? a : b)

double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr );
TH1D* DoRelativeUncertainty(TH1D* hMaxDev = 0x0, TString lWhichSystVar = "V0Radius",TString lV0Type = "Lambda");
TH1D* makeSystPlotsV0s(TString lV0Type = "Lambda", TString lWhichEstimator = "V0M", Double_t lMultBoundLo = 0.0, Double_t lMultBoundHi = 100.0, TString lWhichSystVar = "V0Radius");

void PlotV0s(TString Type = "Lambda", Double_t Low = 0.0, Double_t High = 100.0){
  
  TH1D* hRadius = makeSystPlotsV0s(Type,"V0M",Low,High,"V0Radius");
  TH1D* hNegToPV = makeSystPlotsV0s(Type,"V0M",Low,High,"DCANegToPV");
  TH1D* hPosToPV = makeSystPlotsV0s(Type,"V0M",Low,High,"DCAPosToPV");
  TH1D* hDaught = makeSystPlotsV0s(Type,"V0M",Low,High,"DCAV0Daughters");
  TH1D* hCosPA = makeSystPlotsV0s(Type,"V0M",Low,High,"V0CosPA");
  TH1D* hPLT = makeSystPlotsV0s(Type,"V0M",Low,High,"ProperLifetime");
  TH1D* hNOfCRows = makeSystPlotsV0s(Type,"V0M",Low,High,"NumberOfCrossedRows");
  //TH1D* hNOfCRowsOverFind = makeSystPlotsV0s(Type,"V0M",Low,High,"NumberOfCrossedRowsOverFindable"); 
  TH1D* hTPCNSigmas = makeSystPlotsV0s(Type,"V0M",Low,High,"TPCPIDNSigmas");
  //TH1D* hV0Rej = makeSystPlotsV0s(Type,"V0M",Low,High,"CompetingV0Rejection");
  TH1D* hSigmasForExtr = makeSystPlotsV0s(Type,"V0M",Low,High,"SigmaForSignalExtraction");

  TH1D* hSystTopological = (TH1D*)hRadius->Clone("hSystTopological");
  hSystTopological->Reset();
  for (int i = 1; i<= hRadius->GetNbinsX(); i++){
    hSystTopological->SetBinContent(i, 
      hRadius->GetBinContent(i)  + 
      hNegToPV->GetBinContent(i) +
      hPosToPV->GetBinContent(i) +
      hDaught->GetBinContent(i)  +
      hCosPA->GetBinContent(i)
      );
  }
  TH1D* hSystOthers = (TH1D*)hRadius->Clone("hSystOthers");
  hSystOthers->Reset();
  for (int i = 1; i<= hRadius->GetNbinsX(); i++){
      hSystOthers->SetBinContent(i, 
      hPLT->GetBinContent(i)  + 
      hNOfCRows->GetBinContent(i) +
      //hNOfCRowsOverFind->GetBinContent(i) +
      hTPCNSigmas->GetBinContent(i)  +
      hSigmasForExtr->GetBinContent(i)  
      //hV0Rej->GetBinContent(i)
      );
  }
  
  new TCanvas;
  hSystOthers->SetLineWidth(2);
  hSystOthers->SetLineColor(kRed);
  hSystTopological->SetLineWidth(2);
  hSystTopological->SetLineColor(kBlue);
  hSystTopological->Draw();
  hSystOthers->Draw("SAME");

  return;
}

//--------------------------------------------------------------
TH1D* DoRelativeUncertainty(TH1D* hMaxDev = 0x0, TString lWhichSystVar = "V0Radius", TString lV0Type = "Lambda")
{

  cout<<"  "<<endl;
  cout<<"-------------------------------------------------------------"<<endl;
  cout<<"        Make plots of relative uncertainty for V0s       "<<endl;
  cout<<"-------------------------------------------------------------"<<endl;
  cout<<endl;

  TH1D* hRelUCut = (TH1D*)hMaxDev->Clone("hRelUCut");
  hRelUCut->Reset();
  TString lDataFilename = Form("Results-%s-13TeV-V0M_000_100_ZDC_000_100_UseMCRatioFD.root",lV0Type.Data());
  TFile* InputFile = new TFile(lDataFilename.Data());  
  //Get histograms from files
  TH1D* lHistPt = (TH1D *) InputFile->Get(Form("fHistPt%s", lV0Type.Data()));

  for (int i = 1; i<= lHistPt->GetNbinsX(); i++){
    hRelUCut->SetBinContent(i, hMaxDev->GetBinContent(i)/lHistPt->GetBinContent(i)); 
    //Rel Uncertainty is Syst Uncer./Def. Corr. Spectra
    cout << "pT = " << lHistPt->GetBinCenter(i) << "   " << hMaxDev->GetBinContent(i) << "   " << lHistPt->GetBinContent(i) << endl;
  }
  hRelUCut->SetTitle(Form("Relative Uncertainty for %s", lWhichSystVar.Data()));
  hRelUCut->SetYTitle("Relative Uncertainty");
  hRelUCut->SetXTitle("p_{T}");
  return hRelUCut;
}

//---------------------------------------------------------------
TH1D* makeSystPlotsV0s(TString lV0Type = "Lambda",
		      TString lWhichEstimator = "V0M",
		      Double_t lMultBoundLo = 0.0,
          Double_t lMultBoundHi = 100.0,
		      TString lWhichSystVar = "V0Radius"){
  
  cout<<"  "<<endl;
  cout<<"-------------------------------------------------------------"<<endl;
  cout<<"        Macro to make plots of systematics for V0s       "<<endl;
  cout<<"-------------------------------------------------------------"<<endl;
  cout<<endl;
  
  //Set data files
  TString lSystFile = "Results-Systematics-V0M_000_100_ZDC_000_100_UseMCRatioFD"; 
  if(lV0Type == "Lambda"    ) lSystFile.Append("-Lambda-");
  if(lV0Type == "AntiLambda") lSystFile.Append("-AntiLambda-");

  int n = 0;
  if (lWhichSystVar == "NumberOfCrossedRowsOverFindable" || lWhichSystVar == "CompetingV0Rejection")  n = 2;
  else if (lWhichSystVar == "NumberOfCrossedRows" || lWhichSystVar == "ProperLifetime")  n = 3;
  else if (lWhichSystVar == "TPCPIDNSigmas" || lWhichSystVar == "SigmaForSignalExtraction") n = 4;
  else n = 5;
  TH1D* hMaxDev = 0x0;

  const int nfiles = n;
  
  TString lDataFilename[nfiles];
  lDataFilename[0] = Form("Results-%s-13TeV-V0M_000_100_ZDC_000_100_UseMCRatioFD.root",lV0Type.Data());
  for (int i = 1; i <= nfiles-1; i++){
    lDataFilename[i] = lSystFile + lWhichSystVar + Form("-%i.root",i);
  }
  
  cout<<endl<<"---> Do Systematics for: "<< lWhichSystVar << endl;
  cout<<endl<<"---> Number of syst files: "<< nfiles-1 << endl;
  
  
  cout << " --- Set Minimum Bias file   :  " << lDataFilename[0].Data() << endl;
  for (int i = 1; i <= nfiles-1; i++){
    cout << " --- Set Systematics file #"<< i << " :  " << lDataFilename[i].Data() << endl;
  }  
  
  TFile* InputFile[nfiles]; 
  for (int i = 0; i < nfiles; i++){
    InputFile[i] = new TFile(lDataFilename[i].Data());
  }
  
  //Get histograms from files
  TH1D* lHistPt[nfiles];
  for (int i = 0; i <= nfiles-1; i++){
    lHistPt[i] = (TH1D *) InputFile[i]->Get(Form("fHistPt%s", lV0Type.Data()));
  }
  
  //Get histograms from files
  TH1D* lHistPtRaw[nfiles];
  for (int i = 0; i <= nfiles-1; i++){
    lHistPtRaw[i] = (TH1D *) InputFile[i]->Get("lInvMassReal/lInvMassRealRawData/fHistPtRaw");
  }

   //################# TPC CR Findable and dE/dx Selection ##########################
  if (lWhichSystVar == "NumberOfCrossedRowsOverFindable" || lWhichSystVar == "CompetingV0Rejection") {

    TH1D* hCut1 = 0x0;

    if (lWhichSystVar == "NumberOfCrossedRowsOverFindable"){
      lHistPt[0]->SetName("MB");
      lHistPt[1]->SetName("tight");
      
      hCut1= (TH1D*)lHistPt[0]->Clone("tight");
    }

    else if (lWhichSystVar == "CompetingV0Rejection"){
      lHistPt[0]->SetName("MB");
      lHistPt[1]->SetName("loose");
     
      hCut1 = (TH1D*)lHistPt[0]->Clone("loose");
    }
        
    for (int i = 1; i <= lHistPt[0]->GetNbinsX(); i++){
      hCut1->SetBinContent(i,lHistPt[1]->GetBinContent(i)/lHistPt[0]->GetBinContent(i));
      hCut1->SetBinError(i,ErrorInRatio(lHistPt[1]->GetBinContent(i),lHistPt[1]->GetBinError(i),lHistPt[0]->GetBinContent(i),lHistPt[0]->GetBinError(i)));
    }

    hMaxDev = (TH1D*)lHistPt[0]->Clone("hMaxDev");
    hMaxDev->Reset();
    for (int i = 1; i <= lHistPt[0]->GetNbinsX(); i++){

      double maximo = TMath::Abs(hCut1->GetBinContent(i)-1.);
      double err = hCut1->GetBinError(i);
      cout << maximo << endl;
      hMaxDev->SetBinContent(i,TMath::Abs(maximo));
      if (TMath::Abs(maximo)==0)hMaxDev->SetBinContent(i,0.00000001);
      hMaxDev->SetBinError(i,err);
    }

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
    hMaxDev->SaveAs(Form("MaxRelDev%s.png",lWhichSystVar.Data()));
    TCanvas* c = new TCanvas("c"," ",1000,800);
    c->SetGridy();
    c->SetGridx();
    hMaxDev->SetStats(kFALSE);  
    hMaxDev->Draw();
    c->SaveAs(Form("MaxRelDev%s.png",lWhichSystVar.Data()));

    
    TCanvas* c1 = new TCanvas("c1"," ",1000,800);
    c1->SetGridy();
    hCut1->GetXaxis()->SetRangeUser(0.7,10.);
    hCut1->GetYaxis()->SetRangeUser(0.9,1.1);
    hCut1->SetYTitle("Yield^{syst-cut} / Yield^{def-cut}");
    hCut1->SetTitle(Form("%s",lWhichSystVar.Data()));
    hCut1->GetYaxis()->SetTitleSize(0.05);
    hCut1->GetYaxis()->SetTitleOffset(0.9);
    hCut1->GetXaxis()->SetTitleSize(0.05);
    hCut1->GetXaxis()->SetTitleOffset(0.8);
    hCut1->SetMarkerStyle(20);
    hCut1->SetMarkerSize(1.1);
    if (lWhichSystVar == "NumberOfCrossedRowsOverFindable"){
      hCut1->SetMarkerColor(kGreen+3);
      hCut1->SetLineColor(kGreen+3);
    }
    else if (lWhichSystVar == "CompetingV0Rejection"){
      hCut1->SetMarkerColor(kBlue+1);
      hCut1->SetLineColor(kBlue+1);
    }
    hCut1->Draw();
    TLegend* legend = new TLegend (0.6,0.8,0.8,0.9);
    if (lWhichSystVar == "NumberOfCrossedRowsOverFindable"){
      legend->AddEntry(hCut1,"tight","LEP");
    }
    else if (lWhichSystVar == "CompetingV0Rejection"){
      legend->AddEntry(hCut1,"loose","LEP");
    }
    legend->SetTextSize(0.05);
    legend->Draw("SAME");
    hCut1->SetStats(kFALSE);    
    c1->SaveAs(Form("%s.png",lWhichSystVar.Data())); 
    delete c;
    delete c1;
  }
    
  //################# TPC Crossed Rows and ProperLifetime ##########################
  if (lWhichSystVar == "NumberOfCrossedRows" || lWhichSystVar == "ProperLifetime") {

    TH1D* hCut1 = 0x0;
    TH1D* hCut2 = 0x0;
    
    if (lWhichSystVar == "NumberOfCrossedRows"){
      lHistPt[0]->SetName("MB");
      lHistPt[1]->SetName("loose");
      lHistPt[2]->SetName("tight");

      hCut1 = (TH1D*)lHistPt[0]->Clone("loose");
      hCut2 = (TH1D*)lHistPt[0]->Clone("tight");
    }

    else if (lWhichSystVar == "ProperLifetime"){
      lHistPt[0]->SetName("MB");
      lHistPt[1]->SetName("loose");
      lHistPt[2]->SetName("tight");

      hCut1 = (TH1D*)lHistPt[0]->Clone("loose");
      hCut2 = (TH1D*)lHistPt[0]->Clone("tight");
    }
    
    for (int i = 1; i<= lHistPt[0]->GetNbinsX(); i++){
      hCut1->SetBinContent(i,lHistPt[1]->GetBinContent(i)/lHistPt[0]->GetBinContent(i));
      hCut1->SetBinError(i,ErrorInRatio(lHistPt[1]->GetBinContent(i),lHistPt[1]->GetBinError(i),lHistPt[0]->GetBinContent(i),lHistPt[0]->GetBinError(i)));
      hCut2->SetBinContent(i,lHistPt[2]->GetBinContent(i)/lHistPt[0]->GetBinContent(i));
      hCut2->SetBinError(i,ErrorInRatio(lHistPt[2]->GetBinContent(i),lHistPt[2]->GetBinError(i),lHistPt[0]->GetBinContent(i),lHistPt[0]->GetBinError(i)));
    }

    hMaxDev = (TH1D*)lHistPt[0]->Clone("hMaxDev");
    hMaxDev->Reset();
    for (int i = 1; i<= lHistPt[0]->GetNbinsX(); i++){
      double maximo = max(TMath::Abs(hCut1->GetBinContent(i)-1),TMath::Abs(hCut2->GetBinContent(i)-1));
      double err = 0.;
      if (maximo == TMath::Abs(hCut1->GetBinContent(i)-1)) err = hCut1->GetBinError(i);
      if (maximo == TMath::Abs(hCut2->GetBinContent(i)-1)) err = hCut2->GetBinError(i);
      hMaxDev->SetBinContent(i,TMath::Abs(maximo));
      if (TMath::Abs(maximo)==0)hMaxDev->SetBinContent(i,0.00000001);
     
      hMaxDev->SetBinError(i,err);
    }

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
    hMaxDev->SaveAs(Form("MaxRelDev%s.png",lWhichSystVar.Data()));
    TCanvas* c = new TCanvas("c"," ",1000,800);
    c->SetGridy();
    c->SetGridx();
    hMaxDev->SetStats(kFALSE);  
    hMaxDev->Draw();
    c->SaveAs(Form("MaxRelDev%s.png",lWhichSystVar.Data()));
    
    TCanvas * c1 = 
    new TCanvas("c1"," ",1000,800);
    c1->SetGridy();
    hCut1->GetXaxis()->SetRangeUser(0.7,8.);
    hCut1->GetYaxis()->SetRangeUser(0.9,1.1);
    hCut1->SetYTitle("Yield^{syst-cut} / Yield^{def-cut}");
    hCut1->SetTitle(Form("%s",lWhichSystVar.Data()));
    hCut1->GetYaxis()->SetTitleSize(0.05);
    hCut1->GetYaxis()->SetTitleOffset(0.9);
    hCut1->GetXaxis()->SetTitleSize(0.05);
    hCut1->GetXaxis()->SetTitleOffset(0.8);
    hCut1->SetMarkerStyle(20);
    hCut2->SetMarkerStyle(20);
    hCut1->SetMarkerSize(1.1);
    hCut2->SetMarkerSize(1.1);
    if (lWhichSystVar == "NumberOfCrossedRows"){
      hCut1->SetMarkerColor(kGreen+3);
      hCut2->SetMarkerColor(kAzure+9);
      hCut1->SetLineColor(kGreen+3);
      hCut2->SetLineColor(kAzure+9);
    }
    else if (lWhichSystVar == "ProperLifetime"){
      hCut1->SetMarkerColor(kBlue+1);
      hCut2->SetMarkerColor(kGreen+3);
      hCut1->SetLineColor(kBlue+1);
      hCut2->SetLineColor(kGreen+3);
    }
    hCut1->Draw();
    hCut2->Draw("SAME");  
    TLegend* legend = new TLegend (0.6,0.7,0.8,0.9);
    if (lWhichSystVar == "NumberOfCrossedRows"){
      legend->AddEntry(hCut1,"tight","LEP");
      legend->AddEntry(hCut2,"very tight","LEP");
    }
    else if (lWhichSystVar == "ProperLifetime"){
      legend->AddEntry(hCut1,"loose","LEP");
      legend->AddEntry(hCut2,"tight","LEP");
    }
    legend->Draw("SAME");
    legend->SetTextSize(0.05);
    hCut1->SetStats(kFALSE);   
    c1->SaveAs(Form("%s.png",lWhichSystVar.Data()));
    delete c;
    delete c1;
  
  }

  //################# TPCdEdx Selection and Sig Extr Bin Counting ##########################
  else if (lWhichSystVar == "TPCPIDNSigmas" || lWhichSystVar == "SigmaForSignalExtraction") {
    
    TH1D* hCut1 = 0x0;
    TH1D* hCut2 = 0x0;
    TH1D* hCut3 = 0x0;

    if (lWhichSystVar == "TPCPIDNSigmas"){
      lHistPt[0]->SetName("MB");
      lHistPt[1]->SetName("very loose");
      lHistPt[2]->SetName("loose");
      lHistPt[3]->SetName("tight");
      
      hCut1 = (TH1D*)lHistPt[0]->Clone("very loose");
      hCut2 = (TH1D*)lHistPt[0]->Clone("loose");
      hCut3 = (TH1D*)lHistPt[0]->Clone("tight");
    }
    
    else if (lWhichSystVar == "SigmaForSignalExtraction"){
      lHistPt[0]->SetName("MB");
      lHistPt[1]->SetName("loose");
      lHistPt[2]->SetName("tight");
      lHistPt[3]->SetName("very tight");
      
      hCut1 = (TH1D*)lHistPt[0]->Clone("loose");
      hCut2 = (TH1D*)lHistPt[0]->Clone("tight");
      hCut3 = (TH1D*)lHistPt[0]->Clone("very tight");
    }
            
    for (int i = 1; i<= lHistPt[0]->GetNbinsX(); i++){
      hCut1->SetBinContent(i,lHistPt[1]->GetBinContent(i)/lHistPt[0]->GetBinContent(i));
      hCut1->SetBinError(i,ErrorInRatio(lHistPt[1]->GetBinContent(i),lHistPt[1]->GetBinError(i),lHistPt[0]->GetBinContent(i),lHistPt[0]->GetBinError(i)));
      hCut2->SetBinContent(i,lHistPt[2]->GetBinContent(i)/lHistPt[0]->GetBinContent(i));
      hCut2->SetBinError(i,ErrorInRatio(lHistPt[2]->GetBinContent(i),lHistPt[2]->GetBinError(i),lHistPt[0]->GetBinContent(i),lHistPt[0]->GetBinError(i)));
      hCut3->SetBinContent(i,lHistPt[3]->GetBinContent(i)/lHistPt[0]->GetBinContent(i));
      hCut3->SetBinError(i,ErrorInRatio(lHistPt[3]->GetBinContent(i),lHistPt[3]->GetBinError(i),lHistPt[0]->GetBinContent(i),lHistPt[0]->GetBinError(i)));
    }

    hMaxDev = (TH1D*)lHistPt[0]->Clone("hMaxDev");
    for (int i = 1; i<= lHistPt[0]->GetNbinsX(); i++){
      double maximo = max(max(TMath::Abs(hCut1->GetBinContent(i)-1) , TMath::Abs(hCut2->GetBinContent(i)-1)),TMath::Abs(hCut3->GetBinContent(i)-1));
      double err = 0.;
      if (maximo == TMath::Abs(hCut1->GetBinContent(i)-1)) err = hCut1->GetBinError(i);
      if (maximo == TMath::Abs(hCut2->GetBinContent(i)-1)) err = hCut2->GetBinError(i);
      if (maximo == TMath::Abs(hCut3->GetBinContent(i)-1)) err = hCut3->GetBinError(i);
      hMaxDev->SetBinContent(i,TMath::Abs(maximo));
      if (TMath::Abs(maximo)==0)hMaxDev->SetBinContent(i,0.00000001);
     
      hMaxDev->SetBinError(i,err);
    }

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
    hMaxDev->SaveAs(Form("MaxRelDev%s.png",lWhichSystVar.Data()));
    TCanvas* c = new TCanvas("c"," ",1000,800);
    c->SetGridy();
    c->SetGridx();
    hMaxDev->SetStats(kFALSE);  
    hMaxDev->Draw();
    c->SaveAs(Form("MaxRelDev%s.png",lWhichSystVar.Data()));
    
    TCanvas * c1 = new TCanvas("c1"," ",1000,800);
    c1->SetGridy();
    hCut1->GetXaxis()->SetRangeUser(0.7,8.);
    hCut1->GetYaxis()->SetRangeUser(0.9,1.1);
    hCut1->SetYTitle("Yield^{syst-cut} / Yield^{def-cut}");
    hCut1->SetTitle(Form("%s",lWhichSystVar.Data()));
    hCut1->GetYaxis()->SetTitleSize(0.05);
    hCut1->GetYaxis()->SetTitleOffset(0.9);
    hCut1->GetXaxis()->SetTitleSize(0.05);
    hCut1->GetXaxis()->SetTitleOffset(0.8);
    hCut1->SetMarkerStyle(20);
    hCut2->SetMarkerStyle(20);
    hCut3->SetMarkerStyle(20);
    hCut1->SetMarkerSize(1.1);
    hCut2->SetMarkerSize(1.1);
    hCut3->SetMarkerSize(1.1);
    if (lWhichSystVar == "TPCPIDNSigmas"){
      hCut1->SetMarkerColor(kRed);
      hCut2->SetMarkerColor(kBlue+1);
      hCut3->SetMarkerColor(kGreen+3);
      hCut1->SetLineColor(kRed);
      hCut2->SetLineColor(kBlue+1);
      hCut3->SetLineColor(kGreen+3);
    }
    else if (lWhichSystVar == "SigmaForSignalExtraction"){
      hCut3->SetMarkerColor(kAzure+9);
      hCut1->SetMarkerColor(kBlue+1);
      hCut2->SetMarkerColor(kGreen+3);
      hCut3->SetLineColor(kAzure+9);
      hCut1->SetLineColor(kBlue+1);
      hCut2->SetLineColor(kGreen+3);
    }
    hCut1->Draw();
    hCut2->Draw("SAME");  
    hCut3->Draw("SAME");
    TLegend* legend = new TLegend (0.6,0.7,0.9,0.9);
    if (lWhichSystVar == "TPCPIDNSigmas"){
      legend->AddEntry(hCut1,"very loose","LEP");
      legend->AddEntry(hCut2,"loose","LEP");
      legend->AddEntry(hCut3,"tight","LEP");
    }
    else if (lWhichSystVar == "SigmaForSignalExtraction"){
      legend->AddEntry(hCut1,"loose","LEP");
      legend->AddEntry(hCut2,"tight","LEP");
      legend->AddEntry(hCut3,"very tight","LEP");
    }
    legend->SetTextSize(0.05);
    legend->Draw("SAME");
    hCut1->SetStats(kFALSE); 
    c1->SaveAs(Form("%s.png",lWhichSystVar.Data()));
    delete c;
    delete c1;
  
  }
  
  
  else{
    //set names
    cout << "ME HERE" << endl;
    lHistPt[0]->SetName("MB");
    lHistPt[1]->SetName("very loose");
    lHistPt[2]->SetName("loose");
    lHistPt[3]->SetName("tight");
    lHistPt[4]->SetName("very tight");
    
    TH1D* hCut1 = (TH1D*)lHistPt[0]->Clone("very loose");
    TH1D* hCut2 = (TH1D*)lHistPt[0]->Clone("loose");
    TH1D* hCut3 = (TH1D*)lHistPt[0]->Clone("tight");
    TH1D* hCut4 = (TH1D*)lHistPt[0]->Clone("very tight");
    
    for (int i = 1; i<= lHistPt[0]->GetNbinsX(); i++){
      hCut1->SetBinContent(i,lHistPt[1]->GetBinContent(i)/lHistPt[0]->GetBinContent(i));
      hCut1->SetBinError(i,ErrorInRatio(lHistPt[1]->GetBinContent(i),lHistPt[1]->GetBinError(i),lHistPt[0]->GetBinContent(i),lHistPt[0]->GetBinError(i)));
      hCut2->SetBinContent(i,lHistPt[2]->GetBinContent(i)/lHistPt[0]->GetBinContent(i));
      hCut2->SetBinError(i,ErrorInRatio(lHistPt[2]->GetBinContent(i),lHistPt[2]->GetBinError(i),lHistPt[0]->GetBinContent(i),lHistPt[0]->GetBinError(i)));
      hCut3->SetBinContent(i,lHistPt[3]->GetBinContent(i)/lHistPt[0]->GetBinContent(i));
      hCut3->SetBinError(i,ErrorInRatio(lHistPt[3]->GetBinContent(i),lHistPt[3]->GetBinError(i),lHistPt[0]->GetBinContent(i),lHistPt[0]->GetBinError(i)));
      hCut4->SetBinContent(i,lHistPt[4]->GetBinContent(i)/lHistPt[0]->GetBinContent(i));
      hCut4->SetBinError(i,ErrorInRatio(lHistPt[4]->GetBinContent(i),lHistPt[4]->GetBinError(i),lHistPt[0]->GetBinContent(i),lHistPt[0]->GetBinError(i)));
    }

    hMaxDev = (TH1D*)lHistPt[0]->Clone("hMaxDev");
    hMaxDev->Reset();
    for (int i = 1; i<= lHistPt[0]->GetNbinsX(); i++){
      double maximo = max(
         max(
            TMath::Abs(hCut1->GetBinContent(i)-1) , TMath::Abs(hCut2->GetBinContent(i)-1)),
          max(
           TMath::Abs(hCut3->GetBinContent(i)-1),TMath::Abs(hCut4->GetBinContent(i)-1))
                      );
      double err = 0.;
      if (maximo == TMath::Abs(hCut1->GetBinContent(i)-1)) err = hCut1->GetBinError(i);
      if (maximo == TMath::Abs(hCut2->GetBinContent(i)-1)) err = hCut2->GetBinError(i);
      if (maximo == TMath::Abs(hCut3->GetBinContent(i)-1)) err = hCut3->GetBinError(i);
      if (maximo == TMath::Abs(hCut4->GetBinContent(i)-1)) err = hCut4->GetBinError(i);
      hMaxDev->SetBinContent(i,maximo);
      if (TMath::Abs(maximo)==0)hMaxDev->SetBinContent(i,0.000000000001);
     
      hMaxDev->SetBinError(i,err);
    }

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
    hMaxDev->SaveAs(Form("MaxRelDev%s.png",lWhichSystVar.Data()));
    TCanvas* c = new TCanvas("c"," ",1000,800);
    c->SetGridy();
    c->SetGridx();
    hMaxDev->SetStats(kFALSE);  
    hMaxDev->Draw();
    c->SaveAs(Form("MaxRelDev%s.png",lWhichSystVar.Data()));

    TCanvas * c1 =  new TCanvas("c1"," ",1000,800);
    c1->SetGridy();
    hCut1->GetXaxis()->SetRangeUser(0.7,8.);
    hCut1->GetYaxis()->SetRangeUser(0.9,1.1);
    hCut1->SetYTitle("Yield^{syst-cut} / Yield^{def-cut}");
    hCut1->SetTitle(Form("%s",lWhichSystVar.Data()));
    hCut1->GetYaxis()->SetTitleSize(0.05);
    hCut1->GetYaxis()->SetTitleOffset(0.9);
    hCut1->GetXaxis()->SetTitleSize(0.05);
    hCut1->GetXaxis()->SetTitleOffset(0.8);
    hCut1->SetMarkerStyle(20);
    hCut2->SetMarkerStyle(20);
    hCut3->SetMarkerStyle(20);
    hCut4->SetMarkerStyle(20);
    hCut1->SetMarkerSize(1.1);
    hCut2->SetMarkerSize(1.1);
    hCut3->SetMarkerSize(1.1);
    hCut4->SetMarkerSize(1.1);
    hCut1->SetMarkerColor(kRed);
    hCut2->SetMarkerColor(kBlue+1);
    hCut3->SetMarkerColor(kGreen+3);
    hCut4->SetMarkerColor(kAzure+9);
    hCut1->SetLineColor(kRed);
    hCut2->SetLineColor(kBlue+1);
    hCut3->SetLineColor(kGreen+3);
    hCut4->SetLineColor(kAzure+9);
    hCut1->Draw();
    hCut2->Draw("SAME");	
    hCut3->Draw("SAME");
    hCut4->Draw("SAME");
    TLegend* legend = new TLegend (0.6,0.7,0.9,0.9);
    legend->AddEntry(hCut1,"very loose","LEP");
    legend->AddEntry(hCut2,"loose","LEP");
    legend->AddEntry(hCut3,"tight","LEP");
    legend->AddEntry(hCut4,"very tight","LEP");
    legend->SetTextSize(0.05);
    legend->Draw("SAME");
    hCut1->SetStats(kFALSE); 
    c1->SaveAs(Form("%s.png",lWhichSystVar.Data()));
    delete c;
    delete c1;
   
  }
 
  return hMaxDev;  
  cout<<endl<<"---> Plotted Systematics for: "<< lWhichSystVar << endl;
           


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
