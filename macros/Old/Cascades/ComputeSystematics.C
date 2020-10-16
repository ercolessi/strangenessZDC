#include <iostream>
#include <fstream>
//#define max(a,b) (a>b ? a : b)

double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr );
TH1D* makeSystPlotsV0s(TString lV0Type = "XiMinus", TString lWhichEstimator = "V0M", Double_t lMultBoundLo = 0.0, Double_t lMultBoundHi = 100.0, TString lWhichSystVar = "V0Radius");

void ComputeSystematics(TString Type = "XiMinus", Double_t Low = 10.0, Double_t High = 20.0){
  
  TH1D* hV0Radius = makeSystPlotsV0s(Type,"V0M",Low,High,"V0Radius");
  TH1D* hCascRadius = makeSystPlotsV0s(Type,"V0M",Low,High,"CascRadius");
  TH1D* hNegToPV = makeSystPlotsV0s(Type,"V0M",Low,High,"DCANegToPV");
  TH1D* hPosToPV = makeSystPlotsV0s(Type,"V0M",Low,High,"DCAPosToPV");
  TH1D* hV0ToPV = makeSystPlotsV0s(Type,"V0M",Low,High,"DCAV0ToPV");
  TH1D* hBachToPV = makeSystPlotsV0s(Type,"V0M",Low,High,"DCABachToPV");
  TH1D* hV0Daught = makeSystPlotsV0s(Type,"V0M",Low,High,"DCAV0Daughters");
  TH1D* hCascDaught = makeSystPlotsV0s(Type,"V0M",Low,High,"DCACascDaughters");
  TH1D* hV0CosPA = makeSystPlotsV0s(Type,"V0M",Low,High,"V0CosPA");
  TH1D* hCascCosPA = makeSystPlotsV0s(Type,"V0M",Low,High,"CascCosPA");
  TH1D* hV0Mass = makeSystPlotsV0s(Type,"V0M",Low,High,"V0Mass");
  TH1D* hCompetingSpecies = 0x0;
  if (Type == "OmegaMinus" || Type == "OmegaPlus")
    hCompetingSpecies = makeSystPlotsV0s(Type,"V0M",Low,High,"CompetingSpecies");
  TH1D* hPLT = makeSystPlotsV0s(Type,"V0M",Low,High,"ProperLifetime");
  TH1D* hTPCNSigmas = makeSystPlotsV0s(Type,"V0M",Low,High,"TPCPIDNSigmas");
  TH1D* hTPCNClusters = makeSystPlotsV0s(Type,"V0M",Low,High,"TPCNClusters");
  TH1D* hSigExtBinCount = makeSystPlotsV0s(Type,"V0M",Low,High,"SigExtBinCount");

  //Other contributions:
  Double_t MaterialBudget = 0.04;
  Double_t MultIndipEffic = 0.02;
  Double_t MultIndipFeedDown = 0.02;

  //Topological Systematics
  TH1D* hSystTopological = (TH1D*)hV0Radius->Clone("hSystTopological");
  hSystTopological->Reset();
  Double_t V0Radius   = 0;
  Double_t CascRadius = 0;
  Double_t NegToPV    = 0;
  Double_t PosToPV    = 0;
  Double_t BachToPV   = 0;
  Double_t V0ToPV     = 0;
  Double_t CascDaught = 0;
  Double_t CascCosPA  = 0;
  Double_t V0Daught   = 0;
  Double_t V0CosPA    = 0;

  for (int i = 1; i <= hV0Radius->GetNbinsX(); i++){
    if (TMath::Abs(hV0Radius->GetBinContent(i)/hV0Radius->GetBinError(i)) > 2) V0Radius = hV0Radius->GetBinContent(i); else V0Radius = 0;
    if (TMath::Abs(hCascRadius->GetBinContent(i)/hCascRadius->GetBinError(i)) > 2) CascRadius = hCascRadius->GetBinContent(i); else CascRadius = 0;   
    if (TMath::Abs(hNegToPV->GetBinContent(i)/hNegToPV->GetBinError(i)) > 2) NegToPV= hNegToPV->GetBinContent(i); else  NegToPV = 0;
    if (TMath::Abs(hPosToPV->GetBinContent(i)/hPosToPV->GetBinError(i)) > 2) PosToPV = hPosToPV->GetBinContent(i); else PosToPV = 0;
    if (TMath::Abs(hBachToPV->GetBinContent(i)/hBachToPV->GetBinError(i)) > 2) BachToPV = hBachToPV->GetBinContent(i); else BachToPV = 0;
    if (TMath::Abs(hV0Daught->GetBinContent(i)/hV0Daught->GetBinError(i)) > 2) V0Daught = hV0Daught->GetBinContent(i); else V0Daught = 0;
    if (TMath::Abs(hV0CosPA->GetBinContent(i)/hV0CosPA->GetBinError(i)) > 2) V0CosPA = hV0CosPA->GetBinContent(i); else V0CosPA = 0;
    if (TMath::Abs(hV0ToPV->GetBinContent(i)/hV0ToPV->GetBinError(i)) > 2) V0ToPV = hV0ToPV->GetBinContent(i); else V0ToPV = 0;
    if (TMath::Abs(hCascDaught->GetBinContent(i)/hCascDaught->GetBinError(i)) > 2) CascDaught = hCascDaught->GetBinContent(i); else CascDaught = 0;
    if (TMath::Abs(hCascCosPA->GetBinContent(i)/hCascCosPA->GetBinError(i)) > 2) CascCosPA = hCascCosPA->GetBinContent(i); else CascCosPA = 0;
    
    hSystTopological->SetBinContent(i, 
      TMath::Sqrt( V0Radius*V0Radius + NegToPV*NegToPV + PosToPV*PosToPV + BachToPV*BachToPV +
                   V0Daught*V0Daught + CascDaught*CascDaught + V0CosPA*V0CosPA + CascCosPA*CascCosPA +
                   V0CosPA*V0CosPA + V0ToPV*V0ToPV )    
    );
  }

  //Other Selection Cuts
  TH1D* hSystOthers = (TH1D*)hV0Radius->Clone("hSystOthers");
  hSystOthers->Reset();
  Double_t V0Mass = 0;
  Double_t PLT = 0;
  Double_t CompetingSpecies = 0;
  Double_t TPCNClusters = 0;
  Double_t TPCNSigmas = 0;
  Double_t SigExtBinCount = 0;

  for (int i = 1; i<= hV0Radius->GetNbinsX(); i++){
      if (TMath::Abs(hV0Mass->GetBinContent(i)/hV0Mass->GetBinError(i)) > 2) V0Mass = hV0Mass->GetBinContent(i); else V0Mass = 0;
      if (TMath::Abs(hPLT->GetBinContent(i)/hPLT->GetBinError(i)) > 2) PLT = hPLT->GetBinContent(i); else  PLT= 0;
      if (Type == "OmegaMinus" || Type == "OmegaPlus") {
        if (TMath::Abs(hCompetingSpecies->GetBinContent(i)/hCompetingSpecies->GetBinError(i)) > 2) CompetingSpecies = hCompetingSpecies->GetBinContent(i); else CompetingSpecies = 0;
        } else CompetingSpecies = 0;
      if (TMath::Abs(hTPCNClusters->GetBinContent(i)/hTPCNClusters->GetBinError(i)) > 2) TPCNClusters = hTPCNClusters->GetBinContent(i); else TPCNClusters = 0;
      if (TMath::Abs(hTPCNSigmas->GetBinContent(i)/hTPCNSigmas->GetBinError(i)) > 2) TPCNSigmas = hTPCNSigmas->GetBinContent(i); else TPCNSigmas = 0;
      if (TMath::Abs(hSigExtBinCount->GetBinContent(i)/hSigExtBinCount->GetBinError(i)) > 2) SigExtBinCount = hSigExtBinCount->GetBinContent(i); else SigExtBinCount = 0;
      
      hSystOthers->SetBinContent(i, 
      TMath::Sqrt(V0Mass*V0Mass + PLT*PLT + CompetingSpecies*CompetingSpecies + TPCNClusters*TPCNClusters +
       TPCNSigmas*TPCNSigmas + SigExtBinCount*SigExtBinCount)
    );
  }

  TH1D* hSystTot = (TH1D*)hV0Radius->Clone("hSystTot");  
  hSystTot->Reset();
  for (int i = 1; i<= hV0Radius->GetNbinsX(); i++){
      hSystTot->SetBinContent(i, 
        TMath::Sqrt(hSystOthers->GetBinContent(i)*hSystOthers->GetBinContent(i) + 
          hSystTopological->GetBinContent(i)*hSystTopological->GetBinContent(i)) 
      );
  }

  //Contribution Total
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
  hSystTot->GetYaxis()->SetRangeUser(-0.005,0.1);
  hSystTot->GetYaxis()->SetTitleOffset(1.);
  hSystTot->SetTitle(Form("Systematics contributions for %s",Type.Data()));
  hSystOthers->Draw("SAME");
  hSystTopological->Draw("SAME");
  legend->Draw("SAME");
  cn->SaveAs(Form("images/%s-Systematics-Displayed_%03.0f_%03.0f.png",Type.Data(),Low, High));

  //Contribution Topological 
  TCanvas* cTop = new TCanvas("cTop","",1300,900);
  TLegend* l = new TLegend (0.2,0.55,0.66,0.9);
  l->AddEntry(hV0Radius,"V0 Radius","L");
  l->AddEntry(hCascRadius,"Casc Radius","L");
  l->AddEntry(hNegToPV,"DCA Neg To PV","L");
  l->AddEntry(hPosToPV,"DCA Pos To PV","L");
  l->AddEntry(hBachToPV,"DCA Bach To PV","L");
  l->AddEntry(hV0Daught,"DCA V0 Daughters","L");
  l->AddEntry(hCascDaught,"DCA Casc Daughters","L");
  l->AddEntry(hV0CosPA,"V0 Cosine of PA","L");
  l->AddEntry(hCascCosPA,"Casc Cosine of PA","L");
  l->AddEntry(hV0ToPV,"V0 To PV","L");
  l->SetTextSize(0.03);
  hV0Radius->SetLineColor(kRed+1);
  hCascRadius->SetLineColor(kRed);
  hNegToPV->SetLineColor(kOrange+7);
  hPosToPV->SetLineColor(kOrange);
  hBachToPV->SetLineColor(kYellow);
  hV0Daught->SetLineColor(kSpring+10);
  hCascDaught->SetLineColor(kSpring);
  hV0CosPA->SetLineColor(kAzure+8);
  hCascCosPA->SetLineColor(kAzure);
  hV0ToPV->SetLineColor(kBlue+1);
  hV0Radius->SetMarkerStyle(1);
  hNegToPV->SetMarkerStyle(1);
  hCascRadius->SetMarkerStyle(1);
  hPosToPV->SetMarkerStyle(1);
  hBachToPV->SetMarkerStyle(1);
  hV0Daught->SetMarkerStyle(1);
  hCascDaught->SetMarkerStyle(1);
  hV0CosPA->SetMarkerStyle(1);
  hCascCosPA->SetMarkerStyle(1);
  hV0ToPV->SetMarkerStyle(1);
  hV0Radius->SetLineWidth(2);
  hNegToPV->SetLineWidth(2);
  hCascRadius->SetLineWidth(2);
  hPosToPV->SetLineWidth(2);
  hBachToPV->SetLineWidth(2);
  hV0Daught->SetLineWidth(2);
  hCascDaught->SetLineWidth(2);
  hV0CosPA->SetLineWidth(2);
  hCascCosPA->SetLineWidth(2);
  hV0ToPV->SetLineWidth(2);
  hV0Radius->SetTitle(Form("Contribution of Topological Variables %s",Type.Data()));
  hV0Radius->Draw("L");
  hCascRadius->Draw("SAME L");
  hNegToPV->Draw("SAME L");
  hPosToPV->Draw("SAME L");
  hBachToPV->Draw("SAME L");
  hV0Daught->Draw("SAME L");
  hCascDaught->Draw("SAME L");
  hV0CosPA->Draw("SAME L");
  hCascCosPA->Draw("SAME L");
  hV0ToPV->Draw("SAME L");
  l->Draw("SAME");
  cTop->SaveAs(Form("images/%s-Systematics-Topological_%03.0f_%03.0f.png",Type.Data(),Low,High));

  //Selection contributions
  TCanvas* cSel = new TCanvas("cSel","",1300,900);
  TLegend* l1 = new TLegend (0.14,0.65,0.65,0.9);
  l1->AddEntry(hPLT,"Proper Life Time","L");
  l1->AddEntry(hV0Mass,"V0 Mass","L");
  if (Type == "OmegaMinus" || Type == "OmegaPlus") l1->AddEntry(hCompetingSpecies,"Competing Species","L");
  l1->AddEntry(hTPCNClusters,"TPC N of Clusters","L");
  l1->AddEntry(hTPCNSigmas,"TPC N Sigmas","L");
  l1->AddEntry(hSigExtBinCount,"Sigmas for Sgn Extraction","L");
  l1->SetTextSize(0.03);
  hPLT->SetLineColor(kRed+1);
  hV0Mass->SetLineColor(kOrange+1);
  hSigExtBinCount->SetLineColor(kSpring-1);
  hTPCNClusters->SetLineColor(kAzure-4);
  hTPCNSigmas->SetLineColor(kBlue+2);
  if (Type == "OmegaMinus" || Type == "OmegaPlus") hCompetingSpecies->SetLineColor(kBlack);
  hPLT->SetMarkerStyle(1);
  hPLT->SetLineWidth(2);
  hV0Mass->SetMarkerStyle(1);
  hV0Mass->SetLineWidth(2);
  if (Type == "OmegaMinus" || Type == "OmegaPlus") {
    hCompetingSpecies->SetMarkerStyle(1);
    hCompetingSpecies->SetLineWidth(2);
  }
  hTPCNClusters->SetMarkerStyle(1);
  hTPCNClusters->SetLineWidth(2);
  hTPCNSigmas->SetMarkerStyle(1);
  hTPCNSigmas->SetLineWidth(2);
  hSigExtBinCount->SetMarkerStyle(1);
  hSigExtBinCount->SetLineWidth(2);
  hPLT->SetTitle("Contribution of Selection Variables");
  hPLT->Draw("L");
  hV0Mass->Draw("SAME L");
  if (Type == "OmegaMinus" || Type == "OmegaPlus") hCompetingSpecies->Draw("SAME L");
  hTPCNClusters->Draw("SAME L");
  hTPCNSigmas->Draw("SAME L");
  hSigExtBinCount->Draw("SAME L");
  l1->Draw("SAME");
  cSel->SaveAs(Form("images/%s-Systematics-Selection_%03.0f_%03.0f.png",Type.Data(),Low,High));

  //Corrected Spectra
  TString lDataFilename;
  lDataFilename = Form("Files_%03.0f-%03.0f/Results-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_000_100.root", Low, High, Type.Data(), Low, High);
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
  lHistPt->GetXaxis()->SetRangeUser(1.,6.5);
  hSystPt->GetXaxis()->SetRangeUser(1.,6.5);
  hSystPt->GetYaxis()->SetTitleOffset(4.);
  hSystPt->Draw("E2");
  lHistPt->Draw("SAME");
  s->SaveAs(Form("images/%s-Spettri-Systematics_%03.0f_%03.0f.png",Type.Data(),Low,High));


 
}

//---------------------------------------------------------------
TH1D* makeSystPlotsV0s(
          TString lCascType = "XiMinus",
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
  const int nfiles = 5;  
  TH1D* hMaxDev;
  TH1D* lHistPt[nfiles];
  TH1D* lHistPtRaw[nfiles];
  TH1D* hCut[nfiles-1];
  TString lDataFilename[nfiles];
  TFile* InputFile[nfiles]; 
  

  //Set data files
  TString lSystFile = Form("Files_%03.0f-%03.0f/Results-Systematics",lMultBoundLo,lMultBoundHi);
  lSystFile.Append( Form("-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_000_100-",lCascType.Data(),lMultBoundLo,lMultBoundHi) );
  
  //Files
  lDataFilename[0] = Form("Files_%03.0f-%03.0f/Results-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_000_100.root",lMultBoundLo,lMultBoundHi,lCascType.Data(),lMultBoundLo,lMultBoundHi);
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
  lHistPt[0] = (TH1D *) InputFile[0]->Get(Form("fHistPt%s", lCascType.Data()));
  lHistPtRaw[0] = (TH1D *) InputFile[0]->Get("lInvMassReal/lInvMassRealRawData/fHistPtRaw");
  hMaxDev = (TH1D*)lHistPt[0]->Clone("hMaxDev");
  hMaxDev->Reset();
  for (int i = 1; i < nfiles; i++){
    lHistPt[i] = (TH1D *) InputFile[i]->Get(Form("fHistPt%s", lCascType.Data()));
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
  hMaxDev->GetXaxis()->SetRangeUser(0.7,6.5);
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
  maxdev->SaveAs(Form("images/%s-MaxRelDev%s_%03.0f_%03.0f.png",lCascType.Data(),lWhichSystVar.Data(),lMultBoundLo,lMultBoundHi));

  //Cut Variation
  TCanvas* cutvar = new TCanvas("cutvar"," ",1000,800);
  cutvar->SetGridy();
  hCut[0]->Draw();
  TLegend* legend = new TLegend (0.6,0.75,0.8,0.9);  

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

  for (int k = 0; k < nfiles-1; k++){
    hCut[k]->GetXaxis()->SetRangeUser(0.7,6.5);
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
  hCut[0]->SetStats(kFALSE);
  for (int k = 1; k < nfiles-1; k++){
    hCut[k]->Draw("SAME");
  }
  legend->SetTextSize(0.035);
  legend->Draw("SAME");    
  cutvar->SaveAs(Form("images/%s-CutVar-%s_%03.0f_%03.0f.png",lCascType.Data(),lWhichSystVar.Data(),lMultBoundLo,lMultBoundHi));


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