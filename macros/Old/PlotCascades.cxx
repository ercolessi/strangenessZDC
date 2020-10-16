#include <iostream>
#include <fstream>

double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr );

void makeSystPlotsCascades(TString lCascType = "XiMinus", TString lWhichEstimator = "V0M", Double_t lMultBoundLo = 0.0, Double_t lMultBoundHi = 100.0, TString lWhichSystVar = "V0Radius");

void PlotCascades(TString Type = "XiMinus", Double_t Low = 0.0, Double_t High = 100.0){

  makeSystPlotsCascades(Type,"V0M",Low,High,"V0Radius");
  //makeSystPlotsCascades(Type,"V0M",Low,High,"CascRadius");
  //makeSystPlotsCascades(Type,"V0M",Low,High,"DCABachToPV");
  //makeSystPlotsCascades(Type,"V0M",Low,High,"DCAV0ToPV");
  //makeSystPlotsCascades(Type,"V0M",Low,High,"DCANegToPV");
  //makeSystPlotsCascades(Type,"V0M",Low,High,"DCAPosToPV");
  //makeSystPlotsCascades(Type,"V0M",Low,High,"DCAV0Daughters");
  //makeSystPlotsCascades(Type,"V0M",Low,High,"DCACascDaughters");
  //makeSystPlotsCascades(Type,"V0M",Low,High,"V0CosPA");
  //makeSystPlotsCascades(Type,"V0M",Low,High,"CascCosPA");
  //makeSystPlotsCascades(Type,"V0M",Low,High,"V0Mass");
  //makeSystPlotsCascades(Type,"V0M",Low,High,"TPCPIDNSigmas");
  //makeSystPlotsCascades(Type,"V0M",Low,High,"TPCNClusters");
  //makeSystPlotsCascades(Type,"V0M",Low,High,"ProperLifetime");
  //makeSystPlotsCascades(Type,"V0M",Low,High,"CompetingSpecies");
  //makeSystPlotsCascades(Type,"V0M",Low,High,"SigExtBinCount");
  //makeSystPlotsCascades(Type,"V0M",Low,High,"SigExtFitting"); 
}


//---------------------------------------------------------------
void makeSystPlotsCascades(TString lCascType = "XiMinus",
  TString lWhichEstimator = "V0M",
  Double_t lMultBoundLo = 0.0,
  Double_t lMultBoundHi = 100.0,
  TString lWhichSystVar = "V0Radius"){

  cout<<"  "<<endl;
  cout<<"-------------------------------------------------------------"<<endl;
  cout<<"        Macro to make plots of systematics for Cascades       "<<endl;
  cout<<"-------------------------------------------------------------"<<endl;
  cout<<endl;
  
  //Set data files
  TString lSystFile = "Results-Systematics";
    if( lWhichEstimator.Contains("V0") ) lSystFile.Append( Form("-%s-13TeV-%s-%%s_%03.0f_%03.0f_ZDC_000_100-",lCascType.Data(),lWhichEstimator.Data(),(Int_t)(100.*lMultBoundLo),(Int_t)(100.*lMultBoundHi)) );
    else                                 lSystFile.Append( Form("-%s-13TeV-%s-%s_%03.0f_%03.0f_ZDC_000_100-",lCascType.Data(),lWhichEstimator.Data(),(Int_t)(lMultBoundLo),(Int_t)(lMultBoundHi)) );

  TString lDataFilename[5];
  lDataFilename[0] = "Results-XiMinus-13TeV-V0M-0.0to100.0.root";
  Form("Results-%s-13TeV-%s_%0.1f_%0.1f_ZDC_000_100.root", lCascType.Data(), lWhichEstimator.Data(), (lMultBoundLo), (lMultBoundHi));
  for (int i = 1; i <= 4; i++){
  	lDataFilename[i] = lSystFile + lWhichSystVar + Form("-%i.root",i);
  }
  
  cout << " --- Set Minimum Bias file   :  " << lDataFilename[0].Data() << endl;
  for (int i = 1; i <= 4; i++){
    cout << " --- Set Systematics file #"<< i << " :  " << lDataFilename[i].Data() << endl;
  }

  TFile* InputFile[5];
  for (int i = 0; i < 5; i++){
   InputFile[i] = new TFile(lDataFilename[i].Data());
  }

  //Get histograms from files
  TH1D* lHistPt[5];
  for (int i = 0; i <= 4; i++){
    lHistPt[i] = (TH1D *) InputFile[i]->Get(Form("fHistPt%s", lCascType.Data()));
  }

  //Get histograms from files
  TH1D* lHistPtRaw[5];
  for (int i = 0; i <= 4; i++){
    lHistPtRaw[i] = (TH1D *) InputFile[i]->Get("lInvMassReal/lInvMassRealRawData/fHistPtRaw");
  }

  //set names
  lHistPt[0]->SetName("MB");
  lHistPt[1]->SetName("very loose");
  lHistPt[2]->SetName("loose");
  lHistPt[3]->SetName("tight");
  lHistPt[4]->SetName("very tight");

  TH1D* hCut1 = (TH1D*)lHistPt[0]->Clone("very loose");
  TH1D* hCut2 = (TH1D*)lHistPt[0]->Clone("loose");
  TH1D* hCut3 = (TH1D*)lHistPt[0]->Clone("tight");
  TH1D* hCut4 = (TH1D*)lHistPt[0]->Clone("very tight");

  for (int i = 1; i< lHistPt[0]->GetNbinsX(); i++){
  	hCut1->SetBinContent(i,lHistPt[1]->GetBinContent(i)/lHistPt[0]->GetBinContent(i));
  	hCut1->SetBinError(i,ErrorInRatio(lHistPtRaw[1]->GetBinContent(i),lHistPtRaw[1]->GetBinError(i),lHistPtRaw[0]->GetBinContent(i),lHistPtRaw[0]->GetBinError(i)));
	hCut2->SetBinContent(i,lHistPt[2]->GetBinContent(i)/lHistPt[0]->GetBinContent(i));
  	hCut2->SetBinError(i,ErrorInRatio(lHistPtRaw[2]->GetBinContent(i),lHistPtRaw[2]->GetBinError(i),lHistPtRaw[0]->GetBinContent(i),lHistPtRaw[0]->GetBinError(i)));
  	hCut3->SetBinContent(i,lHistPt[3]->GetBinContent(i)/lHistPt[0]->GetBinContent(i));
  	hCut3->SetBinError(i,ErrorInRatio(lHistPtRaw[3]->GetBinContent(i),lHistPtRaw[3]->GetBinError(i),lHistPtRaw[0]->GetBinContent(i),lHistPtRaw[0]->GetBinError(i)));
  	hCut4->SetBinContent(i,lHistPt[4]->GetBinContent(i)/lHistPt[0]->GetBinContent(i));
  	hCut4->SetBinError(i,ErrorInRatio(lHistPtRaw[4]->GetBinContent(i),lHistPtRaw[4]->GetBinError(i),lHistPtRaw[0]->GetBinContent(i),lHistPtRaw[0]->GetBinError(i)));
  }

  TCanvas * c1 = new TCanvas("c1"," ",1000,800);
  c1->SetGridy();
  hCut1->GetXaxis()->SetRangeUser(0.8,6.5);
  hCut1->GetYaxis()->SetRangeUser(0.8,1.2);
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
  legend->Draw("SAME");
  hCut1->SetStats(kFALSE); 

  cout<<endl<<"---> Plotted Systematics for: "<< lWhichSystVar << endl;
           


}

//---------------------------------------------------------------------------------------------------
double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr ) {
    //Error in a Ratio
    if(B!=0) {
        Double_t errorfromtop = Aerr*Aerr / (B*B) ;
        Double_t errorfrombottom = ((A*A)/(B*B*B*B)) * Berr * Berr;
        return TMath::Sqrt( errorfromtop + errorfrombottom );
    }
    return 1.;
}
