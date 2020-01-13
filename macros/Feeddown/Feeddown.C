#include <TFile.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TClonesArray.h>
#include <TLatex.h>
#include <TMath.h>
#include <TColor.h>

// STL includes
#include <iostream>
#include <fstream>
using std::cout;
using std::endl;

TF1 * LevyTsallis(const Char_t *name, Double_t mass, Double_t n = 5., Double_t C = 0.001, Double_t norm = 0.01);
double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr );

void Feeddown(const char* s = "Lambda"){
  
  printf("\n  Calculating Feed-Down for %s", s);

  //Files----------------------------------------------------------------------------------------------
  const char *outputname = Form("Feeddown%s.root",s);
  const char *MCname = "~/Strangeness/Efficienze/MonteCarloMerged15f.root";
  TString part(s);
  const char* Mother;
  const char *XiSpectraname;   
  const char *LambdaSpectraname;
  
  if (part == "Lambda") {
    Mother = "Xi-";
    XiSpectraname = "~/Strangeness/Spettri/Xi/XiCorrSpectraMB.root";   
    LambdaSpectraname = "~/Strangeness/Spettri/Lambda/LambdaRawSpectraMB.root";
  }
  else if (part == "AntiLambda"){
    Mother = "Xi+";
    XiSpectraname = "~/Strangeness/Spettri/AntiXi/AntiXiCorrSpectraMB.root";   
    LambdaSpectraname = "~/Strangeness/Spettri/AntiLambda/AntiLambdaRawSpectraMB.root";
  }

  //Variable Definition---------------------------------------------------------------------------------
  //pdg
  int pdg, pdgPos, pdgNeg, pdgMotherCharged, pdgMotherNeutral;
  if (part == "Lambda" ) {
    pdg = 3122;//Lambda
    pdgPos = 2212;//proton
    pdgNeg = -211;//pi-
    pdgMotherCharged = 3312;//Xi-
    pdgMotherNeutral = 3322;//Xi0   
  }
  else if (part == "AntiLambda" ) {
    pdg = -3122;//AntiLambda
    pdgNeg = -2212;//antiproton
    pdgPos = +211;//pi+
    pdgMotherCharged = -3312;//Xi+
    pdgMotherNeutral = -3322;//\barXi0   
  }
  
  //kinetic var
  Double_t temppt,ptL, ptMother, m = 1115.683, Ximass = 1321.31;
  
  //Double_t fptbinlimits[18] = { 0.4,  0.6,  0.8,  1,  1.2,  1.4,  1.6,  1.8,  2,  2.2,  2.5,  2.9,  3.4,  4,  5,  6.5,  8, 10.0};
  //Double_t xibinlimits[24] = {0.,  0.2,  0.4,  0.6,  0.8,  0.9,  1,  1.1,  1.2,  1.3,  1.4,  1.5,  1.7,  1.9,  2.2,  2.6,  3.1,  3.9,  4.9,  6,  7.2,  8.5,  10, 12};
  Double_t fptbinlimits[24] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.6,3.0,4.0,5.0,6.5,8.0,10.0}; 
  Double_t xibinlimits[24] = { //Define Xi Binning
    0.00,  0.10,  0.20,  0.30,  0.40,  0.50,
    0.60,  0.70,  0.80,  0.90,  1.00,  1.20,
    1.40,  1.60,  1.80,  2.00,  2.20,  2.60,
    3.00,  4.00,  5.00,  6.50,  8.00,  10.00
    };
  //constants
  const int kTPCrefit = 64;
  const int ctr = 1;
  Double_t cmin [ctr]= {0.0};//, 1.0, 5.0, 10.0, 15.0 , 30.0, 50.0 , 70.0 ,0.0};
  Double_t cmax [ctr]= {100.0};//, 5.0, 10.0 , 15.0, 30.0 , 50.0 , 70.0 , 100.0 , 100.0};
  int cevt [ctr] = {9};
  int cevtL [ctr] = {11};
  Long_t fptbinnumb = sizeof(fptbinlimits)/sizeof(Double_t) - 1; 
  Long_t xibinnumb = sizeof(xibinlimits)/sizeof(Double_t) - 1;
 
  //Histos
  TH1D* hPTXiCorr[ctr];
  TH1D* hPTLambdaRaw[ctr];
  TH1F* fHistPt = new TH1F("fHistPt","Dummy;p_{T} (GeV/c);Counts",fptbinnumb,fptbinlimits);
  TH1F* fHistPtXiReference =  new TH1F("fHistPtXiReference","#Xi candidates uncorrected p_{T};p_{T} (GeV/c);Counts",xibinnumb,xibinlimits);
 
  //File MC
  TFile *fMC = new TFile(MCname);
  TTree *tMC = (TTree *) fMC->Get("PWGLF_StrVsMult_MC/fTreeV0");
  TTree *tMCc = (TTree *) fMC->Get("PWGLF_StrVsMult_MC/fTreeCascade");
  TList *l = (TList *) fMC->Get("PWGLF_StrVsMult_MC/cList");
  const int lNCandidatesMCFD = tMC->GetEntries();
  //Get from MC
  TH3D* f3dHistGenPtVsYVsMultV0FeedDown;
  if (part == "Lambda") f3dHistGenPtVsYVsMultV0FeedDown = (TH3D *) l->FindObject("fHistGeneratedPtVsYVsCentralityXiMinus");
  else if (part == "AntiLambda") f3dHistGenPtVsYVsMultV0FeedDown = (TH3D *) l->FindObject("fHistGeneratedPtVsYVsCentralityXiPlus");
  
  //File Spettri
  TFile *fR = new TFile(LambdaSpectraname);
  TFile *fS = new TFile(XiSpectraname);
  //Get from Spettri
  TH1D* hevt = (TH1D*) fS->Get("hevt");
  TH1D* hevtL = (TH1D*) fR->Get("hevt");
  for (int j = 0 ; j < ctr ; j++) {
    hPTXiCorr[j] = (TH1D*) fS->Get(Form("hPT_%f-%f",cmin[j],cmax[j]));
    hPTLambdaRaw[j] = (TH1D*) fR->Get(Form("hPT_%f-%f",cmin[j],cmax[j]));    
  }
   
  //Feeddown Matrix ----------------------------------------------------------------------------------------------
  // Define Feeddown matrix
  Double_t lFeedDownMatrix [100][100]; // FeedDownMatrix [Lambda Bin][Xi Bin];
  for(Int_t ilb = 0; ilb<100; ilb++) {
    for(Int_t ixb = 0; ixb<100; ixb++) {
      lFeedDownMatrix[ilb][ixb]=0;
    }
  }  
  //Feeddown: will use a double-coordinate system:
  // ( lambda bin, xi bin ) all the time!
  Int_t lWeAreAtBin = 0;         // Lambda bin
  Int_t lWeAreAtXiBin = 0; // Xi Bin  
  
  cout<<"\n\n--------------- MC Data File Loop, Feeddown  -----------"<<endl;
  Long_t lOneTenthOfNCandidatesMCFD = ((double)(lNCandidatesMCFD) / 10. );
   
  //Loop over V0 events
  for(int i=0;i < lNCandidatesMCFD;i++){
    
    if( i % lOneTenthOfNCandidatesMCFD == 0 )
      cout<<" Currently at candidate........: "<<i<<" / "<<lNCandidatesMCFD<<" ( "<<(long)(((double)(i)/(double)(lNCandidatesMCFD))*(100.+1e-3))<<"% )"<<endl;
    
    tMC->GetEvent(i);
    
    //MC PDG check    
    if(
       tMC->GetLeaf("fTreeVariablePIDMother")->GetValue()  == pdgMotherCharged ||
       tMC->GetLeaf("fTreeVariablePIDMother")->GetValue() == pdgMotherNeutral
       ) {
      if(tMC->GetLeaf("fTreeVariablePID")->GetValue() != pdg) continue;
      if(tMC->GetLeaf("fTreeVariablePIDPositive")->GetValue() != pdgPos) continue;
      if(tMC->GetLeaf("fTreeVariablePIDNegative")->GetValue() != pdgNeg) continue;
      //Is Mother Primary?
      //if(tMC->GetLeaf("fTreeVariablePrimaryStatus")->GetValue() != 2) continue;
      //if(tMC->GetLeaf("fTreeVariableRapMother")->GetValue() >= 10.0) continue;
      //if(tMC->GetLeaf("fTreeVariableRapMother")->GetValue() <= -10.0) continue;
      if(tMC->GetLeaf("fTreeVariablePrimaryStatusMother")->GetValue() != 1) continue;
      
      //MOTHER Xi
      ptMother = tMC->GetLeaf("fTreeVariablePtMother")->GetValue();
      
      //Topological 
      if (tMC->GetLeaf("fTreeVariableV0Radius")->GetValue() <= 0.5) continue;
      if (tMC->GetLeaf("fTreeVariableDcaNegToPrimVertex")->GetValue() <= 0.06) continue;
      if (tMC->GetLeaf("fTreeVariableDcaPosToPrimVertex")->GetValue() <= 0.06 ) continue;
      if (tMC->GetLeaf("fTreeVariableV0CosineOfPointingAngle")->GetValue() <= 0.97  ) continue;
      if (tMC->GetLeaf("fTreeVariableDcaV0Daughters")->GetValue() >= 1.0  ) continue;
      //MC Rapidity
      if (tMC->GetLeaf("fTreeVariableRapMC")->GetValue() >= 0.5) continue;
      if (tMC->GetLeaf("fTreeVariableRapMC")->GetValue() <= -0.5) continue;
      //Proper Lifetime
      Double_t L = tMC->GetLeaf("fTreeVariableV0Radius")->GetValue();
      Double_t pt = tMC->GetLeaf("fTreeVariablePtMC")->GetValue();
      Double_t y = tMC->GetLeaf("fTreeVariableRapLambda")->GetValue();
      Double_t p_m = TMath::Sqrt(pt*pt*(TMath::CosH(y))*(TMath::CosH(y))+ (TMath::SinH(y))*(TMath::SinH(y)));      
      if ((L/p_m) >= 30) continue; 
      //Other cuts
      if (tMC->GetLeaf("fTreeVariableNegEta")->GetValue() >= 0.8) continue;
      if (tMC->GetLeaf("fTreeVariableNegEta")->GetValue() <= -0.8) continue;
      if (tMC->GetLeaf("fTreeVariablePosEta")->GetValue() >= 0.8) continue;
      if (tMC->GetLeaf("fTreeVariablePosEta")->GetValue() <= -0.8) continue;
      if (tMC->GetLeaf("fTreeVariableLeastNbrCrossedRows")->GetValue() < 70) continue;
      if (tMC->GetLeaf("fTreeVariableLeastRatioCrossedRowsOverFindable")->GetValue() < 0.8) continue;
      //kTPCrefit
      if ((int(tMC->GetLeaf("fTreeVariablePosTrackStatus")->GetValue())&kTPCrefit) == 0){
	if ((int(tMC->GetLeaf("fTreeVariableNegTrackStatus")->GetValue())&kTPCrefit) == 0){

	  ptL = tMC->GetLeaf("fTreeVariablePtMC")->GetValue();

	  lWeAreAtBin   = fHistPt->FindBin(ptL)-1;
	  if(lWeAreAtBin == -1) lWeAreAtBin   = 99; //UnderFlow, special treatment
	  lWeAreAtXiBin = fHistPtXiReference->FindBin(ptMother)-1;
	  if(lWeAreAtXiBin == -1) lWeAreAtXiBin = 99; //UnderFlow, special treatment
	  //Populate Feeddown Matrix
	  //cout<<" Populate at coordinates "<<lWeAreAtBin<<" vs "<<lWeAreAtXiBin<<" ....."<<endl;
	  lFeedDownMatrix[lWeAreAtBin][lWeAreAtXiBin]++;
	}
      }
    }
  }
  cout<<"--------------- Loop Completed -------------------------"<<endl;
  cout<<endl;

  //Xi Gen Histo Projected
  TH1D *fHistDummyV0FeedDown = f3dHistGenPtVsYVsMultV0FeedDown->ProjectionX("fHistDummyV0FeedDown",
			       f3dHistGenPtVsYVsMultV0FeedDown->GetYaxis()->FindBin(-0.5+1e-2),  
                               f3dHistGenPtVsYVsMultV0FeedDown->GetYaxis()->FindBin(+0.5-1e-2), 
			       f3dHistGenPtVsYVsMultV0FeedDown->GetZaxis()->FindBin(-1),
			       f3dHistGenPtVsYVsMultV0FeedDown->GetZaxis()->FindBin(101));
									   
  TH1F *fHistMCCountbyptXiFeedDown = new TH1F("fHistMCCountbyptXiFeedDown","#Xi MC count;p_{T} (GeV/c);Counts",xibinnumb,xibinlimits);
  //Double_t temppt, y; --- declared before
  for(long i = 1; i<fHistDummyV0FeedDown->GetNbinsX()+1; i++) {
    temppt = fHistDummyV0FeedDown->GetXaxis()->GetBinCenter(i);
    for(long filling = 0; filling<fHistDummyV0FeedDown->GetBinContent(i); filling++) {
      fHistMCCountbyptXiFeedDown->Fill(temppt);
    }
  }
  
  cout<<"------------ Computing actual feeddown matrix------------"<<endl;
  
  TF1* lLevyFitXiFeedDown = LevyTsallis("LevyTsallis",Ximass);
  hPTXiCorr[0]->Fit(lLevyFitXiFeedDown,"","",0.7,8.);
  
  /* new TCanvas;
  hPTXiCorr[0]->GetXaxis()->SetRangeUser(0.6,8.);
  hPTXiCorr[0]->Draw();
  lLevyFitXiFeedDown->Draw("same");*/

  printf("\n------------- LevyTsallis Fit on %s (real-corrected) ------\n",Mother);
  
  Double_t lFeedDownEfficiency[100][100];
  Double_t lFeedDownEfficiencyError[100][100];
  Double_t lFeedDownEfficiencyXiCorrected[100][100];

  ofstream XiPlus("XiPlus.txt");
  
  for(Int_t ilb = 0; ilb<fptbinnumb; ilb++) {
    for(Int_t ixb = 0; ixb<xibinnumb; ixb++) {
      if( fHistMCCountbyptXiFeedDown->GetBinContent(ixb+1)!= 0 ) { //avoid div by zero
	lFeedDownEfficiency[ilb][ixb]=((double)(lFeedDownMatrix[ilb][ixb])) / ((double)( fHistMCCountbyptXiFeedDown->GetBinContent(ixb+1) )) ;
	lFeedDownEfficiencyError[ilb][ixb]=ErrorInRatio(
                                                           ((double)(lFeedDownMatrix[ilb][ixb])),
                                                           TMath::Sqrt((double)(lFeedDownMatrix[ilb][ixb])),
                                                           ((double)( fHistMCCountbyptXiFeedDown->GetBinContent(ixb+1) )),
                                                           TMath::Sqrt((double)( fHistMCCountbyptXiFeedDown->GetBinContent(ixb+1) ))
							);
	lFeedDownEfficiencyXiCorrected[ilb][ixb] = ((double)(lFeedDownMatrix[ilb][ixb]))
	                                           /((double)( fHistMCCountbyptXiFeedDown->GetBinContent(ixb+1)))
	                                           *hPTXiCorr[0]->GetBinContent(ixb+1)*hPTXiCorr[0]->GetBinWidth(ixb+1)
	                                           //*(lLevyFitXiFeedDown->Integral(xibinlimits[ixb],xibinlimits[ixb+1]))
	                                           *hevt->GetBinContent(cevt[0]);
	
	if (((double)(lFeedDownMatrix[ilb][ixb]))!=0) XiPlus << " XiPlus " << " Matrix Value " << ((double)(lFeedDownMatrix[ilb][ixb])) << "\n Xi Gen " << ((double)( fHistMCCountbyptXiFeedDown->GetBinContent(ixb+1))) << " \n PTXiCorr " <<  hPTXiCorr[0]->GetBinContent(ixb+1)*hPTXiCorr[0]->GetBinWidth(ixb+1) << endl;
      } else {
	lFeedDownEfficiency[ilb][ixb] = 0;
	lFeedDownEfficiencyError[ilb][ixb] = 0;
	lFeedDownEfficiencyXiCorrected[ilb][ixb] = 0;
      }
    }
  }

  //Create FD TH2Fs for storing
  TH2F* f2dFeedDownMatrix = new TH2F("f2dFeedDownMatrix","",fptbinnumb,fptbinlimits,xibinnumb,xibinlimits);
  TH2F* f2dFeedDownEfficiency = new TH2F("f2dFeedDownEfficiency","",fptbinnumb,fptbinlimits,xibinnumb,xibinlimits);
  TH2F* f2dFeedDownEfficiencyXiCorrected = new TH2F("f2dFeedDownEfficiencyXiCorrected","",fptbinnumb,fptbinlimits,xibinnumb,xibinlimits);
  for(Int_t ilb = 0; ilb<fptbinnumb; ilb++) {
    for(Int_t ixb = 0; ixb<xibinnumb; ixb++) {
      f2dFeedDownMatrix->SetBinContent(ilb+1,ixb+1,lFeedDownMatrix[ilb][ixb]);
      f2dFeedDownEfficiency->SetBinContent(ilb+1,ixb+1,lFeedDownEfficiency[ilb][ixb]);
      f2dFeedDownEfficiency->SetBinError(ilb+1,ixb+1,lFeedDownEfficiencyError[ilb][ixb]);
      f2dFeedDownEfficiencyXiCorrected->SetBinContent(ilb+1,ixb+1,lFeedDownEfficiencyXiCorrected[ilb][ixb]);
    }
  }
  
  new TCanvas;
  f2dFeedDownMatrix->Draw("colz");

  new TCanvas;
  f2dFeedDownEfficiency->Draw("colz");
  
  new TCanvas;
  f2dFeedDownEfficiencyXiCorrected->Draw("colz");

  cout<<"\n\n" <<"--------------- Generated Xi- Dump (real-corrected) ----"<<endl;

  Double_t lProducedXi[100];
  
  for(Int_t ixb = 0; ixb<xibinnumb; ixb++) {
    lProducedXi[ixb] = hevt->GetBinContent(cevt[0])*hPTXiCorr[0]->GetBinContent(ixb+1)*hPTXiCorr[0]->GetBinWidth(ixb+1);
      //*(lLevyFitXiFeedDown->Integral(xibinlimits[ixb],xibinlimits[ixb+1]))	                                          
  }
  
  //=== Computing actual feeddown numbers... ===
  Double_t lFeedDownToSubtract[100];
  Double_t lFeedDownToSubtractError[100];
  for(Int_t ilb = 0; ilb < fptbinnumb; ilb++) {
    lFeedDownToSubtract[ilb] = 0;
    lFeedDownToSubtractError[ilb] = 0;
    for(Int_t ixb = 0; ixb < xibinnumb; ixb++) {
      lFeedDownToSubtract[ilb] += lProducedXi[ixb]*lFeedDownEfficiency[ilb][ixb];
      lFeedDownToSubtractError[ilb] += TMath::Power( (lProducedXi[ixb]*lFeedDownEfficiencyError[ilb][ixb]), 2 );
    }
    lFeedDownToSubtractError[ilb] = TMath::Sqrt(lFeedDownToSubtractError[ilb]);
  }  
  TH1F* fHistFeeddownSubtraction = new TH1F("fHistFeeddownSubtraction","FD Subtraction;p_{T} (GeV/c);Ratio removed",fptbinnumb,fptbinlimits);
  fHistFeeddownSubtraction->SetDirectory(0);
  Double_t LambdaRaw[100], LambdaRawError[100];
  for(Int_t ilb = 0; ilb<fptbinnumb; ilb++) {
    LambdaRaw[ilb] = hPTLambdaRaw[0]->GetBinContent(ilb+1)*hPTLambdaRaw[0]->GetBinWidth(ilb+1)*hevtL->GetBinContent(cevtL[0]);
    LambdaRawError[ilb] = hPTLambdaRaw[0]->GetBinError(ilb+1)*hPTLambdaRaw[0]->GetBinWidth(ilb+1)*hevtL->GetBinContent(cevtL[0]);
    
    fHistFeeddownSubtraction->SetBinContent( ilb+1, ((double)(lFeedDownToSubtract[ilb]))
					     / ((double)LambdaRaw[ilb]) );
    fHistFeeddownSubtraction->SetBinError( ilb+1 , ErrorInRatio(((double)(lFeedDownToSubtract[ilb])),
                                                     ((double)(lFeedDownToSubtractError[ilb])),
                                                     ((double)(LambdaRaw[ilb]) ),
                                                     ((double)(LambdaRawError[ilb]) ))
                                                     );    
  }
      
  new TCanvas;
  fHistFeeddownSubtraction->Draw();
    
  //Write
  TFile* Write = new TFile (outputname, "recreate");
  f2dFeedDownMatrix->Write();
  f2dFeedDownEfficiency->Write();
  f2dFeedDownEfficiencyXiCorrected->Write();
  fHistFeeddownSubtraction->Write();
  Write->Close();
  




  /*
  TH1D* f2dFeedDownEfficiencyXiCorrectedProj;
  TH1D* hFinalFeeddown;
  f2dFeedDownEfficiencyXiCorrectedProj = f2dFeedDownEfficiencyXiCorrected->ProjectionX();
  new TCanvas;
  f2dFeedDownEfficiencyXiCorrectedProj->Draw();
  f2dFeedDownEfficiencyXiCorrectedProj->SetTitle(Form("f2dFeedDownEfficiencyXiCorrectedProj%i",0));
  f2dFeedDownEfficiencyXiCorrectedProj->Sumw2(); 
  
  hFinalFeeddown = new TH1D(Form("hFinalFeeddown%i",0)," ",fptbinnumb,fptbinlimits);
    
  for (Int_t ilb = 0; ilb<fptbinnumb; ilb++) {   
  double entry = f2dFeedDownEfficiencyXiCorrectedProj->GetBinContent(ilb+1)
        /(hPTLambdaRaw[0]->GetBinContent(ilb+1)
	*hPTLambdaRaw[0]->GetBinWidth(ilb+1)
	*hevtL->GetBinContent(cevtL[0])
	);
	
   hFinalFeeddown->SetBinContent(ilb+1,entry);    
  }
  new TCanvas;
  hFinalFeeddown->GetXaxis()->SetRangeUser(0.6,10.);
  hFinalFeeddown->Draw();*/
}

//--------------------------------------------------------------------------------------------------------

Double_t LevyTsallis_Func(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */
  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t n = p[1];
  Double_t C = p[2];
  Double_t norm = p[3];

  Double_t part1 = (n - 1.) * (n - 2.);
  Double_t part2 = n * C * (n * C + mass * (n - 2.));
  Double_t part3 = part1 / part2;
  Double_t part4 = 1. + (mt - mass) / n / C;
  Double_t part5 = TMath::Power(part4, -n);
  return pt * norm * part3 * part5;
}

TF1 * LevyTsallis(const Char_t *name, Double_t mass, Double_t n = 5., Double_t C = 0.001, Double_t norm = 0.01)
{
  
  TF1 *fLevyTsallis = new TF1(name, LevyTsallis_Func, 0., 10., 4);
  fLevyTsallis->SetParameters(mass, n, C, norm);
  fLevyTsallis->SetParNames("mass", "n", "C", "norm");
  fLevyTsallis->FixParameter(0, mass);
  fLevyTsallis->SetParLimits(1, 1.e-3, 1.e3);
  fLevyTsallis->SetParLimits(2, 1.e-4, 1.e4);
  fLevyTsallis->SetParLimits(3, 1.e-4, 1.e4);
  return fLevyTsallis;
}
  
//---------------------------------------------------------------------------------------------------------------

double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr ) {
    //Error in a Ratio
    if(B!=0) {
        Double_t errorfromtop = Aerr*Aerr / (B*B) ;
        Double_t errorfrombottom = ((A*A)/(B*B*B*B)) * Berr * Berr;
        return TMath::Sqrt( errorfromtop + errorfrombottom );
    }
    return 1.;
}
