/***********************************************

  Lambda Run Analysis Module
  ----------------------

It constructs corrected Lambda, AntiLambda 

This version: 2020
Private version of Francesca Ercolessi 

--- David Dobrigkeit Chinellato
    daviddc@ifi.unicamp.br

***********************************************/

void runV0Analysis(
  TString lV0Type = "Lambda", 
  Double_t ptCutProt = 0.0, 
  TString filename="", 
  Bool_t doSystematics = kTRUE, 
  Bool_t kMultDepAnalysis = kTRUE, 
  TString multEstimator = "V0M", 
  Double_t multBin1=0.0, 
  Double_t multBin2=100., 
  TString fd = "UseMCRatio"){
  
  cout<<"Macro to test V0 analysis module"<<endl;

  cout<<"----------------------------------------------------"<<endl;
  cout<<"               V0 Analysis Macro "<<endl;
  cout<<"----------------------------------------------------"<<endl;
  cout<<endl;
  cout<<"----------------------------------------------------"<<endl;
  cout<<" ---> Compiling needed class, please wait... "<<endl;
  //Compile Macro
  Int_t workedornot = gSystem->CompileMacro("AliV0Module.cxx","-kfo");
  cout<<"----------------------------------------------------"<<endl;
  cout<<endl;
   if( workedornot == 0 ){ 
      cout<<"*************************************"<<endl;
      cout<<" AliV0Module.cxx compilation failed! "<<endl;
      cout<<"*************************************"<<endl;
      return;
   }

  //Load Class
  gSystem->Load("AliV0Module_cxx");

  //Initialize Analysis Object
  AliV0Module *v0 = new AliV0Module(lV0Type);

  //Binning
  Double_t ptbinlimits[] = {0.4, 0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.5, 2.9, 3.4, 4, 5, 6.5, 8, 10}; 
  Long_t ptbinnumb = sizeof(ptbinlimits)/sizeof(Double_t) - 1;
  //Set Pt Bins Used in Analysis
  v0->SetPtBinLimits( ptbinnumb, ptbinlimits );

  /// multiplicity
  v0->SetPerformMultiplicityStudy(kMultDepAnalysis);
  if(kMultDepAnalysis){
  v0->SetMultEstimator(multEstimator);
  v0->SetLowMultValue(multBin1);
  v0->SetHighMultValue(multBin2);
  }
  //Set Default Cuts - note: Particle dependent
  v0->SetDefaultCuts(); 
  v0->SetITSTOFrequest(-1);
  
  TString feeddownTreat = "UseMCRatio";
  //Feeddown treatment switch (applies only to Lambda, AntiLambda) 
  v0->SetFeeddownTreatment             (fd);

  //Set CINT1B/INEL Ratio for normalization
  v0->SetCINT1BoverINEL                ( 0.851 );

  v0->SetGeant3FlukaCorr(kFALSE);
  //if(lV0Type == "AntiLambda") v0->SetGeant3FlukaCorr(kTRUE);

  v0->SetUseIntegratedEff(kTRUE);
  v0->SetRealDataFile("~/Strangeness15f.root");
  v0->SetListName(filename);
  v0->SetMCDataFile("../Efficienze/MonteCarloMerged15f.root" );
  v0->SetFeedDownDataFile("../Efficienze/MonteCarloMerged15f.root");

  if(kMultDepAnalysis) {
    if(!v0->GetPosRap() && !v0->GetNegRap()){
      if(lV0Type.Contains("Lambda")){
        if(fd == "NoFD") v0->SetOutputFile( Form("Results-%s-13TeV-%s_%03.0f_%03.0f_ZDC_000_100%s.root",lV0Type.Data(), multEstimator.Data(),multBin1,multBin2,filename.Data()) );
        else v0->SetOutputFile( Form("Results-%s-13TeV-%s_%03.0f_%03.0f_ZDC_000_100_%sFD%s.root",lV0Type.Data(), multEstimator.Data(),multBin1,multBin2,fd.Data(),filename.Data()) );
      }   
    } 
    else{
      TString rFin = "pos"; 
      if(v0->GetNegRap()) rFin = "neg";
      if(lV0Type.Contains("Lambda")){
        if(fd == "NoFD") v0->SetOutputFile( Form("Results-%s-13TeV-%s_%03.0f_%03.0f_ZDC_000_100%s.root",lV0Type.Data(), multEstimator.Data(),multBin1,multBin2,rFin.Data()) );
        else v0->SetOutputFile( Form("Results-%s-13TeV-%s_%1.1f_%1.1f_ZDC_000_100_%sFD%s.root",lV0Type.Data(), multEstimator.Data(),multBin1,multBin2,fd.Data(),rFin.Data()) );
      } 
    }
  }
  else{
    if ( lV0Type == "Lambda"     ) v0->SetOutputFile( "Results-Lambda-13TeV.root"     );
    if ( lV0Type == "AntiLambda" ) v0->SetOutputFile( "Results-AntiLambda-13TeV.root");
  }

  //Run the analysis  
  v0->DoAnalysis();

  //Perform Systematics if demanded
  if ( doSystematics ){ 
    //------------------------------------------------------------
    // Topological Selection Variables Systematics
    TString lSystFile = "Results-Systematics"; 
    if(lV0Type == "Lambda"    ) lSystFile.Append("-Lambda-");
    if(lV0Type == "AntiLambda") lSystFile.Append("-AntiLambda-");
    if(lV0Type == "K0Short"   ) lSystFile.Append("-K0Short-");

    cout<<endl<<"---> Performing Systematics Studies: V0 Decay Radius"<<endl;
    v0->SetDefaultCuts(); 
    v0->SetCutV0Radius(0.300); v0->SetOutputFile( lSystFile + "V0Radius-1.root" ); v0->DoAnalysis(); 
    v0->SetCutV0Radius(0.400); v0->SetOutputFile( lSystFile + "V0Radius-2.root" ); v0->DoAnalysis(); 
    v0->SetCutV0Radius(0.600); v0->SetOutputFile( lSystFile + "V0Radius-3.root" ); v0->DoAnalysis(); 
    v0->SetCutV0Radius(0.700); v0->SetOutputFile( lSystFile + "V0Radius-4.root" ); v0->DoAnalysis(); 

    cout<<endl<<"---> Performing Systematics Studies: DCA Neg track to PV"<<endl;
    v0->SetDefaultCuts(); 
    v0->SetCutDCANegToPV(0.050); v0->SetOutputFile( lSystFile + "DCANegToPV-1.root" ); v0->DoAnalysis(); 
    v0->SetCutDCANegToPV(0.055); v0->SetOutputFile( lSystFile + "DCANegToPV-2.root" ); v0->DoAnalysis(); 
    v0->SetCutDCANegToPV(0.070); v0->SetOutputFile( lSystFile + "DCANegToPV-3.root" ); v0->DoAnalysis(); 
    v0->SetCutDCANegToPV(0.080); v0->SetOutputFile( lSystFile + "DCANegToPV-4.root" ); v0->DoAnalysis(); 

    cout<<endl<<"---> Performing Systematics Studies: DCA Pos track to PV"<<endl;
    v0->SetDefaultCuts(); 
    v0->SetCutDCAPosToPV(0.050); v0->SetOutputFile( lSystFile + "DCAPosToPV-1.root" ); v0->DoAnalysis(); 
    v0->SetCutDCAPosToPV(0.055); v0->SetOutputFile( lSystFile + "DCAPosToPV-2.root" ); v0->DoAnalysis(); 
    v0->SetCutDCAPosToPV(0.070); v0->SetOutputFile( lSystFile + "DCAPosToPV-3.root" ); v0->DoAnalysis(); 
    v0->SetCutDCAPosToPV(0.080); v0->SetOutputFile( lSystFile + "DCAPosToPV-4.root" ); v0->DoAnalysis(); 

    cout<<endl<<"---> Performing Systematics Studies: V0 Cos PA" <<endl;
    v0->SetDefaultCuts(); 
    v0->SetCutV0CosPA(0.993); v0->SetOutputFile( lSystFile + "V0CosPA-1.root" ); v0->DoAnalysis(); 
    v0->SetCutV0CosPA(0.994); v0->SetOutputFile( lSystFile + "V0CosPA-2.root" ); v0->DoAnalysis(); 
    v0->SetCutV0CosPA(0.996); v0->SetOutputFile( lSystFile + "V0CosPA-3.root" ); v0->DoAnalysis(); 
    v0->SetCutV0CosPA(0.997); v0->SetOutputFile( lSystFile + "V0CosPA-4.root" ); v0->DoAnalysis(); 
  
    cout<<endl<<"---> Performing Systematics Studies: DCA V0 Daughters"<<endl;
    v0->SetDefaultCuts(); 
    v0->SetCutDCAV0Daughters(1.50); v0->SetOutputFile( lSystFile + "DCAV0Daughters-1.root" ); v0->DoAnalysis(); 
    v0->SetCutDCAV0Daughters(1.25); v0->SetOutputFile( lSystFile + "DCAV0Daughters-2.root" ); v0->DoAnalysis(); 
    v0->SetCutDCAV0Daughters(0.75); v0->SetOutputFile( lSystFile + "DCAV0Daughters-3.root" ); v0->DoAnalysis(); 
    v0->SetCutDCAV0Daughters(0.50); v0->SetOutputFile( lSystFile + "DCAV0Daughters-4.root" ); v0->DoAnalysis();

    cout<<endl<<"---> Performing Systematics Studies: Proper Lifetime"<<endl;
    v0->SetDefaultCuts(); 
    v0->SetCutProperLifetime(40); v0->SetOutputFile( lSystFile + "ProperLifetime-1.root" ); v0->DoAnalysis(); 
    v0->SetCutProperLifetime(20); v0->SetOutputFile( lSystFile + "ProperLifetime-2.root" ); v0->DoAnalysis(); 
   
    cout<<endl<<"---> Performing Systematics Studies: Number Of CrossedRows("<<endl;
    v0->SetDefaultCuts(); 
    v0->SetCutLeastNumberOfCrossedRows(75); v0->SetOutputFile( lSystFile + "NumberOfCrossedRows-1.root" ); v0->DoAnalysis(); 
    v0->SetCutLeastNumberOfCrossedRows(80); v0->SetOutputFile( lSystFile + "NumberOfCrossedRows-2.root" ); v0->DoAnalysis(); 
    
    cout<<endl<<"---> Performing Systematics Studies: Number Of Crossed Rows Over Findable"<<endl;
    v0->SetDefaultCuts(); 
    v0->SetCutLeastNumberOfCrossedRowsOverFindable(0.95); v0->SetOutputFile( lSystFile + "NumberOfCrossedRowsOverFindable-1.root" ); v0->DoAnalysis(); 
    
    cout<<endl<<"---> Performing Systematics Studies: TPC PID NSigmas"<<endl;
    v0->SetDefaultCuts(); 
    v0->SetCutTPCPIDNSigmas(7); v0->SetOutputFile( lSystFile + "TPCPIDNSigmas-1.root" ); v0->DoAnalysis(); 
    v0->SetCutTPCPIDNSigmas(6); v0->SetOutputFile( lSystFile + "TPCPIDNSigmas-2.root" ); v0->DoAnalysis(); 
    v0->SetCutTPCPIDNSigmas(4); v0->SetOutputFile( lSystFile + "TPCPIDNSigmas-3.root" ); v0->DoAnalysis(); 
    
    cout<<endl<<"---> Performing Systematics Studies: Competing V0 Rejection"<<endl;
    v0->SetDefaultCuts(); 
    v0->SetCutCompetingV0Rejection(0); v0->SetOutputFile( lSystFile + "CompetingV0Rejection-1.root" ); v0->DoAnalysis(); 

    cout<<endl<<"---> Performing Systematics Studies: Sigma For Signal Extraction"<<endl;
    v0->SetDefaultCuts(); 
    v0->SetCutSigmaForSignalExtraction(7); v0->SetOutputFile( lSystFile + "SigmaForSignalExtraction-1.root" ); v0->DoAnalysis();
    v0->SetCutSigmaForSignalExtraction(5); v0->SetOutputFile( lSystFile + "SigmaForSignalExtraction-2.root" ); v0->DoAnalysis();
    v0->SetCutSigmaForSignalExtraction(4); v0->SetOutputFile( lSystFile + "SigmaForSignalExtraction-3.root" ); v0->DoAnalysis(); 
    //------------------------------------------------------------

    //other systematics: could be placed here in the future...    
  }

  return;
}
