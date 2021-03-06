
////////////////////////////////////////////////////////////
//
// Default macro for calibrating minimum bias pp data.
//
////////////////////////////////////////////////////////////

#include "AliMultEstimator.h"
#include "AliMultSelection.h"
#include "AliMultSelectionCuts.h"

#include "AliMultSelectionCalibrator.h"
#include "Rtypes.h"
#include "TString.h"
#include "TSystem.h"

void CalibratePeriodPP(TString lPeriodName = "LHC15f", TString lWhichData = "pt2MB", Long_t lRunToUseAsDefault = 0)
{

  //Load ALICE stuff
  TString gLibs[] = { "STEER", "ANALYSIS", "ANALYSISalice", "ANALYSIScalib", "OADB" };
  TString thislib = "lib";
  for (Int_t ilib = 0; ilib < 5; ilib++) {
    thislib = "lib";
    thislib.Append(gLibs[ilib].Data());
    cout << "Will load " << thislib.Data() << endl;
    gSystem->Load(thislib.Data());
  }
  gSystem->SetIncludePath("-I$ROOTSYS/include  -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");

  cout << "Alive! " << endl;

  //All fine, let's try the calibrator
  AliMultSelectionCalibrator* lCalib = new AliMultSelectionCalibrator("lCalib");

  lCalib->SetRunToUseAsDefault(lRunToUseAsDefault);

  //============================================================
  // --- Definition of Boundaries ---
  //============================================================

  //Set Adaptive Percentile Boundaries, adjust if finer selection desired
  Double_t lDesiredBoundaries[1000];
  Long_t lNDesiredBoundaries = 0;
  lDesiredBoundaries[0] = 100;
  //From Low To High Multiplicity
  for (Int_t ib = 1; ib < 91; ib++) {
    lNDesiredBoundaries++;
    lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries - 1] - 1.0;
  }
  for (Int_t ib = 1; ib < 91; ib++) {
    lNDesiredBoundaries++;
    lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries - 1] - 0.1;
  }
  for (Int_t ib = 1; ib < 91; ib++) {
    lNDesiredBoundaries++;
    lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries - 1] - 0.01;
  }
  for (Int_t ib = 1; ib < 101; ib++) {
    lNDesiredBoundaries++;
    lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries - 1] - 0.001;
  }
  lNDesiredBoundaries++;
  lDesiredBoundaries[lNDesiredBoundaries] = 0;

  lCalib->SetBoundaries(lNDesiredBoundaries, lDesiredBoundaries);
  cout << "Boundaries set. Will attempt calibration now... " << endl;

  if (lPeriodName.Contains("LHC10h")) {
    cout << "Setting event selection criteria for Pb-Pb..." << endl;
    lCalib->GetEventCuts()->SetVzCut(10.0);
    lCalib->GetEventCuts()->SetTriggerCut(kTRUE);
    lCalib->GetEventCuts()->SetINELgtZEROCut(kFALSE);
    lCalib->GetEventCuts()->SetTrackletsVsClustersCut(kFALSE);
    lCalib->GetEventCuts()->SetRejectPileupInMultBinsCut(kFALSE);
    lCalib->GetEventCuts()->SetVertexConsistencyCut(kFALSE);
    lCalib->GetEventCuts()->SetNonZeroNContribs(kTRUE);
  }

  if (lPeriodName.Contains("LHC15m") || lPeriodName.Contains("LHC15o")) {
    cout << "Setting event selection criteria for Pb-Pb..." << endl;
    lCalib->GetEventCuts()->SetVzCut(10.0);
    lCalib->GetEventCuts()->SetTriggerCut(kTRUE);
    lCalib->GetEventCuts()->SetINELgtZEROCut(kFALSE);
    lCalib->GetEventCuts()->SetTrackletsVsClustersCut(kFALSE);
    lCalib->GetEventCuts()->SetRejectPileupInMultBinsCut(kFALSE);
    lCalib->GetEventCuts()->SetVertexConsistencyCut(kFALSE);
    lCalib->GetEventCuts()->SetNonZeroNContribs(kTRUE);
  }

  //Additional selections for pp: incompleteDAQ and asymmetric vzero

  //============================================================
  // --- Definition of Input Variables ---
  //============================================================

  lCalib->SetupStandardInput();

  //Changes in new version: create AliMultSelection here
  AliMultSelection* lMultSel = new AliMultSelection();
  lCalib->SetMultSelection(lMultSel);

  //============================================================
  // --- Definition of Estimators ---
  //============================================================

  Double_t lDefaultV0MAnchor = 0;
  Double_t lDefaultV0MPercentile = 0;

  Double_t lDefaultCL0Anchor = 0;
  Double_t lDefaultCL0Percentile = 0;

  Double_t lDefaultCL1Anchor = 0;
  Double_t lDefaultCL1Percentile = 0;

  if (lPeriodName.Contains("LHC15o")) {
    lDefaultV0MAnchor = 133.5;
    lDefaultV0MPercentile = 90.007;
    lDefaultCL0Anchor = 33.5;
    lDefaultCL0Percentile = 90.64;
    lDefaultCL1Anchor = 30.5;
    lDefaultCL1Percentile = 90.485;
  }
  if (lPeriodName.Contains("LHC15m")) {
    lDefaultV0MAnchor = 115.0;
    lDefaultV0MPercentile = 87.5;
    lDefaultCL0Anchor = 39.5;
    lDefaultCL0Percentile = 88.9;
    lDefaultCL1Anchor = 40.5;
    lDefaultCL1Percentile = 88.1;
  }

  AliMultEstimator* fEstV0M = new AliMultEstimator("V0M", "", "(fAmplitude_V0A)+(fAmplitude_V0C)");
  AliMultEstimator* fEstV0A = new AliMultEstimator("V0A", "", "(fAmplitude_V0A)");
  AliMultEstimator* fEstV0C = new AliMultEstimator("V0C", "", "(fAmplitude_V0C)");

  AliMultEstimator* fEstOnlineV0M = new AliMultEstimator("OnlineV0M", "", "(fAmplitude_OnlineV0A)+(fAmplitude_OnlineV0C)");
  AliMultEstimator* fEstOnlineV0A = new AliMultEstimator("OnlineV0A", "", "(fAmplitude_OnlineV0A)");
  AliMultEstimator* fEstOnlineV0C = new AliMultEstimator("OnlineV0C", "", "(fAmplitude_OnlineV0C)");

  AliMultEstimator* fEstADM = new AliMultEstimator("ADM", "", "(fMultiplicity_ADA)+(fMultiplicity_ADC)");
  AliMultEstimator* fEstADA = new AliMultEstimator("ADA", "", "(fMultiplicity_ADA)");
  AliMultEstimator* fEstADC = new AliMultEstimator("ADC", "", "(fMultiplicity_ADC)");

  //Integer estimators
  AliMultEstimator* fEstnSPDClusters = new AliMultEstimator("SPDClusters", "", "(fnSPDClusters)");
  fEstnSPDClusters->SetIsInteger(kTRUE);
  AliMultEstimator* fEstnSPDTracklets = new AliMultEstimator("SPDTracklets", "", "(fnTracklets)");
  fEstnSPDTracklets->SetIsInteger(kTRUE);
  AliMultEstimator* fEstRefMultEta5 = new AliMultEstimator("RefMult05", "", "(fRefMultEta5)");
  fEstRefMultEta5->SetIsInteger(kTRUE);
  AliMultEstimator* fEstRefMultEta8 = new AliMultEstimator("RefMult08", "", "(fRefMultEta8)");
  fEstRefMultEta8->SetIsInteger(kTRUE);

  //ZDC-based estimators
  /* From mail exchange with Cvetan
     -> ZNApp
    "-fZnaFired * fZnaTower + !fZnaFired * 1e6"

    -> ZNCpp
    "-fZncFired * fZncTower + !fZncFired * 1e6"

    -> ZNACpp
    "-0.89 * fZnaFired * fZnaTower - fZncFired * fZncTower + !fZnaFired * !fZncFired * 1e6"
    */

  //ZN
  AliMultEstimator* fEstZNApp = new AliMultEstimator("ZNApp", "", "-(fZnaFired) * (fZnaTower) + !(fZnaFired) * 1e6");
  AliMultEstimator* fEstZNCpp = new AliMultEstimator("ZNCpp", "", "-(fZncFired) * (fZncTower) + !(fZncFired) * 1e6");
  AliMultEstimator* fEstZNACpp = new AliMultEstimator("ZNACpp", "", "-0.89 * (fZnaFired) * (fZnaTower) - (fZncFired) * (fZncTower) + !(fZnaFired) * !(fZncFired) * 1e6");
  //ZP
  AliMultEstimator* fEstZPApp = new AliMultEstimator("ZPApp", "", "-(fZpaFired) * (fZpaTower) + !(fZpaFired) * 1e6");
  AliMultEstimator* fEstZPCpp = new AliMultEstimator("ZPCpp", "", "-(fZpcFired) * (fZpcTower) + !(fZpcFired) * 1e6");
  AliMultEstimator* fEstZPACpp = new AliMultEstimator("ZPACpp", "", "-0.89 * (fZpaFired) * (fZpaTower) - (fZpcFired) * (fZpcTower) + !(fZpaFired) * !(fZpcFired) * 1e6");
  //ZP+ZN
  AliMultEstimator* fEstZPNApp = new AliMultEstimator("ZPNApp", "", "-(fZpaFired) * (fZpaTower) + !(fZpaFired) * 1e6 - (fZnaFired) * (fZnaTower) + !(fZnaFired) * 1e6");
  AliMultEstimator* fEstZPNCpp = new AliMultEstimator("ZPNCpp", "", "-(fZpcFired) * (fZpcTower) + !(fZpcFired) * 1e6 - (fZncFired) * (fZncTower) + !(fZncFired) * 1e6");
  AliMultEstimator* fEstZPNACpp = new AliMultEstimator("ZPNACpp", "", "-0.89 * (fZpaFired) * (fZpaTower) - (fZpcFired) * (fZpcTower) + !(fZpaFired) * !(fZpcFired) * 1e6 - 0.89 * (fZnaFired) * (fZnaTower) - (fZncFired) * (fZncTower) + !(fZnaFired) * !(fZncFired) * 1e6");

  //ZDC percentile like estimators
  
  //ZN
  AliMultEstimator* fEstZNAppTower = new AliMultEstimator("ZNAppTower", "", "-(fZnaTower)");
  AliMultEstimator* fEstZNCppTower = new AliMultEstimator("ZNCppTower", "", "-(fZncTower)");
  AliMultEstimator* fEstZNACppTower = new AliMultEstimator("ZNACppTower", "", "- (fZnaTower) - (fZncTower) ");
  //ZP
  AliMultEstimator* fEstZPAppTower = new AliMultEstimator("ZPAppTower", "", "-(fZpaTower)");
  AliMultEstimator* fEstZPCppTower = new AliMultEstimator("ZPCppTower", "", "-(fZpcTower)");
  AliMultEstimator* fEstZPACppTower = new AliMultEstimator("ZPACppTower", "", "-(fZpaTower) - (fZpcTower)");
  //ZP+ZN
  AliMultEstimator* fEstZPNAppTower = new AliMultEstimator("ZPNAppTower", "", "- (fZpaTower) - (fZnaTower) ");
  AliMultEstimator* fEstZPNCppTower = new AliMultEstimator("ZPNCppTower", "", "- (fZpcTower) - (fZncTower)");
  AliMultEstimator* fEstZPNACppTower = new AliMultEstimator("ZPNACppTower", "", " - (fZpaTower) - (fZpcTower) - (fZnaTower) - (fZncTower)"); 


  //Universal: V0
  lCalib->GetMultSelection()->AddEstimator(fEstV0M);
  lCalib->GetMultSelection()->AddEstimator(fEstV0A);
  lCalib->GetMultSelection()->AddEstimator(fEstV0C);
  //Only do this in run 2, AD didn't exist in Run 1
  //Will also save space in the OADB for old datasets!
  lCalib->GetMultSelection()->AddEstimator(fEstOnlineV0M);
  lCalib->GetMultSelection()->AddEstimator(fEstOnlineV0A);
  lCalib->GetMultSelection()->AddEstimator(fEstOnlineV0C);
  lCalib->GetMultSelection()->AddEstimator(fEstADM);
  lCalib->GetMultSelection()->AddEstimator(fEstADA);
  lCalib->GetMultSelection()->AddEstimator(fEstADC);

  //Universal: Tracking, etc
  lCalib->GetMultSelection()->AddEstimator(fEstnSPDClusters);
  lCalib->GetMultSelection()->AddEstimator(fEstnSPDTracklets);
  lCalib->GetMultSelection()->AddEstimator(fEstRefMultEta5);
  lCalib->GetMultSelection()->AddEstimator(fEstRefMultEta8);

  lCalib->GetMultSelection()->AddEstimator(fEstZNApp);
  lCalib->GetMultSelection()->AddEstimator(fEstZNCpp);
  lCalib->GetMultSelection()->AddEstimator(fEstZNACpp);
  lCalib->GetMultSelection()->AddEstimator(fEstZPApp);
  lCalib->GetMultSelection()->AddEstimator(fEstZPCpp);
  lCalib->GetMultSelection()->AddEstimator(fEstZPACpp);
  lCalib->GetMultSelection()->AddEstimator(fEstZPNApp);
  lCalib->GetMultSelection()->AddEstimator(fEstZPNCpp);
  lCalib->GetMultSelection()->AddEstimator(fEstZPNACpp);
  lCalib->GetMultSelection()->AddEstimator(fEstZNAppTower);
  lCalib->GetMultSelection()->AddEstimator(fEstZNCppTower);
  lCalib->GetMultSelection()->AddEstimator(fEstZNACppTower);
  lCalib->GetMultSelection()->AddEstimator(fEstZPAppTower);
  lCalib->GetMultSelection()->AddEstimator(fEstZPCppTower);
  lCalib->GetMultSelection()->AddEstimator(fEstZPACppTower);
  lCalib->GetMultSelection()->AddEstimator(fEstZPNAppTower);
  lCalib->GetMultSelection()->AddEstimator(fEstZPNCppTower);
  lCalib->GetMultSelection()->AddEstimator(fEstZPNACppTower);

  //============================================================
  // --- Definition of Input/Output ---
  //============================================================

  lCalib->SetInputFile("AnalysisResults.root");
  lCalib->SetBufferFile(Form("buffer-%s-%s.root", lPeriodName.Data(), lWhichData.Data()));
  lCalib->SetOutputFile(Form("OADB-%s-%s.root", lPeriodName.Data(), lWhichData.Data()));
  lCalib->SetMaxEventsPerRun(1000000);
  lCalib->Calibrate();
}
