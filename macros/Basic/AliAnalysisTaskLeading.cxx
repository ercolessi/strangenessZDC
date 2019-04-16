 #include "TChain.h"
#include "TTree.h"
#include "TLeaf.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisTaskLeading.h"

#include "AliAnalysisFilter.h"
#include "AliESDtrackCuts.h"

#include "TDatabasePDG.h"
#include "AliESDVertex.h"
#include "AliESDpid.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"

#include "TRandom.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"

#include "AliTOFcalibHisto.h"
#include "AliCentrality.h"
#include "AliTOFGeometry.h"

#include "TSpline.h"
#include "AliMultSelection.h"
#include "AliMultSelectionTask.h"

ClassImp(AliAnalysisTaskLeading)

//________________________________________________________________________
AliAnalysisTaskLeading::AliAnalysisTaskLeading(const char *name) 
:     AliAnalysisTaskSE(name), fESD(0), fIsMC(kFALSE),fIsMCdata(kFALSE), fStart(0), fNtrack(0), fPV(0), lV0M(300)

{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());

   for(Int_t i=0;i<fNlist;i++){
      DefineOutput(i, TList::Class());
   } 

   // List for PID performance
   sprintf(fNameList[0],"0-esempio");
   fLexample = new TList();
   fList[0] = fLexample; 
  
  // Output slot histos
  for(Int_t i=0;i<fNth1d+fNth1dMC;i++){
      DefineOutput(i+fNlist, TH1D::Class());
  }

  for(Int_t i=0;i<fNth2d;i++){
      DefineOutput(i+fNth1d+fNth1dMC+fNlist, TH2D::Class());
  }

  for(Int_t i=0;i<fNtree;i++){
      DefineOutput(i+fNth1d+fNth1dMC+fNth2d+fNlist, TTree::Class());
  }
}
//________________________________________________________________________
void AliAnalysisTaskLeading::SetMCdata(Bool_t mcsw){
  fIsMCdata=mcsw;
}
//________________________________________________________________________
void AliAnalysisTaskLeading::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  // Create histograms
  // Called once

  fDB  = TDatabasePDG::Instance();

  // create histos

  //1D
  fHTrigger = new TH1D("hTrigger","Trigger Input;;N",14,0,14);
  fHpv = new TH1D("hpv","Primary Vertex;z (cm);N",200,-20,20);

  // fill with the trigger label
  fHTrigger->Fill("events",0);
  fHTrigger->Fill("triggered",0);
  fHTrigger->Fill("SPD-0SMB",0);
  fHTrigger->Fill("SPD-0SM1",0);
  fHTrigger->Fill("SPD-0SM2",0);
  fHTrigger->Fill("TOF-0OM2",0);
  fHTrigger->Fill("TOF-0OM3",0);
  fHTrigger->Fill("TOF-0OUP",0);
  fHTrigger->Fill("TOF-0OIN",0);
  fHTrigger->Fill("VZERO-0VBA",0);
  fHTrigger->Fill("VZERO-0VBC",0);
  fHTrigger->Fill("*CMBAC-B-NOPF-ALL",0);
  fHTrigger->Fill("*CMS2C-B-NOPF-ALL",0);
  fHTrigger->Fill("*CMS2A-B-NOPF-ALL",0);
  fHTrigger->Fill("*C0OM2-B-NOPF-ALL",0);
  fHTrigger->Fill("AtLeast1TPCtr",0);

  fHth1d[0] = fHTrigger;
  fHth1d[1] = fHpv;

  for(Int_t i=0;i<fNth1d;i++){
    fHth1d[i]->Sumw2();
    fHth1d[i]->GetXaxis()->SetNdivisions(409);
    fHth1d[i]->GetYaxis()->SetNdivisions(409);
    fHth1d[i]->GetXaxis()->SetTitleSize(0.05);
    fHth1d[i]->GetYaxis()->SetTitleSize(0.05);
  }

  /*
  for(Int_t i=0;i<fNth1dMC;i++){
    fHth1dMC[i]->Sumw2();
    fHth1dMC[i]->GetXaxis()->SetNdivisions(409);
    fHth1dMC[i]->GetYaxis()->SetNdivisions(409);
    fHth1dMC[i]->GetXaxis()->SetTitleSize(0.05);
    fHth1dMC[i]->GetYaxis()->SetTitleSize(0.05);
  }
  */

  //2D
  fHTriggerMult = new TH2D("hTriggerMult","Multiplicity vs. Trigger Input;;N_{tracks}^{TPC};N",17,0,17,2000,0,2000);

  // fill with the trigger label
  fHTriggerMult->Fill("triggered",0.,0);
  fHTriggerMult->Fill("SPD-0SMB",0.,0);
  fHTriggerMult->Fill("SPD-0SM1",0.,0);
  fHTriggerMult->Fill("SPD-0SM2",0.,0);
  fHTriggerMult->Fill("TOF-0OM2",0.,0);
  fHTriggerMult->Fill("TOF-0OM3",0.,0);
  fHTriggerMult->Fill("TOF-0OUP",0.,0);
  fHTriggerMult->Fill("TOF-0OIN",0.,0);
  fHTriggerMult->Fill("VZERO-0VBA",0.,0);
  fHTriggerMult->Fill("VZERO-0VBC",0.,0);
  fHTriggerMult->Fill("*CMBAC-B-NOPF-ALL",0.,0);
  fHTriggerMult->Fill("*CMS2C-B-NOPF-ALL",0.,0);
  fHTriggerMult->Fill("*CMS2A-B-NOPF-ALL",0.,0);
  fHTriggerMult->Fill("*C0OM2-B-NOPF-ALL",0.,0);

  fHth2d[0] = fHTriggerMult;

  for(Int_t i=0;i<fNth2d;i++){
    fHth2d[i]->Sumw2();
    fHth2d[i]->GetXaxis()->SetNdivisions(409);
    fHth2d[i]->GetYaxis()->SetNdivisions(409);
    fHth2d[i]->GetXaxis()->SetTitleSize(0.05);
    fHth2d[i]->GetYaxis()->SetTitleSize(0.05);
  }

  //tree
  fTevents = new TTree("tEvents","tEvents");
//   fTevents->Branch("ZDCN1",&fZDCN1,"ZDCN1/F");
//   fTevents->Branch("ZDCN2",&fZDCN2,"ZDCN2/F");
//   fTevents->Branch("ZDCP1",&fZDCP1,"ZDCP1/F");
//   fTevents->Branch("ZDCP2",&fZDCP2,"ZDCP2/F");
  fTevents->Branch("adcZDCN1",fadcZDCN1,"adcZDCN1[5]/F");
  fTevents->Branch("adcZDCN2",fadcZDCN2,"adcZDCN2[5]/F");
  fTevents->Branch("adcZDCP1",fadcZDCP1,"adcZDCP1[5]/F");
  fTevents->Branch("adcZDCP2",fadcZDCP2,"adcZDCP2[5]/F");
  fTevents->Branch("adcEM",&fEZDC_EM,"adcEM/F");
  fTevents->Branch("adcEM1",&fEZDC_EM1,"adcEM1/F");
  fTevents->Branch("nrun",&fRun,"nrun/I");
  fTevents->Branch("ntrack",&fNtrack,"ntrack/I");
  fTevents->Branch("V0M",&lV0M);//,"V0M/F");
  fTevents->Branch("V0A",&lV0A,"V0A/F");
  fTevents->Branch("V0C",&lV0C,"V0C/F");
  fTevents->Branch("CL0",&lCL0,"CL0/F");
  fTevents->Branch("CL1",&lCL1,"CL1/F");
  fTevents->Branch("SPDClusters",&lSPDClusters,"SPDClusters/F");
  fTevents->Branch("SPDTracklets",&lSPDTracklets,"SPDTracklets/F");
  fTevents->Branch("ZNA",&lZNA,"ZNA/F");
  fTevents->Branch("ZNC",&lZNC,"ZNC/F");
  fTevents->Branch("ZNApp",&lZNApp,"ZNApp/F");
  fTevents->Branch("ZNCpp",&lZNCpp,"ZNCpp/F");
  fTevents->Branch("RefMult08_p",&fRefMult08,"RefMult08_p/F");
  fTevents->Branch("RefMult05_p",&fRefMult05,"RefMult05_p/F");
  fTevents->Branch("SPDClusters_p",&fSPDClusters,"SPDClusters_p/F");
  fTevents->Branch("SPDTracklets_p",&fSPDTracklets,"SPDTracklets_p/F");
  fTevents->Branch("V0MEstimator_p",&fV0MEstimator,"V0MEstimator_p/F");
  fTevents->Branch("V0AEstimator_p",&fV0AEstimator,"V0AEstimator_p/F");
  fTevents->Branch("V0CEstimator_p",&fV0CEstimator,"V0CEstimator_p/F");
  fTevents->Branch("RefMult08_abs",&AbsRefMult08,"RefMult08_abs/F");
  fTevents->Branch("RefMult05_abs",&AbsRefMult05,"RefMult05_abs/F");
  fTevents->Branch("SPDClusters_abs",&AbsSPDClusters,"SPDClusters_abs/F");
  fTevents->Branch("SPDTracklets_abs",&AbsSPDTracklets,"SPDTracklets_abs/F");
  fTevents->Branch("V0MEstimator_abs",&AbsV0MEstimator,"V0MEstimator_abs/F");
  fTevents->Branch("V0AEstimator_abs",&AbsV0AEstimator,"V0AEstimator_abs/F");
  fTevents->Branch("V0CEstimator_abs",&AbsV0CEstimator,"V0CEstimator_abs/F");


  if(fIsMC){
    fTevents->Branch("N1pdg",&fN1pdg,"N1pdg/I");
    fTevents->Branch("N1pz",&fN1pz,"N1pz/F");
    fTevents->Branch("N1X",&fN1X,"N1X/F");
    fTevents->Branch("N1Y",&fN1Y,"N1Y/F");
    fTevents->Branch("N1Z",&fN1Z,"N1Z/F");
    fTevents->Branch("P1pdg",&fP1pdg,"P1pdg/I");
    fTevents->Branch("P1pz",&fP1pz,"P1pz/F");
    fTevents->Branch("P1X",&fP1X,"P1X/F");
    fTevents->Branch("P1Y",&fP1Y,"P1Y/F");
    fTevents->Branch("P1Z",&fP1Z,"P1Z/F");
    fTevents->Branch("N2pdg",&fN2pdg,"N2pdg/I");
    fTevents->Branch("N2pz",&fN2pz,"N2pz/F");
    fTevents->Branch("N2X",&fN2X,"N2X/F");
    fTevents->Branch("N2Y",&fN2Y,"N2Y/F");
    fTevents->Branch("N2Z",&fN2Z,"N2Z/F");
    fTevents->Branch("P2pdg",&fP2pdg,"P2pdg/I");
    fTevents->Branch("P2pz",&fP2pz,"P2pz/F");
    fTevents->Branch("P2X",&fP2X,"P2X/F");
    fTevents->Branch("P2Y",&fP2Y,"P2Y/F");
    fTevents->Branch("P2Z",&fP2Z,"P2Z/F");
    fTevents->Branch("N1hits",&fN1hits,"N1hits/I");
    fTevents->Branch("N2hits",&fN2hits,"N2hits/I");
    fTevents->Branch("P1hits",&fP1hits,"P1hits/I");
    fTevents->Branch("P2hits",&fP2hits,"P2hits/I");
  }

  fTree[0] = fTevents;

  // TList example
  //  fList[0]->Add(fHBetheBloch); // example to add an object to the list
  
  for(Int_t i=0;i < fNlist;i++){
    fList[i]->SetOwner();
    PostData(i,fList[i]);
  }
  
  for(Int_t i=0;i < fNth1d;i++){
    PostData(i+fNlist,fHth1d[i]);
  }
  /*
  if(fIsMC){
    for(Int_t i=0;i < fNth1dMC;i++){
      PostData(i+fNth1d+fNlist,fHth1dMC[i]);
    }
  }
  */
  
  for(Int_t i=0;i < fNth2d;i++){
    PostData(i+fNth1d+fNth1dMC+fNlist,fHth2d[i]);
  }
  for(Int_t i=0;i < fNtree;i++){
    PostData(i+fNth1d+fNth1dMC+fNth2d+fNlist,fTree[i]);
  }
}
//________________________________________________________________________
void AliAnalysisTaskLeading::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());

  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }

  // TOF-T0 maker inizialization
  if(!fStart){
    fStart=kTRUE;
  }

  // flag
  Bool_t isSelected = ((AliInputEventHandler*)((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler()))->IsEventSelected();
  
  // return if no selection
  if(!isSelected){
    return;
  } 

  fRun = fESD->GetRunNumber();

  // ZDC info
  fZDCN1 = fESD->GetZDCN1Energy();
  fZDCN2 = fESD->GetZDCN2Energy();
  fZDCP1 = fESD->GetZDCP1Energy();
  fZDCP2 = fESD->GetZDCP2Energy();

  fEZDC_EM = fESD->GetESDZDC()->GetZDCEMEnergy(0);
  fEZDC_EM1 = fESD->GetESDZDC()->GetZDCEMEnergy(1);

  // read adc of ZDC
  const Double_t *aZDCN1 = fESD->GetESDZDC()->GetZN1TowerEnergy();
  fadcZDCN1[0] = aZDCN1[0];
  fadcZDCN1[1] = aZDCN1[1];
  fadcZDCN1[2] = aZDCN1[2];
  fadcZDCN1[3] = aZDCN1[3];
  fadcZDCN1[4] = aZDCN1[4];
  const Double_t *aZDCN2 = fESD->GetESDZDC()->GetZN2TowerEnergy();
  fadcZDCN2[0] = aZDCN2[0];
  fadcZDCN2[1] = aZDCN2[1];
  fadcZDCN2[2] = aZDCN2[2];
  fadcZDCN2[3] = aZDCN2[3];
  fadcZDCN2[4] = aZDCN2[4];
  const Double_t *aZDCP1 = fESD->GetESDZDC()->GetZP1TowerEnergy();
  fadcZDCP1[0] = aZDCP1[0];
  fadcZDCP1[1] = aZDCP1[1];
  fadcZDCP1[2] = aZDCP1[2];
  fadcZDCP1[3] = aZDCP1[3];
  fadcZDCP1[4] = aZDCP1[4];
  const Double_t *aZDCP2 = fESD->GetESDZDC()->GetZP2TowerEnergy();
  fadcZDCP2[0] = aZDCP2[0];
  fadcZDCP2[1] = aZDCP2[1];
  fadcZDCP2[2] = aZDCP2[2];
  fadcZDCP2[3] = aZDCP2[3];
  fadcZDCP2[4] = aZDCP2[4];

   
  AliMultSelection* MultSelection = 0x0;
  MultSelection = (AliMultSelection*)fESD->FindListObject("MultSelection");
  //MultSelection->GetEstimatorList()->ls();
  //Printf("%f, p= %f", MultSelection->GetEstimator("RefMult08")->GetValue(), MultSelection->GetEstimator("RefMult08")->GetPercentile());
  if (!MultSelection) {
    //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
    AliWarning("AliMultSelection object not found!");
  } else {
    
    lV0M = MultSelection->GetMultiplicityPercentile("V0M", true);
    lV0A = MultSelection->GetMultiplicityPercentile("V0A", true);
    lV0C = MultSelection->GetMultiplicityPercentile("V0C", true);
    lCL0 = MultSelection->GetMultiplicityPercentile("CL0", true);
    lCL1 = MultSelection->GetMultiplicityPercentile("CL1", true);
    lSPDClusters  = MultSelection->GetMultiplicityPercentile("SPDClusters", true );
    lSPDTracklets = MultSelection->GetMultiplicityPercentile("SPDTracklets", true);
    lZNA = MultSelection->GetMultiplicityPercentile("ZNA", true);
    lZNC = MultSelection->GetMultiplicityPercentile("ZNC", true);
    lZNApp = MultSelection->GetMultiplicityPercentile("ZNApp", true);
    lZNCpp = MultSelection->GetMultiplicityPercentile("ZNCpp", true);
 
    fRefMult08 = MultSelection->GetEstimator("RefMult08")->GetPercentile();
    fRefMult05 = MultSelection->GetEstimator("RefMult05")->GetPercentile();
    fSPDClusters = MultSelection->GetEstimator("SPDClusters")->GetPercentile();
    fSPDTracklets = MultSelection->GetEstimator("SPDTracklets")->GetPercentile();
    fV0MEstimator = MultSelection->GetEstimator("V0M")->GetPercentile();
    fV0AEstimator = MultSelection->GetEstimator("V0A")->GetPercentile();
    fV0CEstimator = MultSelection->GetEstimator("V0C")->GetPercentile();
    
    AbsRefMult08 = MultSelection->GetEstimator("RefMult08")->GetValue();
    AbsRefMult05 = MultSelection->GetEstimator("RefMult05")->GetValue();
    AbsSPDClusters = MultSelection->GetEstimator("SPDClusters")->GetValue();
    AbsSPDTracklets = MultSelection->GetEstimator("SPDTracklets")->GetValue();
    AbsV0MEstimator = MultSelection->GetEstimator("V0M")->GetValue();
    AbsV0AEstimator = MultSelection->GetEstimator("V0A")->GetValue();
    AbsV0CEstimator = MultSelection->GetEstimator("V0C")->GetValue();
    
  }

  // loop over TPC tracks
  fNtrack = fESD->GetNumberOfTracks();
  /*
  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) { // ESD track loop
   AliESDtrack* track = fESD->GetTrack(iTracks);
   if (!track) {
     Printf("ERROR: Could not receive track %d", iTracks);
     continue;
   }
  }
  */

  fPV = fESD->GetPrimaryVertex()->GetZ();

  // Trigger info
  Bool_t ktrigger = kFALSE;
  AliESDHeader* esdH = fESD->GetHeader(); // taken the header from AliESDEvent
  
  if(fESD->IsPileupFromSPDInMultBins()){
      fHTrigger->Fill("PILE-UP-MULT",1);
  }
  if(fESD->IsPileupFromSPD()){
      fHTrigger->Fill("PILE-UP",1);
  }
  if(esdH->IsTriggerInputFired("0VBA")){
      fHTrigger->Fill("VZERO-0VBC",1);
  }
  if(esdH->IsTriggerInputFired("0VBC")){ 
      fHTrigger->Fill("VZERO-0VBA",1);
  }
  if(esdH->IsTriggerInputFired("0SMB")){ 
      fHTrigger->Fill("SPD-0SMB",1);
  }
  if(esdH->IsTriggerInputFired("0SM1")){ 
      fHTrigger->Fill("SPD-0SM1",1);
  }
  if(esdH->IsTriggerInputFired("0SM2")){ 
      fHTrigger->Fill("SPD-0SM2",1);
  }
  if(esdH->IsTriggerInputFired("0OM2")){ // thr >= 2 (input 19)
      fHTrigger->Fill("TOF-0OM2",1);
     fHTrigger->Fill("*C0OM2-B-NOPF-ALL",1);
     ktrigger=kTRUE;
  }
  if(esdH->IsTriggerInputFired("0OM3")){ // thr >= 3 (input 20)
      fHTrigger->Fill("TOF-0OM3",1);
  }
  
  if(esdH->IsTriggerInputFired("0OUP")){ // Ultra Pheripheral (input 21)
      fHTrigger->Fill("TOF-0OUP",1);
  }
  if(esdH->IsTriggerInputFired("0VBA") && esdH->IsTriggerInputFired("0VBC")){
    fHTrigger->Fill("*CMBAC-B-NOPF-ALL",1);
    ktrigger=kTRUE;
  }
  if(esdH->IsTriggerInputFired("0VBC") && esdH->IsTriggerInputFired("0SM2")){
    fHTrigger->Fill("*CMS2C-B-NOPF-ALL",1);
    ktrigger=kTRUE;
  }
  if(esdH->IsTriggerInputFired("0VBA") && esdH->IsTriggerInputFired("0SM2")){
    fHTrigger->Fill("*CMS2A-B-NOPF-ALL",1);
    ktrigger=kTRUE;
  }

  fHTrigger->Fill("events",1);
  if(ktrigger){
    fHTrigger->Fill("triggered",1);
    fHpv->Fill(fPV);
  }
    
  // if analyzing MC data ------------
  AliMCEvent* mcEvent=NULL;
  AliMCEventHandler* eventHandler=NULL;

  // reset some variables
  fN1pdg=0;
  fP1pdg=0;
  fN2pdg=0;
  fP2pdg=0;
  fN1pz=0;
  fN2pz=0;
  fP1pz=0;
  fP2pz=0;
  fN1X=0;
  fN1Y=0;
  fP1X=0;
  fP1Y=0;
  fN2X=0;
  fN2Y=0;
  fP2X=0;
  fP2Y=0;
  fN1Z=0;
  fN2Z=0;
  fP1Z=0;
  fP2Z=0;
  fN1hits=0;
  fN2hits=0;
  fP1hits=0;
  fP2hits=0;

  // take MC truth for ZDC from Track Refs
  if(fIsMC){
    eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!eventHandler) {
      Printf("ERROR: Could not retrieve MC event handler (%i)",fIsMC);
      return;
    }
    
    mcEvent = eventHandler->MCEvent();

    TTree *treeTR = eventHandler->TreeTR();

    Int_t nevTR = treeTR->GetEntries();
    for(Int_t i=0;i<nevTR;i++){
      treeTR->GetEvent(i);
      Float_t L = treeTR->GetLeaf("TrackReferences.fLength")->GetValue();

      if(L > 10300){
	Float_t pz = treeTR->GetLeaf("TrackReferences.fPz")->GetValue();
	Int_t itr = treeTR->GetLeaf("TrackReferences.fTrack")->GetValue();
	Float_t x = treeTR->GetLeaf("TrackReferences.fX")->GetValue()*10;
	Float_t y = treeTR->GetLeaf("TrackReferences.fY")->GetValue()*10;

	AliMCParticle *partMC = (AliMCParticle *) mcEvent->GetTrack(itr);
	TParticle *part = partMC->Particle();
	Float_t pzStack = part->Pz();
	Int_t pdg = part->GetPdgCode();
	Int_t charge = fDB->GetParticle(pdg)->Charge();

	if(y>-40 && y < 90 && x > 130 && x < 260 && L > 11370 && L < 11480){
	  if(pz > fP2pz){
	    fP2pz = pz;
	    fP2X = x;
	    fP2Y = y;
	    fP2pdg = pdg;
	    fP2Z = L;
	  }
	  else if(-pz > fP1pz){
	    fP1pz = -pz;
	    fP1X = x;
	    fP1Y = y;
	    fP1pdg = pdg;	    
	    fP1Z = L;
	  }
	  if(pz > 0) fP2hits++;
	  if(pz < 0) fP1hits++;
	}
	else if(y>-40 && y < 40 && x > -40 && x < 40 && L > 11300 && L < 11360){
	  if(pz > fN2pz){
	    fN2pz = pz;
	    fN2X = x;
	    fN2Y = y;
	    fN2pdg = pdg;
	    fN2Z = L;
	  }
	  else if(-pz > fN1pz){
	    fN1pz = -pz;
	    fN1X = x;
	    fN1Y = y;
	    fN1pdg = pdg;	    
	    fN1Z = L;
	  }
	  if(pz > 0) fN2hits++;
	  if(pz < 0) fN1hits++;
	}
      }      
    }
  }

  //--------------------------------- ZDC MC track Refs
  
  fHTriggerMult->Fill("triggered",fNtrack,1);
  if(fESD->IsPileupFromSPDInMultBins()){
      fHTriggerMult->Fill("PILE-UP-MULT",fNtrack,1);
  }
  if(fESD->IsPileupFromSPD()){
      fHTriggerMult->Fill("PILE-UP",fNtrack,1);
  }
  if(esdH->IsTriggerInputFired("0VBA")){
      fHTriggerMult->Fill("VZERO-0VBC",fNtrack,1);
  }
  if(esdH->IsTriggerInputFired("0VBC")){ 
      fHTriggerMult->Fill("VZERO-0VBA",fNtrack,1);
  }
  if(esdH->IsTriggerInputFired("0SMB")){ 
      fHTriggerMult->Fill("SPD-0SMB",fNtrack,1);
  }
  if(esdH->IsTriggerInputFired("0SM1")){ 
      fHTriggerMult->Fill("SPD-0SM1",fNtrack,1);
  }
  if(esdH->IsTriggerInputFired("0SM2")){ 
      fHTriggerMult->Fill("SPD-0SM2",fNtrack,1);
  }
  if(esdH->IsTriggerInputFired("0OM2")){ // thr >= 2 (input 19)
      fHTriggerMult->Fill("TOF-0OM2",fNtrack,1);
      fHTriggerMult->Fill("*C0OM2-B-NOPF-ALL",fNtrack,1);
  }
  if(esdH->IsTriggerInputFired("0OM3")){ // thr >= 3 (input 20)
      fHTriggerMult->Fill("TOF-0OM3",fNtrack,1);
  }
  
  if(esdH->IsTriggerInputFired("0OUP")){ // Ultra Pheripheral (input 21)
      fHTriggerMult->Fill("TOF-0OUP",fNtrack,1);
  }

  if(esdH->IsTriggerInputFired("0VBA") && esdH->IsTriggerInputFired("0VBC")){
    fHTriggerMult->Fill("*CMBAC-B-NOPF-ALL",fNtrack,1);
  }
  if(esdH->IsTriggerInputFired("0VBC") && esdH->IsTriggerInputFired("0SM2")){
    fHTriggerMult->Fill("*CMS2C-B-NOPF-ALL",fNtrack,1);
  }
  if(esdH->IsTriggerInputFired("0VBA") && esdH->IsTriggerInputFired("0SM2")){
    fHTriggerMult->Fill("*CMS2A-B-NOPF-ALL",fNtrack,1);
  }


  if(fNtrack){
    fHTrigger->Fill("AtLeast1TPCtr",1);
  }

  fTevents->Fill();

  for(Int_t i=0;i < fNlist;i++){
    PostData(i,fList[i]);
  }

  for(Int_t i=0;i < fNth1d;i++){
      PostData(i+fNlist,fHth1d[i]);
  }

  /*
  if(fIsMC){
    for(Int_t i=0;i < fNth1dMC;i++){
      PostData(i+fNth1d+fNlist,fHth1dMC[i]);
    }
  }
  */
  
  for(Int_t i=0;i < fNth2d;i++){
      PostData(i+fNth1d+fNth1dMC+fNlist,fHth2d[i]);
  }
  for(Int_t i=0;i < fNtree;i++){
      PostData(i+fNth1d+fNth1dMC+fNth2d+fNlist,fTree[i]);
  }


 
}      
//________________________________________________________________________
void AliAnalysisTaskLeading::Terminate(Option_t *) 
{
  for(Int_t i=0;i < fNlist;i++){
      fList[i]=dynamic_cast<TList*> (GetOutputData(i));
      if (!fList[i]) {
	  Printf("ERROR: TList n.%i not available",i+1);
	  return;   
      } 
  }
  for(Int_t i=0;i < fNth1d;i++){
      fHth1d[i]=dynamic_cast<TH1D*> (GetOutputData(i+fNlist));
      if (!fHth1d[i]) {
	  Printf("ERROR: TH1D n.%i not available",i+1);
	  return;   
      } 
  }
  /*
  if(fIsMC){
      for(Int_t i=0;i < fNth1dMC;i++){
	  fHth1dMC[i]=dynamic_cast<TH1D*> (GetOutputData(i+fNth1d+fNlist));
	  if (!fHth1dMC[i]) {
	      Printf("ERROR: TH1D MC n.%i not available",i+1);
	      return;   
	  } 
      }
  }
  */

  for(Int_t i=0;i < fNth2d;i++){
      fHth2d[i]=dynamic_cast<TH2D*> (GetOutputData(i+fNth1d+fNth1dMC+fNlist));
      if (!fHth1d[i]) {
	  Printf("ERROR: TH2D n.%i not available",i+1);
	  return;   
      } 
  }
  for(Int_t i=0;i < fNtree;i++){
      fTree[i]=dynamic_cast<TTree*> (GetOutputData(i+fNth1d+fNth1dMC+fNth2d+fNlist));
      if (!fTree[i]) {
	  Printf("ERROR: TTree n.%i not available",i+1);
	  return;   
      } 
  }

  system("touch ok.job");

}
//________________________________________________________________________
Int_t AliAnalysisTaskLeading::GetMotherArrays(AliMCEvent *mcEvent,Int_t iparticle,Int_t motherPdg[],Int_t motherLabel[]){
  Int_t nm=0;
  const Int_t nmaxmother = 10;
  AliMCParticle *partMC = (AliMCParticle *) mcEvent->GetTrack(iparticle);
  TParticle *part = partMC->Particle();

  while(part->GetFirstMother() != -1 && nm < nmaxmother){
    iparticle = TMath::Abs(part->GetFirstMother());
    partMC = (AliMCParticle *) mcEvent->GetTrack(iparticle);
    part = partMC->Particle();
    motherPdg[nm] = part->GetPdgCode();
    motherLabel[nm] = iparticle;    
    nm++;
  }

  return nm++;
}
