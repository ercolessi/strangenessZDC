#ifndef AliAnalysisTaskLeading_h
#define AliAnalysisTaskLeading_h

class AliESDEvent;
class AliAnalysisFilter;
class TDatabasePDG;
class TTree;
class AliMCEvent;
class TH1D;
class TH2D;

#include "AliAnalysisTaskSE.h"
#include "AliMultSelection.h"
#include "AliMultSelectionTask.h"

class AliAnalysisTaskLeading : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskLeading(const char *name = "AliAnalysisTaskLeading"); 
  virtual ~AliAnalysisTaskLeading() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  Int_t GetNList(){return fNlist;};
  Int_t GetNTh1d(){return (fNth1d+fNth1dMC);};
  Int_t GetNTh2d(){return fNth2d;};
  Int_t GetNTree(){return fNtree;};
  const char *GetNameList(Int_t iList) const {if(iList >= 0 && iList < fNlist) return fNameList[iList]; else return "NoList";}

  //Other methods
  void SetMC(Bool_t mcsw=kTRUE){fIsMC=mcsw;SetMCdata(mcsw);} // set MC switch
  void SetMCdata(Bool_t mcsw=kTRUE); // set MC switch w/o kinematics

 private: 
  Bool_t fIsMC,fIsMCdata; // switch if MC data
  Bool_t fStart; //! switch for Exec initialization
  TDatabasePDG *fDB; //!
  Int_t fNtrack; //!
  Float_t fPV;

  // TList Output
  static const Int_t fNlist = 1;
  TList *fList[fNlist]; //!
  TList *fLexample; //!
  char fNameList[fNlist][40]; //!
  TList *f; //!

  // histos 1D
  static const Int_t fNth1d = 2;
  TH1D *fHth1d[fNth1d]; //!
  TH1D *fHTrigger;//!
  TH1D *fHpv;//!

  static const Int_t fNth1dMC = 0;
  //  TH1D *fHth1dMC[fNth1dMC]; //!
  //  TH1D *fHth1dMC[1]; //!

  // histos 2D
  static const Int_t fNth2d = 1;
  TH2D *fHth2d[fNth2d]; //!
  TH2D *fHTriggerMult; //!
  
  // TTree
  static const Int_t fNtree = 1;
  TTree *fTree[fNtree]; //!
  TTree *fTevents; //!
  
  // ZDC
  Int_t fRun; //!
  Float_t fZDCN1,fZDCN2,fZDCP1,fZDCP2; //!
  Float_t fadcZDCN1[5],fadcZDCN2[5],fadcZDCP1[5],fadcZDCP2[5]; //!
  Float_t fEZDC_EM,fEZDC_EM1; //!
  Float_t fN1pz,fN1X,fN1Y,fN1Z; //!
  Float_t fP1pz,fP1X,fP1Y,fP1Z; //!
  Float_t fN2pz,fN2X,fN2Y,fN2Z; //!
  Float_t fP2pz,fP2X,fP2Y,fP2Z; //!
  Int_t fN1pdg,fP1pdg,fN2pdg,fP2pdg; //!
  Int_t fN1hits,fN2hits,fP1hits,fP2hits; //!

  //Calib
  Float_t lV0M;//!
  Float_t lV0A;//!
  Float_t lV0C;//!
  Float_t lCL0;//!
  Float_t lCL1;//!
  Float_t lSPDClusters;//!
  Float_t lSPDTracklets;//!
  Float_t lZNA;//!
  Float_t lZNC;//!
  Float_t lZNApp;//!
  Float_t lZNCpp;//!
  Float_t fRefMult08;//!
  Float_t fRefMult05;//!
  Float_t fSPDClusters;//!
  Float_t fSPDTracklets;//!
  Float_t fV0MEstimator;//!
  Float_t fV0AEstimator;//!
  Float_t fV0CEstimator;//!
  Float_t AbsRefMult08;//!
  Float_t AbsRefMult05;//!
  Float_t AbsSPDClusters;//!
  Float_t AbsSPDTracklets;//!
  Float_t AbsV0MEstimator;//!
  Float_t AbsV0AEstimator;//!
  Float_t AbsV0CEstimator;//!


  AliESDEvent *fESD;    //! ESD object

  Int_t GetMotherArrays(AliMCEvent *mcEvent,Int_t iparticle,Int_t motherPdg[],Int_t motherLabel[]);

  AliAnalysisTaskLeading(const AliAnalysisTaskLeading&); // not implemented
  AliAnalysisTaskLeading& operator=(const AliAnalysisTaskLeading&); // not implemented
  
  ClassDef(AliAnalysisTaskLeading, 1); // example of analysis
};

#endif
