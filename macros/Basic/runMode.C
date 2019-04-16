//How to run 

//testmode) gridtest=KTRUE

//gridmode) gridtest=KFALSE
//--1)run with runmode=full and setMergeViaJDL(kTRUE) to run on grid
//--2)run with runmode=terminate and SetMergeViaJDL(kTRUE) to merge
//--3)run with runmode=terminate and SetMergeViaJDL(kFALSE) to download in local
//    the merged root file

#if !defined (__CINT__) || defined (__CLING__)
#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisTaskLeading.h"
#include "AliMultSelection.h"
#include "AliMultSelectionTask.h"
#include "AddTaskPhysicsSelection.h"
#endif

void runMode()
  
{
  // if you need to access data remotely from the GRID
  Bool_t grid = 1;  
  // if you run on grid, specify test mode (kTRUE) or full grid model (kFALSE)
  Bool_t gridTest = kFALSE;
  // if the data are MC
  Bool_t isMCdata=kFALSE;
  Bool_t isMC=kFALSE; // if you have the kinematics information
  if(isMC) isMCdata=kTRUE;
  // if it is multiplicity calib. mode
  Bool_t isMultCalib = kFALSE;

  // Load common libraries
  gSystem->Load("libCore.so");  
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");   
   
  // include the path you need to compile and the library you need
  gSystem->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");   
  //Tell root where to look for headers
#if !defined (__CINT__) || defined (__CLING__)
  gInterpreter->ProcessLine(".include $ROOTSYS/include");
  gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
#else
  gROOT->ProcessLine(".include $ROOTSYS/include");
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
#endif

  // Create the analysis manager and ESD handler
  AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");
  AliESDInputHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);
  // add interface to MC if requested
  if(isMC){
    AliMCEventHandler* mcH = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mcH);
  }  

  //PhysicsSelection Configuration
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* ps = AddTaskPhysicsSelection(kFALSE, kTRUE);
  //Signature: Bool_t mCAnalysisFlag, Bool_t applyPileupCuts

  //MultSelection
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  AliMultSelectionTask* ms = AddTaskMultSelection(isMultCalib);
  ms->SetAlternateOADBFullManualBypass("alien:///alice/cern.ch/user/f/fercoles/DataFiles/OADB-LHC12b-MB.root");  
  ms->SetSelectedTriggerClass(AliVEvent::kINT7); // kINT7 is default, this is OK for Run2; in LHC10h you need kMB

  // compile my own code
  if(!isMultCalib){
    gROOT->LoadMacro("AliAnalysisTaskLeading.cxx++g");
    gROOT->LoadMacro("AddTaskLeading.C"); // eventually include the compilation of own task
    AliAnalysisTaskLeading *task = AddTaskLeading();
    if(isMC) task->SetMC();
    else if(isMCdata){
      task->SetMCdata();
    }
  }

  if(!mgr->InitAnalysis()) return;
  mgr->SetDebugLevel(2);
  mgr->PrintStatus();

  if (grid)
    {
      AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
      alienHandler->SetOverwriteMode();

      alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
        
      // Set versions of used packages
      alienHandler->SetAPIVersion("V1.1x");
      //alienHandler->SetROOTVersion("v5-26-00b-6");
      // alienHandler->SetAliROOTVersion("v4-19-21-AN");
      alienHandler->SetAliPhysicsVersion("vAN-20190109-1");   

      // number of files per subjob
      alienHandler->SetSplitMaxInputFileNumber(30);
      //set Executable
      alienHandler->SetExecutable("LeadingESD.sh");
      //specify how many seconds your job may take
      alienHandler->SetTTL(36000);
      //set jdl name
      alienHandler->SetJDLName("LeadingESD.jdl");
      
      alienHandler->SetOutputToRunNo(kTRUE);
      alienHandler->SetKeepLogs(kTRUE);
   
      alienHandler->SetTerminateFiles("event_stat.root");  //to have the output file of the Physics Selection class
      alienHandler->SetInputFormat("xml-single");
      alienHandler->SetPrice(1);      
      // Optionally modify split mode (default 'se')    
      alienHandler->SetSplitMode("se");
       
      // make sure your source files get copied to grid
      alienHandler->SetAdditionalLibs("AliAnalysisTaskLeading.cxx AliAnalysisTaskLeading.h");
      alienHandler->SetAnalysisSource("AliAnalysisTaskLeading.cxx");
      
      //Declare input data to be processed.
      //Method 1: Create automatically XML collections using alien 'find' command.
      
      alienHandler->SetGridDataDir("/alice/data/2012/LHC12b/");
      alienHandler->SetRunPrefix("000");
      // runnumber
      alienHandler->AddRunNumber(178167);
      alienHandler->SetDataPattern("/pass2/*/AliESDs.root");
      alienHandler->SetGridWorkingDir("IICalib12bpass2/Utilizzo/");
             
      // define the output folder
      alienHandler->SetGridOutputDir("OutputDir");
      alienHandler->SetDefaultOutputs();
      
      //number of times to merge (if a lot of data need a higher number)
      alienHandler->SetMaxMergeStages(1);
      // we can specify that we want to, later on, use Grid to also merge
      // our output. to enable this, we will set 'SetMergeViaJDL' to kTRUE
      alienHandler->SetMergeViaJDL(kFALSE);
      //When all your merging jobs have finished running,
      //there is one step to be taken still, which is downloading the output
      //of the merging jobs. This is done by changing SetMergeViaJDL(kFALSE)
      //and running one last time
       
      // Connect plug-in to the analysis manager
      mgr->SetGridHandler(alienHandler);
      
      if(gridTest) {
	// speficy on how many files you want to run
	alienHandler->SetNtestFiles(5);
	// and launch the analysis
	// Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
	alienHandler->SetRunMode("test");
	
	mgr->StartAnalysis("grid");
	
      } else
	{
	  //full grid      
	  alienHandler->SetRunMode("terminate");
	  mgr->StartAnalysis("grid");
	}
      
    }

}

