
////////////////////////////////////////////////////////////////////////////////////
//                                                                                //
// To Run:                                                                        //
// testmode) gridtest=KTRUE                                                       //
// gridmode) gridtest=KFALSE                                                      //
// 1)run with runmode=full and setMergeViaJDL(kTRUE) to run on grid               //
// 2)run with runmode=terminate and SetMergeViaJDL(kTRUE) to merge                //
// 3)run with runmode=terminate and SetMergeViaJDL(kFALSE) to download in local   //
//    the merged root file                                                        //
//                                                                                //
// in aliroot                                                                     //
// >> .L runMode.C                                                                //
// >> runMode(Period)                                                             //
//                         example for LHC12b  Period=122                         //
//                               (year (12) + b second letter of the alphabet)    // 
//                                                                                //
//                                                                                //              
////////////////////////////////////////////////////////////////////////////////////

#if !defined (__CINT__) || defined (__CLING__)
#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisTaskLeading.h"
#include "AliMultSelection.h"
#include "AliMultSelectionTask.h"
#include "AddTaskPhysicsSelection.h"
#include "AliAnalysisTask.h"
#include "AliPhysicsSelectionTask.h"
#include "AliAnalysisTaskPIDResponse.h"
#include "AliAnalysisTaskWeakDecayVertexer.h"
#include "AliAnalysisTaskStrangeness.h"
#endif

void StrangeMode(Int_t Period)
  
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

  //Forward Multiplicity
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/AddTaskForwardMult.C");
  AliAnalysisTask* fw = AddTaskForwardMult(0.);

  //PhysicsSelection Configuration
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* ps = AddTaskPhysicsSelection(kFALSE, kTRUE);
  //Signature: Bool_t mCAnalysisFlag, Bool_t applyPileupCuts
  
  //MultSelection
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  AliMultSelectionTask* ms = AddTaskMultSelection(isMultCalib);
  if(Period==122) ms->SetAlternateOADBFullManualBypass("alien:///alice/cern.ch/user/f/fercoles/DataFiles/OADB-LHC12b_kMB-MB.root");
  if(Period==123) ms->SetAlternateOADBFullManualBypass("alien:///alice/cern.ch/user/f/fercoles/DataFiles/OADB-LHC12c-MB.root");
  if(Period==126) ms->SetAlternateOADBFullManualBypass("alien:///alice/cern.ch/user/f/fercoles/DataFiles/OADB-LHC12f-MB.root");
  if(Period==189) ms->SetAlternateOADBFullManualBypass("alien:///alice/cern.ch/user/f/fercoles/DataFiles/OADB-LHC18i-MB.root");
  if(Period==129) ms->SetAlternateOADBFullManualBypass("alien:///alice/cern.ch/user/f/fercoles/DataFiles/OADB-LHC12i-MB.root");
  if(Period==1710) ms->SetAlternateOADBFullManualBypass("alien:///alice/cern.ch/user/f/fercoles/DataFiles/OADB-LHC17j-MB.root");
  if(Period==168) ms->SetAlternateOADBFullManualBypass("alien:///alice/cern.ch/user/f/fercoles/DataFiles/OADB-LHC16h-MB.root");
  if(Period==158) ms->SetAlternateOADBFullManualBypass("alien:///alice/cern.ch/user/f/fercoles/DataFiles/OADB-LHC15h-MB.root");
  if(Period==156) ms->SetAlternateOADBFullManualBypass("alien:///alice/cern.ch/user/f/fercoles/DataFiles/OADB-LHC15f-MB.root");
  ms->SetSelectedTriggerClass(AliVEvent::kINT7 | AliVEvent::kMB | AliVEvent::kINT5); // kINT7 is default, this is OK for Run2; in LHC10h you need kMB

  //PID Response
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskPIDResponse* pid = AddTaskPIDResponse(0.);
  
  //Weak Decay Vertexer
  gROOT->LoadMacro("AliAnalysisTaskWeakDecayVertexer.cxx++g");
  gROOT->LoadMacro("AddTaskWeakDecayVertexer.C");
  AliAnalysisTaskWeakDecayVertexer * wv = AddTaskWeakDecayVertexer();
  wv-> SetPreselectDedx(kTRUE);
  wv-> SetPreselectDedxLambda(kTRUE);
  wv-> SetExtraCleanup(kTRUE);
  
  //MyTask
  if(!isMultCalib){
    gROOT->LoadMacro("AliAnalysisTaskLeading.cxx++g");
    gROOT->LoadMacro("AddTaskLeading.C"); // eventually include the compilation of own task
    AliAnalysisTaskLeading *task = AddTaskLeading();
    TString n("alien:///alice/cern.ch/user/f/fercoles/DataFiles/EZDC18i.root");
    task->SetZDCFileName(n);
    if(isMC) task->SetMC();
    else if(isMCdata){
      task->SetMCdata();
    }
  }
  
  //Strangness Task
  gROOT->LoadMacro("AliAnalysisTaskStrangeness.cxx++g");
  gROOT->LoadMacro("AddTaskStrangenessVsMultiplicityRun2.C");
  AliAnalysisTaskStrangeness* st = AddTaskStrangenessVsMultiplicityRun2(1,1,1,"ABC");
  st->SetDownScaleV0(kTRUE, 1.0);
  st->SetDownScaleCascade(kTRUE, 1.0);
  st->SetUseExtraEvSels(kFALSE);
  st->SetRunVertexers(kFALSE);
  st->SetUseLightVertexers(kFALSE);
  
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
      alienHandler->SetAnalysisMacro("AnalysisLeading.C"); 

      // number of files per subjob
      alienHandler->SetSplitMaxInputFileNumber(20);
      //set Executable
      alienHandler->SetExecutable("runLeadingESD.sh");
      //specify how many seconds your job may take
      alienHandler->SetTTL(36000);
      //set jdl name
      alienHandler->SetJDLName("runLeadingESD.jdl");
      
      alienHandler->SetOutputToRunNo(kTRUE);
      alienHandler->SetKeepLogs(kTRUE);
   
      alienHandler->SetTerminateFiles("event_stat.root");  //to have the output file of the Physics Selection class
      alienHandler->SetInputFormat("xml-single");
      alienHandler->SetPrice(1);      
      // Optionally modify split mode (default 'se')    
      alienHandler->SetSplitMode("se");
       
      // make sure your source files get copied to grid
      alienHandler->SetAdditionalLibs("AliAnalysisTaskLeading.cxx AliAnalysisTaskLeading.h AliAnalysisTaskStrangeness.cxx AliAnalysisTaskStrangeness.h AliAnalysisTaskWeakDecayVertexer.cxx AliAnalysisTaskWeakDecayVertexer.h");
      alienHandler->SetAnalysisSource("AliAnalysisTaskLeading.cxx AliAnalysisTaskStrangeness.cxx AliAnalysisTaskWeakDecayVertexer.cxx");
      
      //Declare input data to be processed.
      //Method 1: Create automatically XML collections using alien 'find' command.

       //LHC12b-----------------------------------------------
      if (Period==122){
	// select the input data
	alienHandler->SetGridDataDir("/alice/data/2012/LHC12b/");
	alienHandler->SetRunPrefix("000");
	alienHandler->SetDataPattern("/pass2/*/AliESDs.root");
	// runnumber
	Int_t runList[26] = {/*177798, 177799, 177802, 177804, 177805, 177810, 177858, 177860, 177861, 177864, 177866, 177869, 177938, 177942, */178018, 178024, 178025, 178026, 178028, 178029, 178030, 178031, 178052, 178053 /*,178163, 178167*/};
	for (Int_t i= 5;i <10; i++) alienHandler->AddRunNumber(runList[i]);
	
	alienHandler->SetGridWorkingDir("AnalysisLeading2019/Strangeness/12bV0/");
      }

      //LHC12c----------------------------------------------
       if (Period==123){
	// select the input data
	alienHandler->SetGridDataDir("/alice/data/2012/LHC12c/");
	alienHandler->SetRunPrefix("000");
	alienHandler->SetDataPattern("/pass2/*/AliESDs.root");
	// runnumber
	Int_t runList[95] = {179569, 179571, 179584, 179585, 179591, 179618, 179621, 179638, 179639, 179796, 179802, 179803, 179806, 179837, 179858, 179859, 179916, 179917, 179918, 179919, 179920, 180000, 180042, 180044, 180127, 180129, 180130, 180131, 180132, 180133, 180195, 180199, 180200, 180201, 180230, 180500, 180501, 180507, 180510, 180512, 180515, 180517, 180561, 180562, 180564, 180566, 180567, 180569, 180716, 180717, 180719, 180720, 181617, 181618, 181619, 181620, 181652, 181694, 181698, 181701, 181703, 182017, 182018, 182022, 182023, 182110, 182111, 182207, 182289, 182295, 182297, 182299, 182300, 182302, 182322, 182323, 182324, 182325, 182509, 182513, 182624, 182635, 182684, 182686, 182687, 182691, 182692, 182724, 182725, 182728, 182729, 182730, 182740, 182741, 182744};
	for (Int_t i= 80;i <95; i++) alienHandler->AddRunNumber(runList[i]);
	
	alienHandler->SetGridWorkingDir("AnalysisLeading2019/Strangeness/12cpass2/");
      }

        //LHC12f--------------------------------------------
       if (Period==126){
	// select the input data
	alienHandler->SetGridDataDir("/alice/data/2012/LHC12f/");
	alienHandler->SetRunPrefix("000");
	alienHandler->SetDataPattern("/pass2/*/AliESDs.root");
	// runnumber
	Int_t runList[10] = {187749, 187785, 188025, 188027, 188028, 188029, 188093, 188101, 188108, 188123};
	for (Int_t i= 0;i <10; i++) alienHandler->AddRunNumber(runList[i]);
	
	alienHandler->SetGridWorkingDir("AnalysisLeading2019/Strangeness/12fV0/");
      }

       
        //LHC12i--------------------------------------------
       if (Period==129){
	// select the input data
	alienHandler->SetGridDataDir("/alice/data/2012/LHC12i/");
	alienHandler->SetRunPrefix("000");
	alienHandler->SetDataPattern("/pass2/*/AliESDs.root");
	// runnumber
	Int_t runList[30] = {193007, 193008, 193010, 193011, 193014, 193047, 193049, 193051, 193092, 193093, 193094, 193097, 193148, 193150, 193151, 193152, 193155, 193156, 193184, 193187, 193188, 193189, 193192, 193194, 193750, 193751, 193752, 193758, 193759, 193766};
	for (Int_t i= 1;i <30; i++) alienHandler->AddRunNumber(runList[i]);
	
	alienHandler->SetGridWorkingDir("AnalysisLeading2019/Strangeness/12ipass2/");
      }

       
        //LHC17j---------------------------------------------
       if (Period==1710){
	// select the input data
	alienHandler->SetGridDataDir("/alice/data/2017/LHC17j/");
	alienHandler->SetRunPrefix("000");
	alienHandler->SetDataPattern("/pass1/*/AliESDs.root");
	// runnumber
	Int_t runList[10] = {274593,  274594, 274595, 274596, 274601, 274653, 274657, 274667, 274669, 274671};
	for (Int_t i= 0;i<1; i++) alienHandler->AddRunNumber(runList[i]);
	
	alienHandler->SetGridWorkingDir("AnalysisLeading2019/Strangeness/17jpass1_testStrangness/");
      }

           //LHC16h---------------------------------------------
       if (Period==168){
	// select the input data
	alienHandler->SetGridDataDir("/alice/data/2016/LHC16h/");
	alienHandler->SetRunPrefix("000");
	alienHandler->SetDataPattern("/pass1/*/AliESDs.root");
	// runnumber
	Int_t runList[11] = {254378, 254381, 254394, 254395, 254396, 254476, 254479, 255008, 255009, 255010, 255011};
	for (Int_t i= 0;i <11; i++) alienHandler->AddRunNumber(runList[i]);
	
	alienHandler->SetGridWorkingDir("AnalysisLeading2019/Strangeness/16hpass1/");
      }
          
           //LHC15h---------------------------------------------
       if (Period==158){
	// select the input data
	alienHandler->SetGridDataDir("/alice/data/2015/LHC15h/");
	alienHandler->SetRunPrefix("000");
	alienHandler->SetDataPattern("/pass2/*/AliESDs.root");
	// runnumber
	Int_t runList[6] = {233912, 233969, 233971, 233972, 233975, 233976};
	for (Int_t i= 0;i <6; i++) alienHandler->AddRunNumber(runList[i]);
	
	alienHandler->SetGridWorkingDir("AnalysisLeading2019/Strangeness/15hpass2/");
      }
          
           //LHC15f---------------------------------------------
       if (Period==156){
	// select the input data
	alienHandler->SetGridDataDir("/alice/data/2015/LHC15f/");
	alienHandler->SetRunPrefix("000");
	alienHandler->SetDataPattern("/pass2/*/AliESDs.root");
	// runnumber
	Int_t runList[58] = {225016, 225026, 225031, 225035, 225037, 225041, 225043, 225050, 225051, 225052, 225105, 225106, 225305, 225307, 225309, 225310, 225313, 225314, 225315, 225322, 225576, 225578, 225579, 225580, 225582, 225586, 225587, 225589, 225609, 225611, 225705, 225707, 225708, 225709, 225710, 225716, 225717, 225719, 226062, 226085, 226175, 226176, 226177, 226208, 226210, 226212, 226217, 226220, 226225, 226444, 226445, 226452, 226466, 226468, 226472, 226483, 226495, 226500};
	for (Int_t i= 0;i <58; i++) alienHandler->AddRunNumber(runList[i]);
	
	alienHandler->SetGridWorkingDir("AnalysisLeading2019/Strangeness/15fpass2/");
      }
       
       //LHC18i---------------------------------------------
       if (Period==189){
	// select the input data
	alienHandler->SetGridDataDir("/alice/data/2018/LHC18i/");
	alienHandler->SetRunPrefix("000");
	alienHandler->SetDataPattern("/pass1/*/AliESDs.root");
	// runnumber
	Int_t runList[10] = {288861, 288862, 288863, 288864, 288868, 288897, 288902, 288903, 288908, 288909};
	for (Int_t i= 2;i <6; i++) alienHandler->AddRunNumber(runList[i]);
	
	alienHandler->SetGridWorkingDir("AnalysisLeading2019/Strangeness/18ipass1/");
      }
                   
      
      // define the output folder
      alienHandler->SetGridOutputDir("OutputDir");
      alienHandler->SetDefaultOutputs();
      
      //number of times to merge (if a lot of data need a higher number)
      alienHandler->SetMaxMergeStages(5);
      //alienHandler->SetSplitMaxInputFileNumber(150);
      // we can specify that we want to, later on, use Grid to also merge
      // our output. to enable this, we will set 'SetMergeViaJDL' to kTRUE
      alienHandler->SetMergeViaJDL(kTRUE);
      //When all your merging jobs have finished running,
      //there is one step to be taken still, which is downloading the output
      //of the merging jobs. This is done by changing SetMergeViaJDL(kFALSE)
      //and running one last time
       
      // Connect plug-in to the analysis manager
      mgr->SetGridHandler(alienHandler);
      
      if(gridTest) {
	// speficy on how many files you want to run
	alienHandler->SetNtestFiles(1);
	// and launch the analysis
	// Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
	alienHandler->SetRunMode("test");
	mgr->StartAnalysis("grid");
	
      } else
	{
	  //full grid      
	  alienHandler->SetRunMode("full");
	  mgr->StartAnalysis("grid");
	}
      
    }

}

