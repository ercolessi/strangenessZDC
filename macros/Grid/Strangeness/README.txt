Here we create Tree with information from ZDC and Strangeness info.

1) Run the macro StrangeMode.C as:                                                                   
  testmode) gridtest=KTRUE 
  if you want to test it in local before grid  
                                                     
  gridmode) gridtest=KFALSE                                                      
  in aliroot                                                                     
  >> .L StrangeMode.C                                                                
  >> StrangeMode(Period)                                                             
  example for LHC12b  Period=122 (year (12) + b second letter of the alphabet)

  This will create a folder in you aliensh AnalysisLeading2019/Strangeness/**periodfoldername** which contains folders for each run number (/000*runnumber*/) with ZDC and multiplicity information on that period stored in a root file
 

2) Run the macro StrangeMode.C changing the flag "full" to "terminate" in this line (at the very end):  alienHandler->SetRunMode("full"); 

  testmode) gridtest=KTRUE 
  if you want to test it in local before grid      
                                                 
  gridmode) gridtest=KFALSE      
  in aliroot              
  >> .L StrangeMode.C   
  >> StrangeMode(Period)  example for LHC12b  Period=122 (year (12) + b second letter of the alphabet)

   This will merge all root files from the subjobs from phase 1 and will create a unique file called AnalysisLeading.root for each period


3) Run the macro DoMerge.C and merge in one file all root files from all the run numbers for one period and download it in local

 
