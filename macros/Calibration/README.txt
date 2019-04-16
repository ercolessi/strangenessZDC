Here we calibrate and create a OADB file containing the information about multiplicity estimators

1) Run the macro CalibMode.C as:                                                                   
  testmode) gridtest=KTRUE 
  if you want to test it in local before grid  
                                                     
  gridmode) gridtest=KFALSE                                                      
  in aliroot                                     
  >> .L CalibMode.C                                  
  >> CalibMode(Period)                                  
  example for LHC12b  Period=122 (year (12) + b second letter of the alphabet)

  This will create a folder in you aliensh AnalysisLeading2019/Calibration/Calib**periodfoldername** which contains folders for each run number (/000*runnumber*/) with calibration information on that period
 

2) Run the macro CalibMerge.C 
  testmode) gridtest=KTRUE 
  if you want to test it in local before grid      
                                                 
  gridmode) gridtest=KFALSE      
  in aliroot              
  >> .L CalibMerge.C   
  >> CalibMerge(Period)  example for LHC12b  Period=122 (year (12) + b second letter of the alphabet)

   This will merge all root files from the subjobs from phase 1 and will create a unique file called AnalysisResults.root for each period


3) Run the macro DoMerge.C (on ROOT6) and merge in one file all root files from all the run numbers for one period and download it in local, don't forget to set the right path of the period in the macro

4) Run CalibratePeriodPP.cxx (on ROOT6) to create the OADB file you are looking for, don't forget to set the name of the period in the macro
