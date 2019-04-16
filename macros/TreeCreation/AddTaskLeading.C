AliAnalysisTask *AddTaskLeading(){
  char *fileout = "AnalysisLeading.root";

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskLeading", "No analysis manager found.");
    return 0;
  }

  if (!mgr->GetInputEventHandler()->InheritsFrom(AliESDInputHandler::Class())) {
    Error("AddTask_tender_Tender","The analysis tender only works with ESD input!");
    return 0;
  }

  AliAnalysisTaskLeading *task = new AliAnalysisTaskLeading();

  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //Attach input to my tasks
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());

  // Attach output to my tasks
  Int_t nList = task->GetNList();
  printf("number of TList = %i\n",nList);
  for(Int_t i=0;i < nList;i++){
    char nameObj[100];
    sprintf(nameObj,"%s",task->GetNameList(i));
    AliAnalysisDataContainer *cOutputL= mgr->CreateContainer(nameObj,TList::Class(), AliAnalysisManager::kOutputContainer,fileout);
    mgr->ConnectOutput(task, i, cOutputL);
  }

  Int_t nTh1d = task->GetNTh1d();
  printf("number of TH1D = %i\n",nTh1d);
  for(Int_t i=0;i < nTh1d;i++){
    char nome[100];
    sprintf(nome,"cOutputH1_%i",i+1);
    AliAnalysisDataContainer *cOutputH1= mgr->CreateContainer(nome,TH1D::Class(), AliAnalysisManager::kOutputContainer,fileout);
    mgr->ConnectOutput(task, i+nList, cOutputH1);
  }

  Int_t nTh2d = task->GetNTh2d();
  printf("number of TH2D = %i\n",nTh2d);
  for(Int_t i=0;i < nTh2d;i++){
    char nome[100];
    sprintf(nome,"cOutputH2",i+1);
    AliAnalysisDataContainer *cOutputH2= mgr->CreateContainer(nome,TH2D::Class(), AliAnalysisManager::kOutputContainer,fileout);
    mgr->ConnectOutput(task, i+nTh1d+nList, cOutputH2);
  }

  Int_t nTree = task->GetNTree();
  for(Int_t i=0;i < nTree;i++){
    char nome[100];
    sprintf(nome,"cOutputT_%i",i+1);
    AliAnalysisDataContainer *cOutputT= mgr->CreateContainer(nome,TTree::Class(), AliAnalysisManager::kOutputContainer,fileout);
    mgr->ConnectOutput(task, i+nTh1d+nTh2d+nList, cOutputT);
  }

  return task;
}
