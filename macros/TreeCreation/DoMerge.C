#include <TGrid.h>

void DoMerge(const char* period="12fpass2",const char* output="LeadingMerged12f.root", const char *path="/alice/cern.ch/user/f/fercoles/AnalysisLeading2019/TreeGeneration/"){
  system(Form("alien_ls -b %s/%s/OutputDir/*/AnalysisLeading.root |awk \'{print \"alien://\"$2}\' >listtobemerged",path,period));
 
  TGrid::Connect("alien://");
  FILE *f = fopen("listtobemerged","r");

  TFileMerger m(kFALSE);
  m.OutputFile(output);

  Int_t i=0;
  char nome[100];
  while (fscanf(f,"%s",nome)==1) {
    m.AddFile(nome);
    i++;
  }
  if (i)
    m.Merge();
}

