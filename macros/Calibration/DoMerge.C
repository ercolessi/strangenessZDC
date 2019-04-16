#include <TGrid.h>

void DoMerge(const char* period="*bpass2*",const char* output="AnalysisResultsMerged12b.root", const char *path="/alice/cern.ch/user/f/fercoles/AnalysisLeading2019/Calibration/"){
  system(Form("alien_ls -b %s/%s/OutputDir/*/AnalysisResults.root |awk \'{print \"alien://\"$2}\' >listtobemerged",path,period));
 
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

