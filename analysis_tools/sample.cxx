#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <algorithm>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TH1I.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TCanvas.h>

#include "GeomSvc.h"
#include "SRawEvent.h"

using namespace std;

int main(int argc, char *argv[])
{
  GeomSvc* p_geomSvc = GeomSvc::instance();
  p_geomSvc->init(GEOMETRY_VERSION);

  SRawEvent* rawEvent = new SRawEvent();

  TFile* dataFile = new TFile(argv[1], "READ");
  TTree* dataTree = (TTree *)dataFile->Get("save");

  dataTree->SetBranchAddress("rawEvent", &rawEvent);


  for(Int_t i = 0; i < dataTree->GetEntries(); i++)
    {
      dataTree->GetEntry(i);

      rawEvent->clear();
    }
  
  return 1;
}
