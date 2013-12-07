#include <iostream>
#include <fstream>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TH1I.h>
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

  dataTree->SetBranchAddress("rawEvent_new", &rawEvent);

  char histName[100];
  TH1D* hist[8][72];
  for(int i = 0; i < 8; ++i)
    {
      for(int j = 0; j < 72; ++j)
	{
	  sprintf(histName, "%s_%i", p_geomSvc->getDetectorName(i+41).c_str(), j);
	  hist[i][j] = new TH1D(histName, histName, 70, 0., 1400.);
	}
    }

  for(Int_t i = 0; i < dataTree->GetEntries(); i++)
    {
      dataTree->GetEntry(i);
      rawEvent->reIndex("a");
      for(int j = 0; j < rawEvent->getNHitsAll(); ++j)
	{
	  Hit h = rawEvent->getHit(j);
	  if(h.detectorID <= 40) continue;

	  hist[h.detectorID-41][h.elementID-1]->Fill(h.tdcTime);
	}

      rawEvent->clear();
    }
  
  TCanvas* c1[8];
  for(int i = 0; i < 8; ++i)
    {
      sprintf(histName, "plane_%d", i);
      c1[i] = new TCanvas(histName, histName);

      c1[i]->Divide(8, 9);
      for(int j = 0; j < 72; j++)
	{
	  c1[i]->cd(j+1);
	  hist[i][j]->Draw();
	}

      sprintf(histName, "plane_%d.pdf", i);
      c1[i]->SaveAs(histName);
    }
  
  return 1;
}
