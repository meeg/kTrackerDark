#include <iostream>
#include <cmath>
#include <algorithm>
#include <string>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>

#include "GeomSvc.h"
#include "SRawEvent.h"

using namespace std;

int main(int argc, char *argv[])
{
  //Initialize geometry service
  GeomSvc* geometrySvc = GeomSvc::instance();
  geometrySvc->init(GEOMETRY_VERSION);
  geometrySvc->loadCalibration("calibration.txt");

  //Old data
#ifndef MC_MODE
  SRawEvent* rawEvent = new SRawEvent();
#else
  SRawMCEvent* rawEvent = new SRawMCEvent();
#endif

  TFile *dataFile = new TFile(argv[1], "READ");
  TTree *dataTree = (TTree *)dataFile->Get("save");

  dataTree->SetBranchAddress("rawEvent", &rawEvent);
 
  TFile *saveFile = new TFile(argv[2], "recreate");
  TTree *saveTree = dataTree->CloneTree(0);
 
  for(int i = 0; i < dataTree->GetEntries(); i++)
    {
      dataTree->GetEntry(i);

      for(unsigned int j = 0; j < rawEvent->getNHitsAll(); ++j)
	{
	  Hit h = rawEvent->getHit(j);

	  h.pos = geometrySvc->getMeasurement(h.detectorID, h.elementID);
	  if((h.detectorID <= 24 || h.detectorID >= 41) && h.inTime > 0 && geometrySvc->isCalibrationLoaded())
	    {
	      h.driftDistance = geometrySvc->getDriftDistance(h.detectorID, h.tdcTime);
	    }

	  rawEvent->setHit(j, h);
	}

      for(unsigned int j = 0; j < rawEvent->getNTriggerHits(); ++j)
	{
	  Hit h = rawEvent->getTriggerHit(j);
	  h.pos = geometrySvc->getMeasurement(h.detectorID, h.elementID);

	  rawEvent->setTriggerHit(j, h);
	}

      saveTree->Fill();
      rawEvent->clear();
    }

  saveFile->cd();
  saveTree->Write();
  saveFile->Close();

  if(argc > 3)
    {
      char cmd[300];
      sprintf(cmd, "mv %s %s", argv[2], argv[1]);
      cout << cmd << endl;
      system(cmd);
    }

  return 1;
}
