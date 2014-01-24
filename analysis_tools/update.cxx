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
  SRawEvent* event_old = new SRawEvent();
  SRawEvent* event = new SRawEvent();
#else
  SRawMCEvent* event_old = new SRawMCEvent();
  SRawMCEvent* event = new SRawMCEvent();
#endif

  TFile *dataFile = new TFile(argv[1], "READ");
  TTree *dataTree = (TTree *)dataFile->Get("save");

  dataTree->SetBranchAddress("rawEvent", &event_old);
 
  TFile *saveFile = new TFile(argv[2], "recreate");
  TTree *saveTree = new TTree("save", "save");

  saveTree->Branch("rawEvent", &event, 256000, 99);
 
  for(int i = 0; i < dataTree->GetEntries(); i++)
    {
      dataTree->GetEntry(i);

      event->setEventInfo(event_old);
      std::vector<Hit> hits_old = event_old->getAllHits();
      std::vector<Hit> thits_old = event_old->getTriggerHits();

      for(unsigned int j = 0; j < hits_old.size(); j++)
	{
	  Hit h;

	  h.index = hits_old[j].index;
	  h.detectorID = hits_old[j].detectorID;
	  h.elementID = hits_old[j].elementID;
	  h.tdcTime = hits_old[j].tdcTime;
	  h.driftTime = hits_old[j].driftTime;
	  h.inTime = hits_old[j].inTime;
	  h.pos = geometrySvc->getMeasurement(h.detectorID, h.elementID);
	  h.hodoMask = 1;

	  if(h.detectorID <= 24 && h.inTime > 0 && geometrySvc->isCalibrationLoaded())
	    {
	      h.driftDistance = geometrySvc->getDriftDistance(h.detectorID, h.tdcTime);
	    }
	  else
	    {
	      h.driftDistance = hits_old[j].driftDistance;
	    }

	  event->insertHit(h);
	}

      for(unsigned int j = 0; j < thits_old.size(); ++j)
	{
	  event->insertTriggerHit(thits_old[j]);
	}

      saveTree->Fill();
      
      event_old->clear();
      event->clear();
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
