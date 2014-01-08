#include <iostream>
#include <cmath>
#include <algorithm>
#include <string>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>

#include "GeomSvc.h"
#include "SRawEvent.h"

using namespace std;

int main(int argc, char *argv[])
{
  //Initialize geometry service
  GeomSvc* p_geomSvc = GeomSvc::instance();
  p_geomSvc->init(GEOMETRY_VERSION);
  p_geomSvc->loadCalibration("calibration.txt");

  //Old data
  SRawEvent* event_old = new SRawEvent();

  TFile* dataFile = new TFile(argv[1], "READ");
  TTree* dataTree = (TTree *)dataFile->Get("save");

  dataTree->SetBranchAddress("rawEvent", &event_old);

  SRawEvent* event = new SRawEvent();

  TFile *saveFile = new TFile(argv[2], "recreate");
  TTree *saveTree = new TTree("save", "save");
  
  saveTree->Branch("rawEvent", &event, 256000, 99);

  //Look for the intime peak 
  double h_center[16];
  TH1D* hist[16];
  for(int i = 0; i < 16; ++i)
    {
      hist[i] = new TH1D(p_geomSvc->getDetectorName(i+25).c_str(), p_geomSvc->getDetectorName(i+25).c_str(), 200, 800, 1200);
    }

  for(int i = 0; i < dataTree->GetEntries(); ++i)
    {
      dataTree->GetEntry(i);

      std::vector<Hit> hits = event_old->getAllHits();
      for(unsigned int j = 0; j < hits.size(); j++)
	{
	  if(hits[j].detectorID < 25 || hits[j].detectorID > 40) continue;
	  hist[hits[j].detectorID - 25]->Fill(hits[j].tdcTime);
	}
    }

  for(int i = 0; i < 16; ++i)
    {
      h_center[i] = hist[i]->GetBinCenter(hist[i]->GetMaximumBin());
      cout << p_geomSvc->getDetectorName(i+25) << "  " << h_center[i] << endl;
    }

  //Remove the out-of-time hodo hits
  for(int i = 0; i < dataTree->GetEntries(); ++i)
    {
      dataTree->GetEntry(i);

      event->setEventInfo(event_old);
      std::vector<Hit> hits = event_old->getAllHits();
      std::vector<Hit> hits_trigger = event_old->getTriggerHits();

      for(unsigned int j = 0; j < hits.size(); j++)
	{
	  Hit h = hits[j];
	  if(h.detectorID <= 24)
	    {
	      h.driftDistance = p_geomSvc->getDriftDistance(h.detectorID, h.tdcTime);
	      h.inTime = p_geomSvc->isInTime(h.detectorID, h.tdcTime) ? 1 : 0;
	    }
	  if(h.detectorID > 24 && h.detectorID <= 40) h.inTime = fabs(h.tdcTime - h_center[h.detectorID-25]) < 15. ? 1 : 0;
	  if(h.detectorID > 40) h.inTime = h.tdcTime > 450. && h.tdcTime < 1100. ? 1 : 0;
	  
	  event->insertHit(h);
	}

      for(unsigned int j = 0; j < hits_trigger.size(); ++j)
	{
	  Hit h = hits_trigger[j];
	  event->insertTriggerHit(h);
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
