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

  //Old data
  SRawEvent* event_old = new SRawEvent();

  TFile* dataFile = new TFile(argv[1], "READ");
  TTree* dataTree = (TTree *)dataFile->Get("save");

  dataTree->SetBranchAddress("rawEvent", &event_old);

  Int_t topHits, bottomHits;
  Int_t H1T, H2T, H3T, H4T;
  Int_t H1B, H2B, H3B, H4B;

  SRawEvent* event = new SRawEvent();

  TFile *saveFile = new TFile(argv[2], "recreate");
  TTree *saveTree = new TTree("save", "save");
  
  saveTree->Branch("topHits", &topHits, "topHits/I");
  saveTree->Branch("bottomHits", &bottomHits, "bottomHits/I");
  saveTree->Branch("H1T", &H1T, "H1T/I");
  saveTree->Branch("H2T", &H2T, "H2T/I");
  saveTree->Branch("H3T", &H3T, "H3T/I");
  saveTree->Branch("H4T", &H4T, "H4T/I");
  saveTree->Branch("H1B", &H1B, "H1B/I");
  saveTree->Branch("H2B", &H2B, "H2B/I");
  saveTree->Branch("H3B", &H3B, "H3B/I");
  saveTree->Branch("H4B", &H4B, "H4B/I");
  
  saveTree->Branch("rawEvent", &event, 256000, 99);
 
  double h_center[16] = {880., 880., 890., 890., 900., 900., 890., 890., 970., 970., 975., 975., 970., 970., 965., 965.};
  for(int i = 0; i < dataTree->GetEntries(); i++)
    {
      dataTree->GetEntry(i);

      H1B = event_old->getNHitsInDetector(25) > 0 ? 1 : 0;
      H2B = event_old->getNHitsInDetector(31) > 0 ? 1 : 0;
      H3B = event_old->getNHitsInDetector(33) > 0 ? 1 : 0;
      H4B = event_old->getNHitsInDetector(39) > 0 ? 1 : 0;
      H1T = event_old->getNHitsInDetector(26) > 0 ? 1 : 0;
      H2T = event_old->getNHitsInDetector(32) > 0 ? 1 : 0;
      H3T = event_old->getNHitsInDetector(34) > 0 ? 1 : 0;
      H4T = event_old->getNHitsInDetector(40) > 0 ? 1 : 0;
      topHits = H1T + H2T + H3T + H4T;
      bottomHits = H1B + H2B + H3B + H4B;

      event->setEventInfo(event_old->getRunID(), event_old->getSpillID(), event_old->getEventID());
      std::vector<Hit> hits_old = event_old->getAllHits();
    
      for(unsigned int j = 0; j < hits_old.size(); j++)
	{
	  if(hits_old[j].detectorID >= 25 && hits_old[j].detectorID <= 40)
	    {
	      if(fabs(hits_old[j].tdcTime - h_center[hits_old[j].detectorID - 25]) > 15.) continue;
	    }
	  
	  Hit h = hits_old[j];
	  event->insertHit(h);
	}

      event->reIndex("a");
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
