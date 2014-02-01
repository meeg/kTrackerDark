#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>

#include "SRawEvent.h"
#include "GeomSvc.h"
#include "SeedFinder.h"

using namespace std;

int main(int argc, char* argv[])
{
  GeomSvc *geomSvc = GeomSvc::instance();
  geomSvc->init( );

#ifdef MC_MODE
  SRawMCEvent* rawEvent = new SRawMCEvent();
#else
  SRawEvent* rawEvent = new SRawEvent();
#endif

  TFile *dataFile = new TFile(argv[1], "READ");
  TTree *dataTree = (TTree *)dataFile->Get("save");

  dataTree->SetBranchAddress("rawEvent", &rawEvent);

  Int_t runID, spillID, eventID;
  Int_t nSeeds, nSeedsX, nSeedsY;
  Int_t nSeedHitsX[500], nSeedHitsY[500];
  Double_t ax[500], ay[500], bx[500], by[500];
  Double_t xchisq[500], ychisq[500];
  Int_t xIndex[500], yIndex[500];

  TFile *saveFile = new TFile(argv[2], "recreate");
  TTree *saveTree = dataTree->CloneTree(0);

  saveTree->Branch("runID", &runID, "runID/I");
  saveTree->Branch("spillID", &spillID, "spillID/I");
  saveTree->Branch("eventID", &eventID, "eventID/I");

  saveTree->Branch("nSeedsX", &nSeedsX, "nSeedsX/I");
  saveTree->Branch("nSeedHitsX", nSeedHitsX, "nSeedHitsX[nSeedsX]/I");
  saveTree->Branch("ax", ax, "ax[nSeedsX]/D");
  saveTree->Branch("bx", bx, "bx[nSeedsX]/D");
  saveTree->Branch("xchisq", xchisq, "xchisq[nSeedsX]/D");

  saveTree->Branch("nSeedsY", &nSeedsY, "nSeedsY/I");
  saveTree->Branch("nSeedHitsY", nSeedHitsY, "nSeedHitsY[nSeedsY]/I");
  saveTree->Branch("ay", ay, "ay[nSeedsY]/D");
  saveTree->Branch("by", by, "by[nSeedsY]/D");
  saveTree->Branch("ychisq", ychisq, "ychisq[nSeedsY]/D");

  saveTree->Branch("nSeeds", &nSeeds, "nSeeds/I");
  saveTree->Branch("xIndex", xIndex, "xIndex[nSeeds]/I");
  saveTree->Branch("yIndex", yIndex, "yIndex[nSeeds]/I");
 
  Int_t nEventMax = argc > 3 ? atoi(argv[3]) : dataTree->GetEntries();
  SeedFinder *seeder = new SeedFinder();

  Int_t array_masks[] = {33, 34, 35, 36, 37, 38, 39, 40};
  vector<Int_t> masks(array_masks, array_masks + sizeof(array_masks)/sizeof(Int_t));
  
  Int_t array_start_X[] = {43, 44};
  Int_t array_end_X[] = {45, 46};
  Int_t array_all_X[] = {43, 44, 45, 46, 33, 34, 39, 40};
  vector<Int_t> start_X(array_start_X, array_start_X + sizeof(array_start_X)/sizeof(Int_t));
  vector<Int_t> end_X(array_end_X, array_end_X + sizeof(array_end_X)/sizeof(Int_t));
  vector<Int_t> all_X(array_all_X, array_all_X + sizeof(array_all_X)/sizeof(Int_t));
 
  Int_t array_start_Y[] = {41, 42};
  Int_t array_end_Y[] = {47, 48};
  Int_t array_all_Y[] = {41, 42, 47, 48, 35, 36, 37, 38};
  vector<Int_t> start_Y(array_start_Y, array_start_Y + sizeof(array_start_Y)/sizeof(Int_t));
  vector<Int_t> end_Y(array_end_Y, array_end_Y + sizeof(array_end_Y)/sizeof(Int_t));
  vector<Int_t> all_Y(array_all_Y, array_all_Y + sizeof(array_all_Y)/sizeof(Int_t));
 
  seeder->setMaskIDs(masks);
  Int_t nGoodXEvent = 0;
  Int_t nGoodYEvent = 0;
  for(Int_t i = 0; i < nEventMax; i++)
    {
      dataTree->GetEntry(i);

      cout << "Processing event " << i << " with RunID = " << rawEvent->getRunID() << " and eventID = " << rawEvent->getEventID() << ": " << endl; 

      runID = rawEvent->getRunID();
      spillID = rawEvent->getSpillID();
      eventID = rawEvent->getEventID();
      
      rawEvent->reIndex("a");

      seeder->setDetectorIDs(start_X, end_X, all_X);
      seeder->setSeedLimits(0.2, 200.);
      
      nSeedsX = seeder->processOneEvent(rawEvent);
      list<Seed1D> xseeds = seeder->getFinalSeeds();
      ++nGoodXEvent;
      //seeder->print();

      //char buf[300];
      //sprintf(buf, "seed_run%d_event%d_%d_x_z.jpg", event->getRunID(), i, event->getEventID());     
      //seeder->printResults(string(buf)); 

      seeder->setDetectorIDs(start_Y, end_Y, all_Y);
      seeder->setSeedLimits(0.1, 100.);
      
      nSeedsY = seeder->processOneEvent(rawEvent);
      list<Seed1D> yseeds = seeder->getFinalSeeds();
      ++nGoodYEvent;
      //seeder->print();

      //sprintf(buf, "seed_run%d_event%d_%d_y_z.jpg", event->getRunID(), i, event->getEventID());
      //seeder->printResults(string(buf));

      if(nSeedsX < 1 || nSeedsY < 1) continue;

      nSeeds = 0;
      Int_t index_x = -1;
      Int_t index_y = -1;
      for(list<Seed1D>::iterator seedx = xseeds.begin(); seedx != xseeds.end(); ++seedx)
	{
	  index_x++;
	  index_y = -1;
	  for(list<Seed1D>::iterator seedy = yseeds.begin(); seedy != yseeds.end(); ++seedy)
	    {
	      index_y++;
	      cout << "Testing " << index_x << " with " << index_y << endl;
	      if(!seeder->hodoMask(*seedx, *seedy)) continue;
	      xIndex[nSeeds] = index_x;
	      yIndex[nSeeds] = index_y;

	      nSeeds++;
	    }
	}

      if(nSeeds < 1 || nSeeds > 200) continue;    
      
      index_x = 0; 
      for(list<Seed1D>::iterator seedx = xseeds.begin(); seedx != xseeds.end(); ++seedx)
        {
	  nSeedHitsX[index_x] = seedx->nHits;
	  ax[index_x] = seedx->ax;
	  bx[index_x] = seedx->bx;
	  xchisq[index_x] = seedx->chisq;

	  index_x++;
	}	  

      index_y = 0; 
      for(list<Seed1D>::iterator seedy = yseeds.begin(); seedy != yseeds.end(); ++seedy)
        {
	  nSeedHitsY[index_y] = seedy->nHits;
	  ay[index_y] = seedy->ax;
	  by[index_y] = seedy->bx;
	  ychisq[index_y] = seedy->chisq;

	  index_y++;
	}	

      saveTree->Fill();
      rawEvent->clear();
    } 

  cout << "Nevent that passed basic quality cut: " << nGoodXEvent << "  " << nGoodYEvent << endl;
  cout << saveTree->GetEntries() << " dimuons found" << endl;

  saveFile->cd();
  saveTree->Write();
  saveFile->Close();

  return 1;
}
