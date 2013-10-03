#include <iostream>
#include <cmath>
#include <algorithm>
#include <string>
#include <time.h>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TRandom.h>
#include <TMatrixD.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>

#include "GeomSvc.h"
#include "SRawEvent.h"
#include "KalmanUtil.h"
#include "KalmanTrack.h"
#include "KalmanFilter.h"
#include "KalmanFinder.h"
#include "KalmanFitter.h"
#include "VertexFit.h"
#include "SRecEvent.h"

using namespace std;

int main(int argc, char *argv[])
{
  //Initialize geometry service
  Log("Initializing geometry service ... ");
  GeomSvc* geometrySvc = GeomSvc::instance();
  geometrySvc->init(GEOMETRY_VERSION);

  //Retrieve the raw event
  Log("Retrieving the event stored in ROOT file ... ");

  SRawEvent *rawEvent = new SRawEvent();
  int nSeeds, nSeedsX, nSeedsY;
  double ax[500], ay[500], bx[500], by[500];
  int xIndex[500], yIndex[500];

  TFile *dataFile = new TFile(argv[1], "READ");
  TTree *dataTree = (TTree *)dataFile->Get("save");

  dataTree->SetBranchAddress("rawEvent", &rawEvent);
  dataTree->SetBranchAddress("nSeeds", &nSeeds);
  dataTree->SetBranchAddress("nSeedsX", &nSeedsX);
  dataTree->SetBranchAddress("nSeedsY", &nSeedsY);
  dataTree->SetBranchAddress("ax", ax);
  dataTree->SetBranchAddress("bx", bx);
  dataTree->SetBranchAddress("ay", ay);
  dataTree->SetBranchAddress("by", by);
  dataTree->SetBranchAddress("xIndex", xIndex);
  dataTree->SetBranchAddress("yIndex", yIndex);

  double time;

  int nTracks;
  int charge[30];
  int nRecHits[30];
  double chisq[30];
  int nHodoHits[30][3];
  double p_vertex[30], px_vertex[30], py_vertex[30], pz_vertex[30];
  double p_final[30], px_final[30], py_final[30], pz_final[30];
  //double p_seed[30], px_seed[30], py_seed[30], pz_seed[30];
  double z_vertex[30];
  double chisq_vertex[30];
  
  SRecEvent *recEvent = new SRecEvent();

  TFile *saveFile = new TFile(argv[2], "recreate");
  TTree *saveTree = dataTree->CloneTree(0);

  saveTree->Branch("recEvent", &recEvent, 256000, 99);
  saveTree->Branch("time", &time, "time/D");
 
  saveTree->Branch("nTracks", &nTracks, "nTracks/I");
  saveTree->Branch("charge", charge, "charge[nTracks]/I");
  saveTree->Branch("nRecHits", nRecHits, "nRecHits[nTracks]/I");
  saveTree->Branch("chisq", chisq, "chisq[nTracks]/D");
  saveTree->Branch("nHodoHits", nHodoHits, "nHodoHits[nTracks][3]/I");
  saveTree->Branch("p_vertex", p_vertex, "p_vertex[nTracks]/D");
  saveTree->Branch("px_vertex", px_vertex, "px_vertex[nTracks]/D");
  saveTree->Branch("py_vertex", py_vertex, "py_vertex[nTracks]/D");
  saveTree->Branch("pz_vertex", pz_vertex, "pz_vertex[nTracks]/D");
  //saveTree->Branch("p_seed", p_seed, "p_seed[nTracks]/D");
  //saveTree->Branch("px_seed", px_seed, "px_seed[nTracks]/D");
  //saveTree->Branch("py_seed", py_seed, "py_seed[nTracks]/D");
  //saveTree->Branch("pz_seed", pz_seed, "pz_seed[nTracks]/D");
  saveTree->Branch("p_final", p_final, "p_final[nTracks]/D");
  saveTree->Branch("px_final", px_final, "px_final[nTracks]/D");
  saveTree->Branch("py_final", py_final, "py_final[nTracks]/D");
  saveTree->Branch("pz_final", pz_final, "pz_final[nTracks]/D");

  saveTree->Branch("z_vertex", z_vertex, "z_vertex[nTracks]/D");
  saveTree->Branch("chisq_vertex", chisq_vertex, "chisq_vertex[nTracks]/D");

  //Initialize track finder
  Log("Initializing the track finder and kalman filter ... ");
  KalmanFilter *filter = new KalmanFilter();
  KalmanFinder *finder = new KalmanFinder();
  KalmanFitter *fitter = new KalmanFitter();
  VertexFit    *vtxfit = new VertexFit();

  int offset = argc > 3 ? atoi(argv[3]) : 0;
  int nEvtMax = argc > 4 ? atoi(argv[4]) + offset : dataTree->GetEntries();
  if(nEvtMax > dataTree->GetEntries()) nEvtMax = dataTree->GetEntries();
  Log("Running from event " << offset << " through to event " << nEvtMax);
  for(int i = offset; i < nEvtMax; i++)
    {
      dataTree->GetEntry(i);
      Log("Processing event " << i << " with eventID = " << rawEvent->getEventID());

      if((rawEvent->getRunID() != 2166) && (nSeeds < 2 || nSeedsX < 2)) continue;
      if(nSeeds > 50) continue;

      rawEvent->reIndex("oah");
      finder->setEvent(rawEvent);

      std::list<Seed> _seeds;
      for(int j = 0; j < nSeeds; j++)
	{
    	  Seed _seed;
	  _seed.axz = ax[xIndex[j]]; 
	  _seed.ayz = ay[yIndex[j]];
	  _seed.bxz = bx[xIndex[j]];
	  _seed.byz = by[yIndex[j]];
	  _seed.decideCharge();
	  _seed.decideMomentum();

	  _seeds.push_back(_seed);
	}
      Seed::reduceSeedList(_seeds);

      ///Timing stack.getPositions(3, x, y, z);rts ...
      clock_t time_single = clock();

      std::list<KalmanTrack> _tracks;
      _tracks.clear();
      for(std::list<Seed>::iterator iter = _seeds.begin(); iter != _seeds.end(); ++iter)
	{
	  Log("Working on this seed ... ");
	  finder->processOneSeed(*iter);
	  KalmanTrack _track = finder->getBestCandidate();

	  if(_track.isValid())
	    {
	      _tracks.push_back(_track);
	    }
	}

      time_single = clock() - time_single;
      time = double(time_single)/CLOCKS_PER_SEC;
      Log("It takes " << time << " seconds for this event.");
     
      std::list<KalmanTrack> _tracks_final;// = _tracks;
      finder->reduceTrackList(_tracks, _tracks_final);

      Log("Found " << _tracks_final.size() << " tracks!");
      if(_tracks_final.empty()) continue;

      //finder->printResults("test", _tracks_final);

      nTracks = 0;
      recEvent->setRawEvent(rawEvent);
      for(std::list<KalmanTrack>::iterator iter = _tracks_final.begin(); iter != _tracks_final.end(); ++iter)
	{
	  iter->print();

	  SRecTrack recTrack = iter->getSRecTrack();

	  charge[nTracks] = iter->getCharge();
	  chisq[nTracks] = iter->getChisq();
	  nRecHits[nTracks] = iter->getNHits();
	  
	  for(int j = 0; j < 3; j++)
	    {
	      nHodoHits[nTracks][j] = finder->getNHodoHits(j+1, *iter);
	    }
          recTrack.setHodoHits(nHodoHits[nTracks]);
	  recEvent->insertTrack(recTrack);

	  z_vertex[nTracks] = vtxfit->findSingleMuonVertex(iter->getNodeList().front());
	  p_vertex[nTracks] = iter->getMomentumVertex(z_vertex[nTracks], px_vertex[nTracks], py_vertex[nTracks], pz_vertex[nTracks]);
	  chisq_vertex[nTracks] = iter->getChisqVertex();

	  p_final[nTracks] = iter->getMomentumUpstream(px_final[nTracks], py_final[nTracks], pz_final[nTracks]);
	  
	  nTracks++;
	}

      recEvent->reIndex();
      saveTree->Fill();
      
      rawEvent->clear();
      recEvent->clear();
    }

  saveFile->cd();
  saveTree->Write();
  saveFile->Close();

  filter->close();
  delete finder;
  delete fitter;
  delete vtxfit;

  return 1;
}
