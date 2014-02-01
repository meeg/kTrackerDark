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
#include "JobOptsSvc.h"
#include "SRawEvent.h"
#include "SRecEvent.h"
#include "FastTracklet.h"
#include "KalmanFastTracking.h"
#include "KalmanFitter.h"
#include "VertexFit.h"

#include "MODE_SWITCH.h"

using namespace std;

int main(int argc, char *argv[])
{

  //Initialize job options
  LogInfo("Initializing job options service");
  JobOptsSvc* jobOptsSvc = JobOptsSvc::instance();
  jobOptsSvc->init();

  //Initialize geometry service
  LogInfo("Initializing geometry service ... ");
  GeomSvc* geometrySvc = GeomSvc::instance();
  geometrySvc->init(GEOMETRY_VERSION);

  //Retrieve the raw event
  LogInfo("Retrieving the event stored in ROOT file ... ");
#ifdef MC_MODE
  SRawMCEvent* rawEvent = new SRawMCEvent();
#else
  SRawEvent* rawEvent = new SRawEvent();
#endif

  TFile *dataFile = new TFile(argv[1], "READ");
  TTree *dataTree = (TTree *)dataFile->Get("save");

  dataTree->SetBranchAddress("rawEvent", &rawEvent);
 
  //Output definition
  int nTracklets;
  TClonesArray* tracklets = new TClonesArray("Tracklet");
  TClonesArray& arr_tracklets = *tracklets;

  //int nTracklets_back;
  //TClonesArray* tracklets_back = new TClonesArray("Tracklet");
  //TClonesArray& arr_tracklets_back = *tracklets_back;

  double time;
  SRecEvent* recEvent = new SRecEvent();

  TFile* saveFile = new TFile(argv[2], "recreate");
  TTree* saveTree = dataTree->CloneTree(0);

  saveTree->Branch("recEvent", &recEvent, 256000, 99);
  saveTree->Branch("time", &time, "time/D");
  saveTree->Branch("nTracklets", &nTracklets, "nTracklets/I");
  saveTree->Branch("tracklets", &tracklets, 256000, 99);
  //saveTree->Branch("nTracklets_back", &nTracklets_back, "nTracklets_back/I");
  //saveTree->Branch("tracklets_back", &tracklets_back, 256000, 99);
  tracklets->BypassStreamer();
  //tracklets_back->BypassStreamer();

  //Initialize track finder
  LogInfo("Initializing the track finder and kalman filter ... ");
#ifdef _ENABLE_KF
  KalmanFilter* filter = new KalmanFilter();
  KalmanFastTracking* fastfinder = new KalmanFastTracking();
  VertexFit* vtxfit = new VertexFit();
#else
  KalmanFastTracking* fastfinder = new KalmanFastTracking(false);
#endif

  int offset = argc > 3 ? atoi(argv[3]) : 0;
  int nEvtMax = argc > 4 ? atoi(argv[4]) + offset : dataTree->GetEntries();
  if(nEvtMax > dataTree->GetEntries()) nEvtMax = dataTree->GetEntries();
  LogInfo("Running from event " << offset << " through to event " << nEvtMax);
  for(int i = offset; i < nEvtMax; i++)
    {
      dataTree->GetEntry(i);
      cout << "\r Processing event " << i << " with eventID = " << rawEvent->getEventID() << ", ";
      cout << (i - offset + 1)*100/(nEvtMax - offset) << "% finished .. ";

      clock_t time_single = clock();

      rawEvent->reIndex("ao");
      if(!fastfinder->setRawEvent(rawEvent)) continue;

      //Fill the TClonesArray
      arr_tracklets.Clear();
      //arr_tracklets_back.Clear();
      std::list<Tracklet>& rec_tracklets = fastfinder->getFinalTracklets();
      //std::list<Tracklet>& rec_tracklets_back = fastfinder->getBackPartials();
      if(rec_tracklets.empty()) continue;

      /*
      nTracklets_back = 0;
      for(std::list<Tracklet>::iterator iter = rec_tracklets_back.begin(); iter != rec_tracklets_back.end(); ++iter)
	{
	  iter->calcChisq();

	  new(arr_tracklets_back[nTracklets_back]) Tracklet(*iter);
	  ++nTracklets_back;
	}
      */

      nTracklets = 0;
      recEvent->setRawEvent(rawEvent);
      for(std::list<Tracklet>::iterator iter = rec_tracklets.begin(); iter != rec_tracklets.end(); ++iter)
	{
	  iter->calcChisq();
	  //iter->print();
	  new(arr_tracklets[nTracklets]) Tracklet(*iter);
	  ++nTracklets;

#ifndef _ENABLE_KF
	  SRecTrack recTrack = iter->getSRecTrack();
	  recTrack.setHodoHits();      //Any track that can be reconstructed already required hodo and prop. tube masking
	  recEvent->insertTrack(recTrack);
#endif
	}

#ifdef _ENABLE_KF
      std::list<KalmanTrack>& rec_tracks = fastfinder->getKalmanTracks();
      for(std::list<KalmanTrack>::iterator iter = rec_tracks.begin(); iter != rec_tracks.end(); ++iter)
	{
	  //iter->print();
	  SRecTrack recTrack = iter->getSRecTrack();
	  recTrack.setHodoHits(); 
	  recTrack.setZVertex(vtxfit->findSingleMuonVertex(recTrack));
          recEvent->insertTrack(recTrack);
	}
#endif

      time_single = clock() - time_single;
      time = double(time_single)/CLOCKS_PER_SEC;
      cout << "it takes " << time << " seconds for this event." << flush;

      recEvent->reIndex();
      saveTree->Fill();
      
      recEvent->clear();
      rawEvent->clear();
    }
  cout << endl;
  cout << "kFastTracking ended successfully." << endl;

  saveFile->cd();
  saveTree->Write();
  saveFile->Close();

  delete fastfinder;
#ifdef _ENABLE_KF
  filter->close();
#endif

  return 1;
}
