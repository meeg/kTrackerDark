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
#include "SRecEvent.h"
#include "FastTracklet.h"
#include "KalmanFastTracking.h"
#include "KalmanFitter.h"
#include "VertexFit.h"
#include "TriggerAnalyzer.h"
#include "JobOptsSvc.h"
#include "MODE_SWITCH.h"

using namespace std;

int main(int argc, char *argv[])
{
  if(argc != 2)
    {
      cout << "Usage: " << argv[0] << "  <options file>" << endl;
      exit(0);
    }

  //Initialize job options
  JobOptsSvc* jobOptsSvc = JobOptsSvc::instance();
  jobOptsSvc->init(argv[1]);

  //Retrieve the raw event
  SRawEvent* rawEvent = jobOptsSvc->m_mcMode ? (new SRawMCEvent()) : (new SRawEvent());

  TFile* dataFile = new TFile(jobOptsSvc->m_inputFile.c_str(), "READ");
  TTree* dataTree = (TTree *)dataFile->Get("save");

  dataTree->SetBranchAddress("rawEvent", &rawEvent);
 
  //Output definition
  int nTracklets;
  TClonesArray* tracklets = new TClonesArray("Tracklet");
  TClonesArray& arr_tracklets = *tracklets;

  double time;
  SRecEvent* recEvent = new SRecEvent();

  TFile* saveFile = new TFile(jobOptsSvc->m_outputFile.c_str(), "recreate");
  TTree* saveTree = dataTree->CloneTree(0);

  saveTree->Branch("recEvent", &recEvent, 256000, 99);
  saveTree->Branch("time", &time, "time/D");
  saveTree->Branch("nTracklets", &nTracklets, "nTracklets/I");
  saveTree->Branch("tracklets", &tracklets, 256000, 99);
  tracklets->BypassStreamer();

  //Initialize track finder
  LogInfo("Initializing the track finder and kalman filter ... ");
#ifdef _ENABLE_KF
  KalmanFilter* filter = new KalmanFilter();
  KalmanFastTracking* fastfinder = new KalmanFastTracking();
#else
  KalmanFastTracking* fastfinder = new KalmanFastTracking(false);
#endif

  TriggerAnalyzer* triggerAna = new TriggerAnalyzer();
  if(jobOptsSvc->m_enableTriggerMask)
    {
      triggerAna->init();
      triggerAna->buildTriggerTree();
    }

  const int offset = jobOptsSvc->m_firstEvent;
  int nEvtMax = jobOptsSvc->m_nEvents > 0 ? jobOptsSvc->m_nEvents + offset : dataTree->GetEntries();
  if(nEvtMax > dataTree->GetEntries()) nEvtMax = dataTree->GetEntries();
  LogInfo("Running from event " << offset << " through to event " << nEvtMax);
  for(int i = offset; i < nEvtMax; ++i)
    {
      dataTree->GetEntry(i);
      cout << "\r Processing event " << i << " with eventID = " << rawEvent->getEventID() << ", ";
      cout << (i - offset + 1)*100/(nEvtMax - offset) << "% finished .. ";

      clock_t time_single = clock();

      if(jobOptsSvc->m_enableTriggerMask)
	{
	  triggerAna->trimEvent(rawEvent);
	  rawEvent->reIndex("aoct");
	}
      else
	{
	  rawEvent->reIndex("aoc");
	}
      if(!fastfinder->setRawEvent(rawEvent)) continue;

      //Fill the TClonesArray
      arr_tracklets.Clear();
      std::list<Tracklet>& rec_tracklets = fastfinder->getFinalTracklets();
      if(rec_tracklets.empty()) continue;

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
	  recEvent->insertTrack(recTrack);
#endif
	}

#ifdef _ENABLE_KF
      std::list<KalmanTrack>& rec_tracks = fastfinder->getKalmanTracks();
      for(std::list<KalmanTrack>::iterator iter = rec_tracks.begin(); iter != rec_tracks.end(); ++iter)
	{
	  //iter->print();
	  SRecTrack recTrack = iter->getSRecTrack();
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
  cout << "kFastTracking ends successfully." << endl;

  saveFile->cd();
  saveTree->Write();
  saveFile->Close();

  delete fastfinder;
  delete triggerAna;
#ifdef _ENABLE_KF
  filter->close();
#endif

  return 1;
}
