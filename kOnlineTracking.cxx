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
#include "MySQLSvc.h"
#include "JobOptsSvc.h"
#include "TriggerAnalyzer.h"

#include "MODE_SWITCH.h"

using namespace std;

int main(int argc, char *argv[])
{
  //Initialize job options
  JobOptsSvc* jobOptsSvc = JobOptsSvc::instance();
  jobOptsSvc->init(argv[1]);

  //Initialize geometry service with calibration
  GeomSvc* geometrySvc = GeomSvc::instance();
  geometrySvc->loadCalibration(jobOptsSvc->m_calibrationsFile);

  //Initialize MySQL service and connect to database, e906-db1 by default
  MySQLSvc* p_mysqlSvc = MySQLSvc::instance();
  p_mysqlSvc->setUserPasswd("production", "qqbar2mu+mu-");
  p_mysqlSvc->connect();
  p_mysqlSvc->setWorkingSchema(jobOptsSvc->m_inputFile);
  p_mysqlSvc->bookOutputTables();

  //Data output definition
  int nTracklets;
  SRawEvent* rawEvent = new SRawEvent();
  SRecEvent* recEvent = new SRecEvent();
  TClonesArray* tracklets = new TClonesArray("Tracklet");
  TClonesArray& arr_tracklets = *tracklets;

  /*
  TFile* saveFile = new TFile(jobOptsSvc->m_outputFile.c_str(), "recreate");
  TTree* saveTree = new TTree("save", "save");

  saveTree->Branch("nTracklets", &nTracklets, "nTracklets/I");
  saveTree->Branch("rawEvent", &rawEvent, 256000, 99);
  saveTree->Branch("recEvent", &recEvent, 256000, 99);
  saveTree->Branch("tracklets", &tracklets, 256000, 99);
  tracklets->BypassStreamer();
  */

  //Initialize track finder
  KalmanFastTracking* fastfinder = new KalmanFastTracking(false);
  VertexFit* vtxfit  = new VertexFit();

  //Initialize the trigger analyzer
  TriggerAnalyzer* triggerAna = new TriggerAnalyzer();
  if(jobOptsSvc->m_enableTriggerMask)
    {
      triggerAna->init();
      triggerAna->buildTriggerTree();
    }

  //Quality control numbers and plots
  int nEvents_loaded = 0;
  int nEvents_tracked = 0;
  int nEvents_dimuon = 0;
  int nEvents_dimuon_real = 0;

  //Start tracking
  int nEvents = p_mysqlSvc->getNEvents();
  cout << "There are " << nEvents << " events in " << jobOptsSvc->m_inputFile << endl;
  for(int i = 0; i < nEvents; ++i) 
    {
      //Read data
      if(!p_mysqlSvc->getNextEvent(rawEvent)) continue;
      ++nEvents_loaded;

      //Do the tracking
      cout << "\r Tracking runID = " << rawEvent->getRunID() << " eventID = " << rawEvent->getEventID() << ", " << (i+1)*100/nEvents << "% finished. ";
      cout << nEvents_tracked*100/nEvents_loaded << "% have at least one track, " << nEvents_dimuon*100/nEvents_loaded << "% have at least one dimuon pair, ";
      cout << nEvents_dimuon_real*100/nEvents_loaded << "% have successful dimuon vertex fit.";

      if(jobOptsSvc->m_enableTriggerMask)
	{
	  triggerAna->trimEvent(rawEvent);
	  rawEvent->reIndex("aoct");
	}
      else
	{
	  rawEvent->reIndex("oac");
	}
      if(!fastfinder->setRawEvent(rawEvent)) continue;
      ++nEvents_tracked;

      //Output
      arr_tracklets.Clear();
      std::list<Tracklet>& rec_tracklets = fastfinder->getFinalTracklets();
      if(rec_tracklets.empty()) continue;

      recEvent->setRawEvent(rawEvent);
      nTracklets = 0;
      int nPos = 0;
      int nNeg = 0;
      for(std::list<Tracklet>::iterator iter = rec_tracklets.begin(); iter != rec_tracklets.end(); ++iter)
	{
	  //iter->print();
	  iter->calcChisq();
	  new(arr_tracklets[nTracklets++]) Tracklet(*iter);

	  SRecTrack recTrack = iter->getSRecTrack();
	  recEvent->insertTrack(recTrack);

	  iter->getCharge() > 0 ? ++nPos : ++nNeg;
	}
      if(nPos > 0 && nNeg > 0) ++nEvents_dimuon;
 
      //Perform dimuon vertex fit 
      if(vtxfit->setRecEvent(recEvent)) ++nEvents_dimuon_real;

      if(nTracklets > 0)
	{
	  p_mysqlSvc->writeTrackingRes(recEvent, tracklets);
  	  //saveTree->Fill();
	}	 
      rawEvent->clear();
      recEvent->clear();
    }
  cout << endl;
  cout << "kOnlineTracking ended successfully." << endl;
  cout << "In total " << nEvents_loaded << " events loaded from " << argv[1] << ": " << nEvents_tracked << " events have at least one track, ";
  cout << nEvents_dimuon << " events have at least one dimuon pair, ";
  cout << nEvents_dimuon_real << " events have successful dimuon vertex fit." << endl;

  /*
  saveFile->cd();
  saveTree->Write();
  saveFile->Close();
  */

  delete fastfinder;
  delete vtxfit;
  delete triggerAna;

  return 1;
}
