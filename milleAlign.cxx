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
#include <TMath.h>
#include <TF1.h>
#include <TH1D.h>

#include "GeomSvc.h"
#include "SRawEvent.h"
#include "KalmanUtil.h"
#include "KalmanTrack.h"
#include "KalmanFilter.h"
#include "KalmanFinder.h"
#include "KalmanFitter.h"
#include "VertexFit.h"
#include "SRecEvent.h"
#include "SMillepede.h"
#include "SMillepedeUtil.h"

#include "MODE_SWITCH.h"

int main(int argc, char *argv[])
{
  //Initialize geometry service
  Log("Initializing geometry service ... ");
  GeomSvc* geometrySvc = GeomSvc::instance();
  geometrySvc->init(GEOMETRY_VERSION);

  //Retrieve the raw event
  Log("Retrieving the event stored in ROOT file ... ");
#ifdef _ENABLE_KF
  SRawEvent* rawEvent = new SRawEvent();
  SRecEvent* recEvent = new SRecEvent();
#else
  TClonesArray* tracklets = new TClonesArray("Tracklet");
  tracklets->Clear();
#endif

  TFile *dataFile = new TFile(argv[1], "READ");
  TTree *dataTree = (TTree *)dataFile->Get("save");

#ifdef _ENABLE_KF
  dataTree->SetBranchAddress("rawEvent", &rawEvent);
  dataTree->SetBranchAddress("recEvent", &recEvent);
#else
  dataTree->SetBranchAddress("tracklets", &tracklets);
#endif

  //Initialize track finder
  Log("Initializing the millepede ... ");
  SMillepede* mille = new SMillepede();
  mille->init();
  if(argc > 4)
    {
      mille->bookEvaluationTree(argv[4]);
    }
  else
    {
      mille->bookEvaluationTree("align_eval.root");
    }

  int nEvtMax = dataTree->GetEntries();
  for(int i = 0; i < nEvtMax; i++)
    {    
      dataTree->GetEntry(i);

#ifdef _ENABLE_KF
      mille->setEvent(rawEvent, recEvent);
#else
      mille->setEvent(tracklets);
      tracklets->Clear();
#endif
    }

  mille->fitAlignment();
  mille->printQAPlots();
  mille->printResults(argv[2], argv[3]);

  delete mille;

  return 1;
}
