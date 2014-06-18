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

#include "GeomSvc.h"
#include "SRawEvent.h"
#include "KalmanUtil.h"
#include "KalmanTrack.h"
#include "KalmanFilter.h"
#include "KalmanFitter.h"
#include "VertexFit.h"
#include "SRecEvent.h"
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
  SRecEvent* recEvent = new SRecEvent();

  TFile* dataFile = new TFile(jobOptsSvc->m_inputFile.c_str(), "READ");
  TTree* dataTree = (TTree*)dataFile->Get("save");

  dataTree->SetBranchAddress("recEvent", &recEvent);

  TFile* saveFile = new TFile(jobOptsSvc->m_outputFile.c_str(), "recreate");
  TTree* saveTree = dataTree->CloneTree(0);
  
  //Initialize track finder
  LogInfo("Initializing the track finder and kalman filter ... ");
  VertexFit* vtxfit = new VertexFit();
  vtxfit->enableOptimization();
  if(jobOptsSvc->m_enableEvaluation) 
    {
      string evalFileName = "eval_" + jobOptsSvc->m_outputFile;
      vtxfit->bookEvaluation(evalFileName.c_str());
    }

  const int offset = jobOptsSvc->m_firstEvent;
  int nEvtMax = jobOptsSvc->m_nEvents > 0 ? jobOptsSvc->m_nEvents + offset : dataTree->GetEntries();
  if(nEvtMax > dataTree->GetEntries()) nEvtMax = dataTree->GetEntries();
  LogInfo("Running from event " << offset << " through to event " << nEvtMax);
  for(int i = offset; i < nEvtMax; ++i)
    {
      dataTree->GetEntry(i);
      cout << "\r Processing event " << i << " with eventID = " << recEvent->getEventID() << ", ";
      cout << (i + 1)*100/nEvtMax << "% finished .. " << flush;

      vtxfit->setRecEvent(recEvent);

      if(recEvent->getNDimuons() > 0) saveTree->Fill();
      recEvent->clear();
    }
  cout << endl;
  cout << "kVertex ends successfully." << endl;

  saveFile->cd();
  saveTree->Write();
  saveFile->Close();

  delete vtxfit;

  return 1;
}
