#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TLorentzVector.h>

#include "SRawEvent.h"
#include "GeomSvc.h"
#include "MySQLSvc.h"
#include "JobOptsSvc.h"

using namespace std;

int main(int argc, char **argv)
{
  cout << "Exporting Run: " << argv[1] << " to ROOT file: " << argv[2] << endl;

  ///Initialize the job option service
  JobOptsSvc* p_jobOptsSvc = JobOptsSvc::instance();

  ///Initialize the geometry service and output file 
  GeomSvc* p_geomSvc = GeomSvc::instance();
  p_geomSvc->loadCalibration(p_jobOptsSvc->m_calibrationsFile);

  MySQLSvc* p_mysqlSvc = MySQLSvc::instance();
  if(argc > 4)
    {
      p_mysqlSvc->connect(argv[3], atoi(argv[4]));
    }
  else
    {
      p_mysqlSvc->connect();
    }
  p_mysqlSvc->setWorkingSchema(argv[1]);

  SRawEvent* rawEvent = new SRawEvent();

  TFile *saveFile = new TFile(argv[2], "recreate");
  TTree *saveTree = new TTree("save", "save");

  saveTree->Branch("rawEvent", &rawEvent, 256000, 99);

  int nEvents = p_mysqlSvc->getNEventsFast();
  cout << "Totally " << nEvents << " events in this run" << endl;
  
  if(argc > 5) nEvents = atoi(argv[5]);
  for(int i = 0; i < nEvents; ++i)
    {
      if(!p_mysqlSvc->getNextEvent(rawEvent)) continue;
      cout << "\r Converting event " << rawEvent->getEventID() << ", " << (i+1)*100/nEvents << "% finished." << flush;

      saveTree->Fill();
    }
  cout << endl;
  cout << "sqlDataReader ends successfully." << endl;

  saveFile->cd();
  saveTree->Write();
  saveFile->Close();

  delete p_mysqlSvc;
  delete p_geomSvc;
  delete p_jobOptsSvc;

  return 1;
}
