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
#include <TClonesArray.h>

#include "SRecEvent.h"
#include "JobOptsSvc.h"
#include "GeomSvc.h"
#include "MySQLSvc.h"

using namespace std;

int main(int argc, char **argv)
{
  cout << "Uploading file: " << argv[1] << " to sql schema " << argv[2] << endl;

  ///Initialize job option service
  JobOptsSvc* p_jobOptsSvc = JobOptsSvc::instance();

  ///Initialize the geometry service and output file 
  GeomSvc* p_geomSvc = GeomSvc::instance();

  ///Intialize the mysql service
  MySQLSvc* p_mysqlSvc = MySQLSvc::instance();
  p_mysqlSvc->setUserPasswd("production", "qqbar2mu+mu-");
  p_mysqlSvc->connectOutput();
  p_mysqlSvc->setOutputSchema(argv[2]);
  if(!p_mysqlSvc->initWriter()) exit(EXIT_FAILURE);

  ///Retrieve data from file
  TClonesArray* tracklets = new TClonesArray("Tracklet");
  SRecEvent* recEvent = new SRecEvent();

  TFile* dataFile = new TFile(argv[1], "READ");
  TTree* dataTree = (TTree*)dataFile->Get("save");

  dataTree->SetBranchAddress("recEvent", &recEvent);
  dataTree->SetBranchAddress("tracklets", &tracklets);

  int nEvents = dataTree->GetEntries();
  for(int i = 0; i < nEvents; ++i)
    {
      dataTree->GetEntry(i);
      cout << "\r Uploading event " << recEvent->getEventID() << ", " << (i+1)*100/nEvents << "% finished." << flush;

      p_mysqlSvc->writeTrackingRes(recEvent, tracklets);
    }
  cout << endl;
  cout << "sqlResWriter ends successfully." << endl;

  delete p_mysqlSvc;
  delete p_geomSvc;
  delete p_jobOptsSvc;

  return EXIT_SUCCESS;
}
