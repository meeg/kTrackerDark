#include <iostream>
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

using namespace std;

int main(int argc, char **argv)
{
  cout << "Exporting Run: " << argv[1] << " to ROOT file: " << argv[2] << endl;

  ///Initialize the geometry service and output file 
  GeomSvc* p_geomSvc = GeomSvc::instance();
  p_geomSvc->init(GEOMETRY_VERSION);

  MySQLSvc* p_mysqlSvc = MySQLSvc::instance();
  p_mysqlSvc->connect();
  p_mysqlSvc->setWorkingSchema(argv[1]);

  SRawEvent* rawEvent = new SRawEvent();

  TFile *saveFile = new TFile(argv[2], "recreate");
  TTree *saveTree = new TTree("save", "save");

  saveTree->Branch("rawEvent", &rawEvent, 256000, 99);

  int nEvents = p_mysqlSvc->getNEventsFast();
  cout << "Totally " << nEvents << " events in this run" << endl;
  for(int i = 0; i < nEvents; ++i)
    {
      if(!p_mysqlSvc->getNextEvent(rawEvent)) continue;
      if(rawEvent->getEventID() % 1000 == 0)
	{
	  cout << "Converting event " << rawEvent->getEventID() << endl;
	  //rawEvent->print();
	}

      saveTree->Fill();
    }

  saveFile->cd();
  saveTree->Write();
  saveFile->Close();

  delete p_mysqlSvc;
  delete p_geomSvc;

  return 1;
}
