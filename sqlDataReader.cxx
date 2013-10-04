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

#include "SRawEvent.h"
#include "GeomSvc.h"

using namespace std;

int main(int argc, char **argv)
{
  cout << "Exporting Run: " << argv[1] << " to ROOT file: " << argv[2] << endl;

  ///Initialize the geometry service and output file 
  GeomSvc *p_geomSvc = GeomSvc::instance();
  p_geomSvc->init("geometry_R997");

  SRawEvent *event = new SRawEvent();

  TFile *saveFile = new TFile(argv[2], "recreate");
  TTree *saveTree = new TTree("save", "save");

  saveTree->Branch("rawEvent", &event, 256000, 99);

  ///Connect to the SQL databases
  //TSQLServer *con = TSQLServer::Connect("mysql://e906-gat2.fnal.gov", "seaguest","qqbar2mu+mu-");
  TSQLServer *con = TSQLServer::Connect("mysql://seaquel.physics.illinois.edu", "seaguest","qqbar2mu+mu-");
  char query[300]; 
  const char *buf = "SELECT runID,spillID,eventID,hitID,elementID,tdcTime,driftTime,driftDistance,detectorName,"
    "inTime,masked,wirePosition FROM %s.Hit WHERE (detectorName LIKE 'D%%' OR detectorName LIKE 'H%%' OR detectorName LIKE 'P%%')"
    " ORDER BY eventID,inTime DESC LIMIT 10000";
  sprintf(query, buf, argv[1]); cout << query << endl;
  TSQLResult *res = con->Query(query);

  UInt_t nEntries = res->GetRowCount();
  cout << "Totally " << nEntries << " hits in this run" << endl;
  Int_t eventID_curr, eventID;
  eventID_curr = -1;
  for(UInt_t i = 0; i < nEntries; i++)
    {
      if(i % 1000 == 0)
	{
	  cout << "Converting hit " << i << endl;
	}

      TSQLRow *row = res->Next();

      eventID = atoi(row->GetField(2));
      if(eventID_curr < eventID || i == nEntries - 1)
	{
	  if(eventID_curr >= 0)
	    {
	      event->setEventInfo(atoi(row->GetField(0)), atoi(row->GetField(1)), eventID_curr);
              event->reIndex();
	      saveTree->Fill();
	    }

	  event->clear();
	  eventID_curr = eventID;
	}

      string detectorName(row->GetField(8));
      Int_t elementID = atoi(row->GetField(4));
      p_geomSvc->toLocalDetectorName(detectorName, elementID);
      
      Hit h;

      h.index = atoi(row->GetField(3));
      h.detectorID = p_geomSvc->getDetectorID(detectorName);
      h.elementID = elementID;
      h.tdcTime = atof(row->GetField(5));
      h.inTime = atoi(row->GetField(9));
      
      if(row->GetField(6) != NULL)
	{
	  h.driftTime = atof(row->GetField(6));
	}
      else
	{
	  h.driftTime = 0.;
	}

      if(row->GetField(7) != NULL)
	{
	  h.driftDistance = atof(row->GetField(7));
	}
      else
	{
	  h.driftDistance = 0.;
	}

      if(row->GetField(10) != NULL)
	{
	  h.hodoMask = atoi(row->GetField(10)); 
	}	  
      else
	{
	  h.hodoMask = 1;
	}

      if(row->GetField(11) != NULL)
	{
	  h.pos = atof(row->GetField(11));
	}
      else
	{
	  h.pos = p_geomSvc->getMeasurement(h.detectorID, h.elementID);
	}
      
      event->insertHit(h);
      delete row;
    }

  delete res;
  delete con;

  saveFile->cd();
  saveTree->Write();
  saveFile->Close();

  return 1;
}
