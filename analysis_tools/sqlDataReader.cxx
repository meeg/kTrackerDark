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
#include <TStopwatch.h>

#include "SRawEvent.h"
#include "GeomSvc.h"
#include "MySQLSvc.h"
#include "JobOptsSvc.h"
#include "ThresholdSvc.h"

using namespace std;
using Threshold::live;

int main(int argc, char **argv)
{
  cout << "Exporting Run: " << argv[1] << " to ROOT file: " << argv[2] << endl;
  TStopwatch timer;
  timer.Start();


  //NOTE:
  //Comment out this line if you want this to print at every event.
  //This is a hack until this program takes a job options file as the argument
  Threshold::ThresholdSvc::Get().SetLive( false );

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
  p_mysqlSvc->initReader();

  SRawEvent* rawEvent = new SRawEvent();

  TFile *saveFile = new TFile(argv[2], "recreate");
  TTree *saveTree = new TTree("save", "save");

  saveTree->Branch("rawEvent", &rawEvent, 256000, 99);

  //will force save every 1000 events
  int saveFreq = 1000;

  int nEvents = p_mysqlSvc->getNEventsFast();
  cout << "Totally " << nEvents << " events in this run" << endl;
  
  //if user specified number of events, then do that many
  //   but make sure they can't request more than exist
  if(argc > 5) nEvents = std::min( nEvents, atoi(argv[5]) );

  //plan to print progress at each %1 (unless in live mode)
  const int printFreq = (nEvents/100);

  timer.Start();
  for(int i = 0; i < nEvents; ++i)
    {
      if(!p_mysqlSvc->getNextEvent(rawEvent)) continue;

      //print progress every event or every 1%
      const int fracDone = (i+1)*100/nEvents;
      if( live() ) 
        cout << "\r Converting event " << rawEvent->getEventID() << ", " << fracDone << "% finished." << flush;
      else if( 0 == i % printFreq ) 
        {
          timer.Stop();
          cout << Form( "Converting Event %d, %d%% finished.  Time to process last %d events shown below:", rawEvent->getEventID(), fracDone, printFreq ) << endl;
          timer.Print();
          timer.Start();
        }

      saveTree->Fill();

      //save every n events to avoid large losses
      if(0 == i % saveFreq)
        saveTree->AutoSave("SaveSelf");
    }
  cout << endl;
  cout << "sqlDataReader ends successfully." << endl;

  saveFile->cd();
  saveTree->Write();
  saveFile->Close();

  delete p_mysqlSvc;
  delete p_geomSvc;
  delete p_jobOptsSvc;

  return EXIT_SUCCESS;
}
