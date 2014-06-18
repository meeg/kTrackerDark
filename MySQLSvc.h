/*
IO manager to handle fast extraction of data from database or upload data to database

Author: Kun Liu, liuk@fnal.gov
Created: 2013.9.29
*/

#ifndef _MYSQLSVC_H
#define _MYSQLSVC_H

#include <iostream>
#include <vector>
#include <string>
#include <list>
#include <algorithm>

#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TRandom.h>
#include <TClonesArray.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include "GeomSvc.h"
#include "SRawEvent.h"
#include "FastTracklet.h"
#include "TriggerAnalyzer.h"
#include "JobOptsSvc.h"

//#define OUT_TO_SCREEN
//#define USE_M_TABLES

class MySQLSvc
{
public:
  MySQLSvc();
  ~MySQLSvc();
  static MySQLSvc* instance();
  
  //Connect to the server
  bool connect(std::string mysqlServer = "", int mysqlPort = -1);

  //Set username/password
  void setUserPasswd(std::string user_input, std::string passwd_input) { user = user_input; passwd = passwd_input; }

  //check if the run is stopped
  bool isRunStopped();

  //Get run info
  int getNEventsFast();
  int getNEvents();

  //Gets
  bool getEvent(SRawEvent* rawEvent, int eventID);
  bool getLatestEvt(SRawEvent* rawEvent);
  bool getRandomEvt(SRawEvent* rawEvent);
  bool getNextEvent(SRawEvent* rawEvent);
  bool getNextEvent(SRawMCEvent* rawEvent);

  //Check if the event has been loaded
  bool isEventLoaded(int eventID) { return std::find(eventIDs.begin(), eventIDs.end(), eventID) != eventIDs.end(); } 

  //Get the event header
  bool getEventHeader(SRawEvent* rawEvent, int eventID);
  bool getMCGenInfo(SRawMCEvent* mcEvent, int eventID);

  //Output to database/txt file/screen
  void bookOutputTables();
  void writeTrackingRes(SRecEvent* recEvent, TClonesArray* tracklets);
  void writeTrackTable(int trackID, SRecTrack* recTrack);
  void writeTrackHitTable(int trackID, Tracklet* tracklet);
  void writeDimuonTable(int dimuonID, SRecDimuon dimuon);

  //Set the data schema
  void setWorkingSchema(std::string schema);
  void setLoggingSchema(std::string schema) { logSchema = schema; } 
  void enableQIE(bool opt) { readQIE = opt; }
  void enableTargetPos(bool opt) { readTargetPos = opt; }
  void enableTriggerHits(bool opt) { readTriggerHits = opt; }

  //Memory-safe sql queries
  int makeQuery();
  bool nextEntry();
  
  int getInt(int id, int default_val = 0);
  double getDouble(int id, double default_val = 0.);
  std::string getString(int id, std::string default_val = "");

private:
  //pointer to the only instance
  static MySQLSvc* p_mysqlSvc;

  //Username and password
  std::string user;
  std::string passwd;

  //pointer to the geometry service
  GeomSvc* p_geomSvc;

  //pointer to trigger analyzer
  TriggerAnalyzer* p_triggerAna;

  //SQL server
  TSQLServer* server;
  TSQLResult* res;
  TSQLRow* row;

  //Test if QIE/TriggerHits table exists
  bool readQIE;
  bool readTriggerHits;
  bool readTargetPos;
  bool setTriggerEmu;

  //Random generator
  TRandom rndm;

  //run-related info
  int runID;
  int spillID;
  int nEvents;
  std::list<int> eventIDs;

  //Query string used in all clause
  char query[2000];

  //name of the production schema working on
  std::string dataSchema;
  std::string logSchema;

  //Internal counter of tracks and dimuons
  int nTracks;
  int nDimuons;
};

#endif
