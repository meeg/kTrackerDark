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

#define OUT_TO_SCREEN

class MySQLSvc
{
public:
  MySQLSvc();
  static MySQLSvc* instance();
  
  //Connect to the server
  bool connect(std::string sqlServer = "seaquel.physics.illinois.edu");

  //check if the new event is available
  bool isNewEvtAvailable();

  //check if the run is stopped
  bool isRunStopped();

  //Get run info
  int getNEventsFast();
  int getNEvents();

  //Get the latest event
  bool getEventByID(SRawEvent* rawEvent, int eventID);
  bool getLatestEvt(SRawEvent* rawEvent);
  bool getRandomEvt(SRawEvent* rawEvent);

  //Output to database/txt file/screen
  void writeTrackingRes(SRecEvent* recEvent, TClonesArray* tracklets);
  void writeTrackTable(int trackID, SRecTrack* recTrack);
  void writeTrackHitTable(int trackID, Tracklet* tracklet);
  void writeDimuonTable(int dimuonID, int idx_positive, int idx_negative);

  //Set the data schema
  void setWorkingSchema(std::string schema) { dataSchema = schema; } 
  void setLoggingSchema(std::string schema) { logSchema = schema; } 

  //Memory-safe sql queries
  bool makeQuery();
  bool nextEntry();

private:
  //pointer to the only instance
  static MySQLSvc* p_mysqlSvc;

  //pointer to the geometry service
  GeomSvc* p_geomSvc;

  //SQL server
  TSQLServer* server;
  TSQLResult* res;
  TSQLRow* row;

  //Random generator
  TRandom rndm;

  //Last eventID used
  int eventID_last;

  //run-related info
  int runID;
  int spillID;
  int nEvents;

  //Query string used in all clause
  char query[500];

  //name of the production schema working on
  std::string dataSchema;
  std::string logSchema;

  //Internal counter of tracks and dimuons
  int nTracks;
  int nDimuons;

  std::vector<TLorentzVector> mom_vertex;
  std::vector<TVector3> pos_vertex;
};

#endif
