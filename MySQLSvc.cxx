/*
MySQLSvc.cxx
	    
Implementation of class MySQLSvc.
	       	         
Author: Kun Liu, liuk@fnal.gov		    
Created: 9-29-2013
*/

#include <TLorentzVector.h>

#include "FastTracklet.h"
#include "MySQLSvc.h"

MySQLSvc* MySQLSvc::p_mysqlSvc = NULL;

MySQLSvc::MySQLSvc()
{
  runID = -1;
  spillID = -1;
  eventID_last = 0;
  
  server = NULL;
  res = NULL;
  row = NULL;

  p_geomSvc = GeomSvc::instance();

  dataSchema = "run_002173_R003";
  logSchema = "log";

  nTracks = 0;
  nDimuons = 0;

  rndm.SetSeed(0);
}

MySQLSvc* MySQLSvc::instance()
{
  if(p_mysqlSvc == NULL)
    {
      p_mysqlSvc = new MySQLSvc();
    }

  return p_mysqlSvc;
}

bool MySQLSvc::connect(std::string sqlServer)
{
  char address[300];
  sprintf(address, "mysql://%s", sqlServer.c_str());

  server = TSQLServer::Connect(address, "seaguest", "qqbar2mu+mu-");
  
  if(server == NULL)
    {
      Log("Connection to database " << sqlServer.c_str() << " failed!");
      return false;
    }
  return true;
}

void MySQLSvc::setWorkingSchema(std::string schema)
{
  dataSchema = schema;
  sprintf(query, "USE %s", dataSchema.c_str());

  server->Exec(query);
}

bool MySQLSvc::isNewEvtAvailable()
{
  sprintf(query, "SELECT eventID FROM Event ORDER BY eventID DESC LIMIT 2");
  if(makeQuery() < 2) return false;
  
  nextEntry(); nextEntry();
  int eventID_curr = getInt(row->GetField(0));
  if(eventID_curr == eventID_last) return false;
    
  eventID_last = eventID_curr;
  return true;
}

bool MySQLSvc::isRunStopped()
{
  sprintf(query, "SELECT productionEnd from %s.production WHERE runID=%d", logSchema.c_str(), runID);
  if(makeQuery() != 1) return false;

  nextEntry();
  if(row->GetField(0) == NULL)
    {
      return true;
    }
  return false;
}

bool MySQLSvc::getLatestEvt(SRawEvent* rawEvent)
{
  if(!isNewEvtAvailable()) return false;

  rawEvent->clear();
  if(!getEventHeader(rawEvent, eventID_last))
    {
      return false;
    }
  return getEvent(rawEvent, eventID_last);
}

bool MySQLSvc::getRandomEvt(SRawEvent* rawEvent)
{
  eventID_last = int(rndm.Rndm()*getNEventsFast());
  
  rawEvent->clear();
  if(!getEventHeader(rawEvent, eventID_last))
    {
      return false;
    }
  return getEvent(rawEvent, eventID_last);
}

bool MySQLSvc::getNextEvent(SRawEvent* rawEvent)
{
  ++eventID_last;
  
  rawEvent->clear();
  if(!getEventHeader(rawEvent, eventID_last))
    {
      return false;
    }
  return getEvent(rawEvent, eventID_last);
}

bool MySQLSvc::getNextEvent(SRawMCEvent* mcEvent)
{
  ++eventID_last;

  mcEvent->clear();
  if(!getEventHeader(mcEvent, eventID_last))
    {
      return false;
    }
  return getEvent(mcEvent, eventID_last);
}

bool MySQLSvc::getEvent(SRawEvent* rawEvent, int eventID)
{
#ifdef USE_M_TABLES
  sprintf(query, "SELECT mHitID,elementID,tdcTime,driftTime,driftDistance,detectorName,1,1 FROM mHit WHERE (detectorName LIKE 'D%%'"
	  " OR detectorName LIKE 'H%%' OR detectorName LIKE 'P%%') AND eventID=%d", eventID);
#else
  sprintf(query, "SELECT hitID,elementID,tdcTime,driftTime,driftDistance,detectorName,inTime,masked FROM Hit WHERE (detectorName LIKE 'D%%'"
	  " OR detectorName LIKE 'H%%' OR detectorName LIKE 'P%%') AND eventID=%d", eventID);
#endif
  if(makeQuery() == 0) return false;

  int nHits = res->GetRowCount();
  if(nHits < 1) return false;

  for(int i = 0; i < nHits; ++i)
    {
      nextEntry();

      std::string detectorName(row->GetField(5));
      int elementID = getInt(row->GetField(1));
      p_geomSvc->toLocalDetectorName(detectorName, elementID);
       
      Hit h;
      h.index = getInt(row->GetField(0));
      h.detectorID = p_geomSvc->getDetectorID(detectorName);
      h.elementID = elementID;
      h.tdcTime = getDouble(row->GetField(2));
      h.inTime = getInt(row->GetField(6), 1);
      h.pos = p_geomSvc->getMeasurement(h.detectorID, h.elementID);
      h.driftTime = getDouble(row->GetField(3));
      h.driftDistance = getDouble(row->GetField(4));
      h.hodoMask = getInt(row->GetField(7), 1);
      
      rawEvent->insertHit(h);
    }

  eventID_last = eventID;
  rawEvent->reIndex();
  return true;
}

int MySQLSvc::getNEvents()
{
  sprintf(query, "SELECT MAX(eventID) FROM Event");
  if(makeQuery() != 1) return 0;

  nextEntry();
  int nTotal = getInt(row->GetField(0));
  return nTotal;
}

bool MySQLSvc::getEventHeader(SRawEvent* rawEvent, int eventID)
{
  //Get the event header
  sprintf(query, "SELECT runID,spillID,NIM1,NIM2,NIM3,NIM4,NIM5,MATRIX1,MATRIX2,MATRIX3,MATRIX4,MATRIX5 FROM Event WHERE eventID=%d", eventID);
  if(makeQuery() != 1) return false;

  nextEntry();
  runID = getInt(row->GetField(0));
  spillID = getInt(row->GetField(1));
  rawEvent->setEventInfo(runID, spillID, eventID);

  int triggers[10];
  for(int i = 0; i < 10; ++i)
    {
      triggers[i] = getInt(row->GetField(i+2));
    }
  rawEvent->setTriggerBits(triggers);

  return true;
}

bool MySQLSvc::getEventHeader(SRawMCEvent* mcEvent, int eventID)
{
  sprintf(query, "SELECT mTrackID1,mTrackID2,sigWeight,mass,xF,xB,xT,dx,dy,dz,dpx,dpy,runID,spillID FROM mDimuon WHERE acceptHodoAll=1 AND acceptDriftAll=1 AND eventID=%d", eventID);
  if(makeQuery() != 1) return false;
  nextEntry();

  runID = getInt(row->GetField(12));
  spillID = getInt(row->GetField(13));
  mcEvent->setEventInfo(runID, spillID, eventID);

  int trackID[2] = { getInt(row->GetField(0)), getInt(row->GetField(1)) };
  mcEvent->weight = getDouble(row->GetField(2));
  mcEvent->mass = getDouble(row->GetField(3));
  mcEvent->xF = getDouble(row->GetField(4));
  mcEvent->x1 = getDouble(row->GetField(5));
  mcEvent->x2 = getDouble(row->GetField(6));
  mcEvent->vtx.SetXYZ(getDouble(row->GetField(7)), getDouble(row->GetField(8)), getDouble(row->GetField(9)));

  double px = getDouble(row->GetField(10));
  double py = getDouble(row->GetField(11));
  mcEvent->pT = sqrt(px*px + py*py);

  for(int i = 0; i < 2; ++i)
    {
      //At vertex
      sprintf(query, "SELECT px0,py0,pz0 FROM mTrack WHERE mTrackID=%d", trackID[i]);
      if(makeQuery() != 1) return false;
      
      nextEntry();
      mcEvent->p_vertex[i].SetXYZ(getDouble(row->GetField(0)), getDouble(row->GetField(1)), getDouble(row->GetField(2)));
    
      //At station 1,2,3,4
      sprintf(query, "SELECT hpx,hpy,hpz,hx,hy,hz FROM mGeantHit WHERE geantName RLIKE 'O[1-4]' AND mTrackID=%d", trackID[i]);
      if(makeQuery() != 4) return false;

      nextEntry();
      mcEvent->p_station1[i].SetXYZ(getDouble(row->GetField(0)), getDouble(row->GetField(1)), getDouble(row->GetField(2)));
      mcEvent->v_station1[i].SetXYZ(getDouble(row->GetField(3)), getDouble(row->GetField(4)), getDouble(row->GetField(5)));
    
      nextEntry();
      mcEvent->p_station2[i].SetXYZ(getDouble(row->GetField(0)), getDouble(row->GetField(1)), getDouble(row->GetField(2)));
      mcEvent->v_station2[i].SetXYZ(getDouble(row->GetField(3)), getDouble(row->GetField(4)), getDouble(row->GetField(5)));

      nextEntry();
      mcEvent->p_station3[i].SetXYZ(getDouble(row->GetField(0)), getDouble(row->GetField(1)), getDouble(row->GetField(2)));
      mcEvent->v_station3[i].SetXYZ(getDouble(row->GetField(3)), getDouble(row->GetField(4)), getDouble(row->GetField(5)));

      nextEntry();
      mcEvent->p_station4[i].SetXYZ(getDouble(row->GetField(0)), getDouble(row->GetField(1)), getDouble(row->GetField(2)));
      mcEvent->v_station4[i].SetXYZ(getDouble(row->GetField(3)), getDouble(row->GetField(4)), getDouble(row->GetField(5)));
    }

  return true;
}

void MySQLSvc::bookOutputTables()
{
  //Clear all the tables if exist
  std::string tableNames[3] = {"kTrack", "kHit", "kDimuon"};
  for(int i = 0; i < 3; ++i)
    {
      sprintf(query, "DROP TABLE IF EXISTS %s", tableNames[i].c_str());
#ifndef OUT_TO_SCREEN
      server->Exec(query);
#else
      std::cout << __FUNCTION__ << ": " << query << std::endl;
#endif
    }

  //Book kTrack table
  sprintf(query, "CREATE TABLE kTrack ("
	  "trackID     INTEGER,"
	  "runID       INTEGER,"
	  "spillID     INTEGER,"
	  "eventID     INTEGER,"
	  "charge      INTEGER,"
	  "numHits     INTEGER,"
	  "chisq       DOUBLE, "
	  "x0          DOUBLE, "
	  "y0          DOUBLE, "
	  "z0          DOUBLE, "
	  "px0         DOUBLE, "
	  "py0         DOUBLE, "
	  "pz0         DOUBLE, "
	  "x1          DOUBLE, "
	  "y1          DOUBLE, "
	  "z1          DOUBLE, "
	  "px1         DOUBLE, "
	  "py1         DOUBLE, "
	  "pz1         DOUBLE, "
	  "x3          DOUBLE, "
	  "y3          DOUBLE, "
	  "z3          DOUBLE, "
	  "px3         DOUBLE, "
	  "py3         DOUBLE, "
	  "pz3         DOUBLE, "
          "PRIMARY KEY(runID, spillID, eventID), "
	  "INDEX(eventID), INDEX(charge))");
#ifndef OUT_TO_SCREEN
  server->Exec(query);
#else
  std::cout << __FUNCTION__ << ": " << query << std::endl;
#endif

  //Book kHit table
  sprintf(query, "CREATE TABLE kHit ("
	  "runID       SMALLINT,"
	  "eventID     INTEGER, "
	  "trackID     INTEGER, "
	  "hitID       BIGINT,  "
	  "driftSign   SMALLINT,"
	  "residual    DOUBLE,  "
	  "PRIMARY KEY(runID, trackID, hitID))");
#ifndef OUT_TO_SCREEN
  server->Exec(query);
#else
  std::cout << __FUNCTION__ << ": " << query << std::endl;
#endif

  //Bool kDimuon table
  sprintf(query, "CREATE TABLE kDimuon ("
	  "dimuonID    INTEGER,"
	  "runID       INTEGER,"
	  "spillID     INTEGER,"
	  "eventID     INTEGER,"
	  "posTrackID  INTEGER,"
	  "negTrackID  INTEGER,"
	  "dx          DOUBLE, "
	  "dy          DOUBLE, "
	  "dz          DOUBLE, "
	  "dpx         DOUBLE, "
	  "dpy         DOUBLE, "
	  "dpz         DOUBLE, "
	  "mass        DOUBLE, "
	  "xF          DOUBLE, "
	  "xB          DOUBLE, "
	  "xT          DOUBLE, "
	  "trackSeparation DOUBLE,"
	  "PRIMARY KEY(runID, dimiunID, eventID))");
#ifndef OUT_TO_SCREEN
  server->Exec(query);
#else
  std::cout << __FUNCTION__ << ": " << query << std::endl;
#endif
}

void MySQLSvc::writeTrackingRes(SRecEvent* recEvent, TClonesArray* tracklets)
{
  //Fill Track table/TrackHit table
  mom_vertex.clear(); 
  pos_vertex.clear();
  int nTracks_local = recEvent->getNTracks();
  for(int i = 0; i < nTracks_local; ++i)
    {
      int trackID = nTracks + i;
      writeTrackTable(trackID, &recEvent->getTrack(i));
      writeTrackHitTable(trackID, (Tracklet*)tracklets->At(i));
    }
  
  //Fill dimuon table
  std::vector<int> idx_plus = recEvent->getChargedTrackIDs(+1);
  std::vector<int> idx_minus = recEvent->getChargedTrackIDs(-1);
  for(unsigned int i = 0; i < idx_plus.size(); ++i)
    {
      for(unsigned int j = 0; j < idx_minus.size(); ++j)
	{
	  writeDimuonTable(nDimuons++, idx_plus[i], idx_minus[j]);
	}
    }

  nTracks += nTracks_local;
}

void MySQLSvc::writeTrackTable(int trackID, SRecTrack* recTrack)
{      
  double px0, py0, pz0, x0, y0, z0; 
  double px1, py1, pz1, x1, y1, z1;
  double px3, py3, pz3, x3, y3, z3;

  int charge;
  int numHits;
  double chisq;

  //Track related
  charge = recTrack->getCharge();
  numHits = recTrack->getNHits();
  chisq = recTrack->getChisq();

  //Vertex point
  x0 = recTrack->getVtxPar(0);
  y0 = recTrack->getVtxPar(1);
  z0 = recTrack->getVtxPar(2);
  recTrack->getMomentumVertex(px0, py0, pz0);

  //At station 1
  z1 = 600.;
  recTrack->getExpPositionFast(z1, x1, y1);
  recTrack->getExpMomentumFast(z1, px1, py1, pz1);

  //At station 3
  z3 = 1900.;
  recTrack->getExpPositionFast(z3, x3, y3);
  recTrack->getExpMomentumFast(z3, px3, py3, pz3);

  //Put the paramters in temporary containers for dimuon combination
  mom_vertex.push_back(TLorentzVector(px0, py0, pz0, sqrt(px0*px0 + py0*py0 + pz0*pz0 + 0.10566*0.10566)));
  pos_vertex.push_back(TVector3(x0, y0, z0));

  //Database output
  sprintf(query, "INSERT INTO kTrack(trackID,runID,spillID,eventID,charge,numHits,chisq,x0,y0,z0,px0,py0,pz0," 
	  "x1,y1,z1,px1,py1,pz1,x3,y3,z3,px3,py3,pz3 VALUE(%d,%d,%d,%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,"
	  "%f,%f,%f,%f,%f,%f,%f,%f,)", trackID, runID, spillID, eventID_last, charge, numHits, chisq, x0, y0, 
	  z0, px0, py0, pz0, x1, y1, z1, px0, py0, pz0, x3, y3, z3, px3, py3, pz3);
#ifndef OUT_TO_SCREEN
  server->Exec(query);
#else
  std::cout << __FUNCTION__ << ": " << query << std::endl;
#endif
}

void MySQLSvc::writeTrackHitTable(int trackID, Tracklet* tracklet)
{
  for(std::list<SignedHit>::iterator iter = tracklet->hits.begin(); iter != tracklet->hits.end(); ++iter)
    {
      if(iter->hit.index < 0) continue;

      sprintf(query, "INSERT INTO kTrackHit(runID,eventID,trackID,hitID,driftSign,residual) VALUES(%d,%d,%d,%d,%d,%f)",
	      runID, eventID_last, trackID, iter->hit.index, iter->sign, tracklet->residual[iter->hit.detectorID-1]);
#ifndef OUT_TO_SCREEN
      server->Exec(query);
#else
      std::cout << __FUNCTION__ << ": " << query << std::endl;
#endif
    }
}

void MySQLSvc::writeDimuonTable(int dimuonID, int idx_positive, int idx_negative)
{
  double mp = 0.938;
  double ebeam = 120.;

  TLorentzVector p_beam(0., 0., sqrt(ebeam*ebeam - mp*mp), ebeam);
  TLorentzVector p_target(0., 0., 0., mp);
  TLorentzVector p_cms = p_beam + p_target;
  TVector3 bv_cms = p_cms.BoostVector();
  double s = p_cms.M2();

  TLorentzVector p_sum = mom_vertex[idx_positive] + mom_vertex[idx_negative];
  TVector3 v_sum = pos_vertex[idx_positive] + pos_vertex[idx_negative];

  double mass, xF, x1, x2, pT, x0, y0, z0, px0, py0, pz0, dz;

  mass = p_sum.M();
  pT = p_sum.Perp();
  px0 = p_sum.Px();
  py0 = p_sum.Py();
  pz0 = p_sum.Pz();
  dz = (pos_vertex[idx_positive] - pos_vertex[idx_negative]).Z();

  p_sum.Boost(-bv_cms);
  xF = 2*p_sum.Pz()/TMath::Sqrt(s);
  double tau = p_sum.M2()/s;
  double y = 0.5*std::log((p_sum.E() + p_sum.Pz())/(p_sum.E() - p_sum.Pz()));

  x1 = TMath::Sqrt(tau)*TMath::Exp(y);
  x2 = TMath::Sqrt(tau)*TMath::Exp(-y);

  x0 = v_sum.X();
  y0 = v_sum.Y();
  z0 = v_sum.Z();

  sprintf(query, "INSERT INTO kDimuon(dimuonID,runID,spillID,eventID,posTrackID,negTrackID,dx,dy,dz,dpx,"
	  "dpy,dpz,mass,xF,xB,xT,trackSeparation) VALUES(%d,%d,%d,%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f)", 
	  dimuonID, runID, spillID, eventID_last, idx_positive+nTracks, idx_negative+nTracks, x0, y0, z0, 
	  px0, py0, pz0, mass, xF, x1, x2, dz);
#ifndef OUT_TO_SCREEN
  server->Exec(query);
#else
  std::cout << __FUNCTION__ << ": " << query << std::endl;
#endif
}

int MySQLSvc::getNEventsFast()
{
  if(nEvents < 1) nEvents = getNEvents();
  return nEvents;
}

int MySQLSvc::makeQuery()
{
  //std::cout << query << std::endl;
  if(server == NULL) return 0;

  if(res != NULL) delete res;
  res = server->Query(query);

  if(res != NULL) return res->GetRowCount();
  return 0;
}

bool MySQLSvc::nextEntry()
{
  if(res == NULL) return false;

  if(row != NULL) delete row;
  row = res->Next();

  if(row != NULL) return true;
  return false;
}
