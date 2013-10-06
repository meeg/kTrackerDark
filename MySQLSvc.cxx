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

bool MySQLSvc::isNewEvtAvailable()
{
  sprintf(query, "SELECT eventID,runID,spillID FROM %s.Event ORDER BY eventID DESC LIMIT 2", dataSchema.c_str());
  if(!makeQuery()) return false;
  if(res->GetRowCount() < 2) return false;
  
  nextEntry(); nextEntry();
  int eventID_curr = atoi(row->GetField(0));
  if(eventID_curr == eventID_last) return false;
    
  eventID_last = eventID_curr;
  runID = atoi(row->GetField(1));
  spillID = atoi(row->GetField(2));

  return true;
}

bool MySQLSvc::isRunStopped()
{
  sprintf(query, "SELECT productionEnd from %s.production WHERE runID=%d", logSchema.c_str(), runID);
  if(!makeQuery()) return false;
  if(res->GetRowCount() != 1) return false;

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
  return getEvent(rawEvent, eventID_last);
}

bool MySQLSvc::getRandomEvt(SRawEvent* rawEvent)
{
  int eventID = int(rndm.Rndm()*getNEventsFast());
  return getEvent(rawEvent, eventID);
}

bool MySQLSvc::getNextEvent(SRawEvent* rawEvent)
{
  return getEvent(rawEvent, eventID_last+1);
}

bool MySQLSvc::getEvent(SRawEvent* rawEvent, int eventID)
{
  rawEvent->clear();
  
  sprintf(query, "SELECT hitID,elementID,tdcTime,driftTime,driftDistance,detectorName,inTime,masked FROM %s.Hit WHERE (detectorName LIKE 'D%%' OR detectorName LIKE 'H%%' OR detectorName LIKE 'P%%') AND inTime=1 AND eventID=%d", dataSchema.c_str(), eventID);
  if(!makeQuery()) return false;

  int nHits = res->GetRowCount();
  if(nHits < 1) return false;
  //Log(nHits);

  rawEvent->setEventInfo(runID, spillID, eventID);
  for(int i = 0; i < nHits; ++i)
    {
      nextEntry();

      std::string detectorName(row->GetField(5));
      int elementID = atoi(row->GetField(1));
      p_geomSvc->toLocalDetectorName(detectorName, elementID);
       
      Hit h;
      h.index = atoi(row->GetField(0));
      h.detectorID = p_geomSvc->getDetectorID(detectorName);
      h.elementID = elementID;
      h.tdcTime = atof(row->GetField(2));
      h.inTime = atoi(row->GetField(6));
      h.pos = p_geomSvc->getMeasurement(h.detectorID, h.elementID);
    
      if(row->GetField(3) == NULL)
	{
	  h.driftTime = 0.;
	}	  
      else
	{
	  h.driftTime = atof(row->GetField(3));
	}

      if(row->GetField(4) == NULL)
	{
	  h.driftDistance = 0.;
	}
      else
	{
	  h.driftDistance = atof(row->GetField(4));
	}

      if(row->GetField(7) == NULL)
	{
	  h.hodoMask = 1;
	}
      else
	{
	  h.hodoMask = atoi(row->GetField(7));
	}

      rawEvent->insertHit(h);
    }

  eventID_last = eventID;
  return true;
}

int MySQLSvc::getNEvents()
{
  sprintf(query, "SELECT COUNT(*) FROM %s.Event", dataSchema.c_str());
  if(!makeQuery()) return 0;

  nextEntry();
  int nTotal = atoi(row->GetField(0));
  return nTotal;
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
  double x0 = recTrack->getVtxPar(0);
  double y0 = recTrack->getVtxPar(1);
  double z0 = recTrack->getVtxPar(2);
  double px0, py0, pz0;
  recTrack->getMomentumVertex(px0, py0, pz0);
  int charge = recTrack->getCharge();

  //Put the paramters in temporary containers for dimuon combination
  mom_vertex.push_back(TLorentzVector(px0, py0, pz0, sqrt(px0*px0 + py0*py0 + pz0*pz0 + 0.10566*0.10566)));
  pos_vertex.push_back(TVector3(x0, y0, z0));

  sprintf(query, "INSERT INTO kTrack(trackID,runID,spillID,eventID,x0,y0,z0,px0,py0,pz0,charge)" 
	  " VALUES(%d,%d,%d,%d,%f,%f,%f,%f,%f,%f,%d)", trackID, runID, spillID, eventID_last, x0,
	  y0, z0, px0, py0, pz0, charge);
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

      sprintf(query, "INSERT INTO kTrackHit(runID,trackID,hitID) VALUES(%d,%d,%d)", runID, trackID, iter->hit.index);
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

  double mass, xF, x1, x2, pT, x0, y0, z0, px0, py0, pz0;

  mass = p_sum.M();
  pT = p_sum.Perp();
  px0 = p_sum.Px();
  py0 = p_sum.Py();
  pz0 = p_sum.Pz();

  p_sum.Boost(-bv_cms);
  xF = 2*p_sum.Pz()/TMath::Sqrt(s);
  double tau = p_sum.M2()/s;
  double y = 0.5*std::log((p_sum.E() + p_sum.Pz())/(p_sum.E() - p_sum.Pz()));

  x1 = TMath::Sqrt(tau)*TMath::Exp(y);
  x2 = TMath::Sqrt(tau)*TMath::Exp(-y);

  x0 = v_sum.X();
  y0 = v_sum.Y();
  z0 = v_sum.Z();

  sprintf(query, "INSERT INTO kDimuon(dimuonID,runID,spillID,eventID,trackID1,trackID2,mass,xF,xB,xT"
	  "dx,dy,dz,dpx,dpy,dpz,dpT,phi_gam,phi_mu,theta_mu) VALUES(%d,%d,%d,%d,%d,%d,%f,%f,%f,%f,%f,"
	  "%f,%f,%f,%f,%f,%f,%f,%f,%f)", dimuonID, runID, spillID, eventID_last, idx_positive+nTracks,
	  idx_negative+nTracks, mass, xF, x1, x2, x0, y0, z0, px0, py0, pz0, pT, 0., 0., 0.);
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

bool MySQLSvc::makeQuery()
{
  //std::cout << query << std::endl;
  if(server == NULL) return false;

  if(res != NULL) delete res;
  res = server->Query(query);

  if(res != NULL) return true;
  return false;
}

bool MySQLSvc::nextEntry()
{
  if(res == NULL) return false;

  if(row != NULL) delete row;
  row = res->Next();

  if(row != NULL) return true;
  return false;
}
