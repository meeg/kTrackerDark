#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>
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

TSQLResult *makeQuery(TSQLServer *con, string query)
{
  //cout << query << endl;
  return con->Query(query.c_str());
}

int main(int argc, char **argv)
{
  cout << "Exporting Run: " << argv[1] << " to ROOT file: " << argv[2] << endl;

  GeomSvc *p_geomSvc = GeomSvc::instance();
  p_geomSvc->init(GEOMETRY_VERSION);

  Double_t mass, xF, x1, x2;

  Int_t nHits[2];
  Int_t charge[2];
  Double_t weight;
  Double_t px0[2], py0[2], pz0[2], p0[2];
  Double_t px1[2], py1[2], pz1[2], p1[2];
  Double_t px2[2], py2[2], pz2[2], p2[2];
  Double_t pxf[2], pyf[2], pzf[2], pf[2];
  Double_t pxm[2], pym[2], pzm[2], pm[2];
  Double_t x0, y0, z0;

  Int_t nSeeds, nSeedsX, nSeedsY;
  Double_t ax[50], ay[50], bx[50], by[50];
  Int_t xIndex[50], yIndex[50];

  SRawEvent *rawEvent = new SRawEvent();

  TFile *saveFile = new TFile(argv[2], "recreate");
  TTree *saveTree = new TTree("save", "save");

  saveTree->Branch("rawEvent", &rawEvent, 256000, 99);

  saveTree->Branch("nHits", nHits, "nHits[2]/I");
  saveTree->Branch("mcharge", charge, "mcharge[2]/I");
  saveTree->Branch("weight", &weight, "weight/D");

  saveTree->Branch("mmass", &mass, "mmass/D");
  saveTree->Branch("mxF", &xF, "mxF/D");
  saveTree->Branch("mx1", &x1, "mx1/D");
  saveTree->Branch("mx2", &x2, "mx2/D");

  saveTree->Branch("x0", &x0, "x0/D");
  saveTree->Branch("y0", &y0, "y0/D");
  saveTree->Branch("z0", &z0, "z0/D");
  
  saveTree->Branch("px0", px0, "px0[2]/D");
  saveTree->Branch("py0", py0, "py0[2]/D");
  saveTree->Branch("pz0", pz0, "pz0[2]/D");
  saveTree->Branch("p0", p0, "p0[2]/D");
 
  saveTree->Branch("px1", px1, "px1[2]/D");
  saveTree->Branch("py1", py1, "py1[2]/D");
  saveTree->Branch("pz1", pz1, "pz1[2]/D");
  saveTree->Branch("p1", p1, "p1[2]/D");
 
  saveTree->Branch("px2", px2, "px2[2]/D");
  saveTree->Branch("py2", py2, "py2[2]/D");
  saveTree->Branch("pz2", pz2, "pz2[2]/D");
  saveTree->Branch("p2", p2, "p2[2]/D");

  saveTree->Branch("pxf", pxf, "pxf[2]/D");
  saveTree->Branch("pyf", pyf, "pyf[2]/D");
  saveTree->Branch("pzf", pzf, "pzf[2]/D");
  saveTree->Branch("pf", pf, "pf[2]/D");

  saveTree->Branch("nSeeds", &nSeeds, "nSeeds/I");
  saveTree->Branch("nSeedsX", &nSeedsX, "nSeedsX/I");  
  saveTree->Branch("nSeedsY", &nSeedsY, "nSeedsY/I");
    
  saveTree->Branch("ax", ax, "ax[nSeedsX]/D");
  saveTree->Branch("bx", bx, "bx[nSeedsX]/D");
  saveTree->Branch("ay", ay, "ay[nSeedsY]/D");
  saveTree->Branch("by", by, "by[nSeedsY]/D");
  saveTree->Branch("xIndex", xIndex, "xIndex[nSeeds]/I");
  saveTree->Branch("yIndex", yIndex, "yIndex[nSeeds]/I");
 
  ///Connect to the SQL databases
  TSQLServer *con = TSQLServer::Connect("mysql://seaquel.physics.illinois.edu", "seaguest","qqbar2mu+mu-");
  const char *buf1 = "SELECT runID,spillID,eventID,mTrackID1,mTrackID2,sigWeight,mass,xF,xB,xT FROM %s.mDimuon WHERE acceptHodoAll=1 AND"
    " acceptDriftAll=1 ORDER BY eventID LIMIT %d";
  const char *buf2 = "SELECT x0,y0,z0,px0,py0,pz0,pxf,pyf,pzf,particleID FROM %s.mTrack WHERE (mTrackID=%d OR mTrackID=%d) ORDER BY particleID"; 
  const char *buf3 = "SELECT a.hpx,a.hpy,a.hpz FROM %s.mGeantHit AS a,%s.mHit AS b WHERE a.geantHitID=b.geantHitID AND a.mTrackID=%d AND b.detectorName"
    " LIKE 'D%%' ORDER BY b.zAtDigiPlane";
  const char *buf4 = "SELECT xAtDigiPlane,yAtDigiPlane,zAtDigiPlane FROM %s.mHit WHERE mTrackID=%d AND detectorName LIKE 'P%%' ORDER BY zAtDigiPlane";
  const char *buf5 = "SELECT hitID,elementID,driftTime,driftDistance,detectorName FROM %s.mHit WHERE runID=%d AND spillID=%d"
    " AND eventID=%d AND (detectorName LIKE 'D%%' OR detectorName LIKE 'H%%' OR detectorName LIKE 'P%%') ORDER BY mTrackID,detectorName";

  char query[500];
  sprintf(query, buf1, argv[1], atoi(argv[3])); 
  TSQLResult *res_dimuon = makeQuery(con, query);

  UInt_t nEntries = res_dimuon->GetRowCount();
  cout << "Totally " << nEntries << " dimuon pairs in this run" << endl;
  for(UInt_t i = 0; i < nEntries; i++)
    {
      if(i % 1 == 0)
	{
	  cout << "Converting dimuon pair " << i << endl;
	}

      TSQLRow *row_dimuon = res_dimuon->Next();

      Int_t runID, spillID, eventID;
      Int_t trackIDs[2];
      runID = atoi(row_dimuon->GetField(0));
      spillID = atoi(row_dimuon->GetField(1));
      eventID = atoi(row_dimuon->GetField(2));
      trackIDs[0] = atoi(row_dimuon->GetField(3));
      trackIDs[1] = atoi(row_dimuon->GetField(4));
      weight = atof(row_dimuon->GetField(5));
      mass = atof(row_dimuon->GetField(6));
      xF = atof(row_dimuon->GetField(7));
      x1 = atof(row_dimuon->GetField(8));
      x2 = atof(row_dimuon->GetField(9));

      sprintf(query, buf2, argv[1], trackIDs[0], trackIDs[1]);
      TSQLResult *res_track = makeQuery(con, query);

      Int_t nTracks = res_track->GetRowCount();
      if(nTracks != 2)
	{
	  cout << "Error: more than 2 tracks in one event! " << endl;
	  cout << eventID << endl; 
	}

      nSeeds = 0;
      nSeedsX = 0;
      nSeedsY = 0;	  
      for(Int_t j = 0; j < nTracks; j++)
	{
	  TSQLRow *row_track = res_track->Next();
	
    	  x0 = atof(row_track->GetField(0));
          y0 = atof(row_track->GetField(1));
          z0 = atof(row_track->GetField(2));
          px0[j] = atof(row_track->GetField(3));
          py0[j] = atof(row_track->GetField(4));
          pz0[j] = atof(row_track->GetField(5));
     	  pxf[j] = atof(row_track->GetField(6));
          pyf[j] = atof(row_track->GetField(7));
          pzf[j] = atof(row_track->GetField(8));
	  charge[j] = atoi(row_track->GetField(9)) > 0 ? 1 : -1;

	  p0[j] = sqrt(px0[j]*px0[j] + py0[j]*py0[j] + pz0[j]*pz0[j]);
	  pf[j] = sqrt(pxf[j]*pxf[j] + pyf[j]*pyf[j] + pzf[j]*pzf[j]);

	  sprintf(query, buf3, argv[1], argv[1], trackIDs[j]);
	  TSQLResult *res_mom = makeQuery(con, query);

	  nHits[j] = res_mom->GetRowCount();
	  for(Int_t k = 0; k < nHits[j]; k++)
	    {
	      TSQLRow *row_mom = res_mom->Next();
	      if(k == 0)
		{
		  px1[j] = atof(row_mom->GetField(0));
		  py1[j] = atof(row_mom->GetField(1));
		  pz1[j] = atof(row_mom->GetField(2));

		  p1[j] = sqrt(px1[j]*px1[j] + py1[j]*py1[j] + pz1[j]*pz1[j]);
		}
	      else if(k == nHits[j]-1)
		{
		  px2[j] = atof(row_mom->GetField(0));
		  py2[j] = atof(row_mom->GetField(1));
		  pz2[j] = atof(row_mom->GetField(2));

		  p2[j] = sqrt(px2[j]*px2[j] + py2[j]*py2[j] + pz2[j]*pz2[j]);
		}

	      delete row_mom;
	    }

	  sprintf(query, buf4, argv[1], trackIDs[j]);
	  TSQLResult *res_seed = makeQuery(con, query);

	  Int_t nSeedHits = res_seed->GetRowCount();
	  if(nSeedHits < 2) break;

	  Double_t x_start, y_start, z_start;
	  Double_t x_end, y_end, z_end;
	  for(Int_t k = 0; k < nSeedHits; k++)
	    {
	      TSQLRow *row_seed = res_seed->Next();

	      if(k == 0)
		{
	     	  x_start = atof(row_seed->GetField(0));
		  y_start = atof(row_seed->GetField(1));
		  z_start = atof(row_seed->GetField(2));
		}
	      else if(k == nSeedHits - 1)
		{
		  x_end = atof(row_seed->GetField(0));
		  y_end = atof(row_seed->GetField(1));
		  z_end = atof(row_seed->GetField(2));
		}

	      delete row_seed;
	    }

	  ax[nSeedsX] = (x_end - x_start)/(z_end - z_start);
	  ay[nSeedsY] = (y_end - y_start)/(z_end - z_start);
	  bx[nSeedsX] = x_end - ax[nSeedsX]*z_end;
	  by[nSeedsY] = y_end - ay[nSeedsY]*z_end;

	  xIndex[nSeeds] = nSeedsX;
	  yIndex[nSeeds] = nSeedsY;

	  ++nSeeds; ++nSeedsX; ++nSeedsY;
	  
	  delete row_track;

	  delete res_seed;
	  delete res_mom; 
	}

      sprintf(query, buf5, argv[1], runID, spillID, eventID);
      TSQLResult *res_hit = makeQuery(con, query);

      rawEvent->setEventInfo(runID, spillID, eventID);
      Int_t nHitsAll = res_hit->GetRowCount();
      for(Int_t j = 0; j < nHitsAll; j++)
	{
	  TSQLRow *row_hit = res_hit->Next();

	  string detectorName(row_hit->GetField(4));
	  Int_t elementID = atoi(row_hit->GetField(1));
	  p_geomSvc->toLocalDetectorName(detectorName, elementID);

	  Hit h;
    	  h.index = j;//atoi(row_hit->GetField(0));
  	  h.detectorID = p_geomSvc->getDetectorID(detectorName);
	  h.elementID = elementID;
	
	  double pos, err;
	  p_geomSvc->getMeasurement(h.detectorID, h.elementID, pos, err);
	  h.pos = pos;

	  if(row_hit->GetField(2) != NULL)
	    {
	      h.driftTime = atof(row_hit->GetField(2));
	    }
	  else
	    {
	      h.driftTime = 0.;
	    }
          h.tdcTime = h.driftTime;

	  if(row_hit->GetField(3) != NULL)
	    {
	      h.driftDistance = fabs(atof(row_hit->GetField(3)));
	    }
	  else
	    {
	      h.driftDistance = 0.;
	    }

	  h.inTime = 1;
	  h.hodoMask = 1;
	  rawEvent->insertHit(h);

	  delete row_hit;
	}

      saveTree->Fill();
      rawEvent->clear();

      delete row_dimuon;
      delete res_track;
      delete res_hit;
    }

  delete res_dimuon;
  delete con;

  saveFile->cd();
  saveTree->Write();
  saveFile->Close();
  
  return 1;
}
