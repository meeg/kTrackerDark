/*
SRawEvent.h

Definition of the class SRawEvent, which essentially works as a 
container of the raw hits information. It also provides serveral 
query utility to retrieve the hit list from a specific plane, etc.

Author: Kun Liu, liuk@fnal.gov
Created: 07-02-2012
*/

#ifndef _SRAWEVENT_H
#define _SRAWEVENT_H

#include "MODE_SWITCH.h"

#include <iostream>
#include <vector>
#include <list>
#include <string>

#include <TObject.h>
#include <TROOT.h>
#include <TVector3.h>

#include "GeomSvc.h"

///Definition of hit structure
class Hit: public TObject
{
public:
  Int_t index;         //unique index for identification
  Int_t detectorID;    //assigned for each detector plane
  Int_t elementID;     
  Double_t tdcTime;    //raw TDC time
  Double_t driftTime;
  Double_t driftDistance; 
  Double_t pos;        //actual measurement in either X, Y, U or V direction

  Int_t inTime;        //In-time flag
  Int_t hodoMask;      //Hodo-mask flag

  //Sign of this hit
  int getSign() { return driftDistance > 0 ? 1 : -1; }

  //overiden comparison operator for track seeding 
  bool operator<(const Hit& elem) const;
  bool operator>(const Hit& elem) const { return detectorID > elem.detectorID; }
  bool operator==(const Hit& elem) const { return detectorID == elem.detectorID; }

  //comparison function for after pulse removal
  static bool sameChannel(const Hit& elem1, const Hit& elem2);

  //Debugging output
  void print() { std::cout << index << " : " << detectorID << " : " << elementID << " : " << pos << " : " << driftDistance << " : " << inTime << " : " << hodoMask << std::endl; }

  ClassDef(Hit, 2)
};

class SRawEvent: public TObject
{
public:
  SRawEvent();
  ~SRawEvent();

  ///Mix events for MC study
  void mixEvent(SRawEvent *event, int nBkgHits = -1);

  ///Gets
  std::list<Int_t> getHitsIndexInDetector(Int_t detectorID);
  std::list<Int_t> getHitsIndexInDetector(Int_t detectorID, Double_t x_exp, Double_t win);
  std::list<Int_t> getHitsIndexInSuperDetector(Int_t detectorID);
  std::list<Int_t> getHitsIndexInDetectors(std::vector<Int_t>& detectorIDs);
  std::list<Int_t> getAdjacentHitsIndex(Hit& _hit);

  Int_t getNHitsAll() { return fNHits[0]; }
  Int_t getNChamberHitsAll();
  Int_t getNHodoHitsAll();
  Int_t getNPropHitsAll();
  Int_t getNHitsInD1();
  Int_t getNHitsInD2();
  Int_t getNHitsInD3p();
  Int_t getNHitsInD3m();

  Int_t getNHitsInDetector(Int_t detectorID) { return fNHits[detectorID]; }
  Int_t getNHitsInSuperDetector(Int_t detectorID) { return fNHits[2*detectorID-1] + fNHits[2*detectorID]; }
  Int_t getNHitsInDetectors(std::vector<Int_t>& detectorIDs);
  
  std::vector<Hit> getAllHits() { return fAllHits; }
  Hit getHit(Int_t index) { return fAllHits[index]; } 

  Int_t getRunID() { return fRunID; }
  Int_t getEventID() { return fEventID; }
  Int_t getSpillID() { return fSpillID; }

  ///Sets
  void setEventInfo(Int_t runID, Int_t spillID, Int_t eventID);

  ///Insert a new hit
  void insertHit(Hit h);

  ///Find a hit
  int findHit(int detectorID, int elementID);

  ///Manipulation of hit list
  void reIndex(std::string option = "");

  ///Type of pair with two adjacent wiree
  typedef std::pair<int, int> hit_pair;
  std::list<SRawEvent::hit_pair> getHitPairsInSuperDetector(Int_t detectorID);
  std::list<SRawEvent::hit_pair> getPartialHitPairsInSuperDetector(Int_t detectorID);  
  std::list<SRawEvent::hit_pair> getHitPairsInSuperDetector(Int_t detectorID, Double_t x_exp, Double_t wind);
  std::list<SRawEvent::hit_pair> getPartialHitPairsInSuperDetector(Int_t detectorID, Double_t x_exp, Double_t wind);  
  
  ///Clear the internal event structure
  void clear();

  ///Print for debugging purposes
  void print(); 

private:
  Int_t fRunID;
  Int_t fEventID;
  Int_t fSpillID;

  ///Hits of this event
  Int_t fNHits[nChamberPlanes+nHodoPlanes+nPropPlanes+1];  //0 for all hits, 1, 2, ..., 24 for number of hits in plane 1, 2, ..., 24
  std::vector<Hit> fAllHits;

  ClassDef(SRawEvent, 2)
};

class SRawMCEvent: public SRawEvent
{
public:
  //sigWeight
  Double_t weight;

  //Dimuon info
  Double_t mass;
  Double_t xF;
  Double_t pT;
  Double_t x1;
  Double_t x2;
  TVector3 vtx;
 
  //Track info, 0 for mu+, 1 for mu-
  Int_t nHits[2];
  TVector3 p_vertex[2];
  TVector3 p_station1[2];
  TVector3 v_station1[2];
  TVector3 p_station2[2];
  TVector3 v_station2[2];
  TVector3 p_station3[2];
  TVector3 v_station3[2];
  TVector3 p_station4[2];
  TVector3 v_station4[2];

  ClassDef(SRawMCEvent, 1) 
};

#endif
