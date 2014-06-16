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

#define triggerBit(n) (1 << (n))

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
  bool operator==(const Hit& elem) const;

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
  void mixEvent(SRawEvent* event, int nBkgHits = -1);

  ///Gets
  std::list<Int_t> getHitsIndexInDetector(Int_t detectorID);
  std::list<Int_t> getHitsIndexInDetector(Int_t detectorID, Double_t x_exp, Double_t win);
  std::list<Int_t> getHitsIndexInSuperDetector(Int_t detectorID);
  std::list<Int_t> getHitsIndexInDetectors(std::vector<Int_t>& detectorIDs);
  std::list<Int_t> getAdjacentHitsIndex(Hit& _hit);

  Int_t getNHitsAll() { return fNHits[0]; }
  Int_t getNTriggerHits() { return fTriggerHits.size(); }
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
  std::vector<Hit> getTriggerHits() { return fTriggerHits; }
  Hit getTriggerHit(Int_t index) { return fTriggerHits[index]; } 

  Hit getHit(Int_t index) { return fAllHits[index]; } 
  Hit getHit(Int_t detectorID, Int_t elementID); 
  void setHit(Int_t index, Hit hit) { fAllHits[index] = hit; }
  void setTriggerHit(Int_t index, Hit hit) { fTriggerHits[index] = hit; }
  void setHitFlag(Int_t index, Int_t flag) { if(index < 0) return; fAllHits[index].inTime = flag; }
  void setHitFlag(Int_t detectorID, Int_t elementID, Int_t flag) { setHitFlag(findHit(detectorID, elementID), flag); }

  Int_t getRunID() { return fRunID; }
  Int_t getEventID() { return fEventID; }
  Int_t getSpillID() { return fSpillID; }

  ///Sets
  void setEventInfo(Int_t runID, Int_t spillID, Int_t eventID);

  ///Insert a new hit
  void insertHit(Hit h);
  void insertTriggerHit(Hit h) { fTriggerHits.push_back(h); }
  
  ///Find a hit -- binary search since hit list is sorted
  Int_t findHit(Int_t detectorID, Int_t elementID);

  ///Manipulation/reduction of hit list
  void reIndex(std::string option = "");
  void deClusterize(std::list<Hit>& hits);
  void processCluster(std::list<Hit>& hits, std::vector<std::list<Hit>::iterator>& cluster);

  ///Type of pair with two adjacent wiree
  typedef std::pair<Int_t, Int_t> hit_pair;
  std::list<SRawEvent::hit_pair> getHitPairsInSuperDetector(Int_t detectorID);
  std::list<SRawEvent::hit_pair> getPartialHitPairsInSuperDetector(Int_t detectorID);  
  std::list<SRawEvent::hit_pair> getHitPairsInSuperDetector(Int_t detectorID, Double_t x_exp, Double_t wind);
  std::list<SRawEvent::hit_pair> getPartialHitPairsInSuperDetector(Int_t detectorID, Double_t x_exp, Double_t wind);  
  
  ///Set/get the trigger types
  Int_t getTriggerBits() { return fTriggerBits; }
  void setTriggerBits(Int_t triggers[]);
  void setTriggerBits(Int_t triggers) { fTriggerBits = triggers; }
  bool isTriggeredBy(Int_t trigger) { return (fTriggerBits & trigger) != 0; }

  //Set/get offline trigger emulation results
  bool isEmuTriggered() { return fTriggerEmu > 0; }
  Int_t getNRoadsPos() { return fNRoads[0] + fNRoads[1]; } 
  Int_t getNRoadsNeg() { return fNRoads[2] + fNRoads[3]; }
  Int_t getNRoadsPosTop() { return fNRoads[0]; } 
  Int_t getNRoadsPosBot() { return fNRoads[1]; } 
  Int_t getNRoadsNegTop() { return fNRoads[2]; } 
  Int_t getNRoadsNegBot() { return fNRoads[3]; } 
  Int_t* getNRoads() { return fNRoads; }
  void setTriggerEmu(bool flag) { fTriggerEmu = flag ? 1 : -1; }
  void setNRoads(Int_t nRoads[]) { for(Int_t i = 0; i < 4; ++i) fNRoads[i] = nRoads[i]; }

  //Set/get the target position
  Int_t getTargetPos() { return fTargetPos; }
  void setTargetPos(Int_t targetPos) { fTargetPos = targetPos; }

  //Set/get the beam info
  Int_t getTurnID() { return fTurnID; }
  Int_t getRFID() { return fRFID; }
  Int_t getIntensity() { return fIntensity[16]; }
  Int_t getIntensity(Int_t i) { return fIntensity[i+16]; }
  Int_t getIntensitySumBefore(Int_t n = 16) { Int_t sum = 0; for(Int_t i = n; i < 16; ++i) sum += fIntensity[i]; return sum; } 
  Int_t getIntensitySumAfter(Int_t n = 16) { Int_t sum = 0; for(Int_t i = 16; i < n+16; ++i) sum += fIntensity[i]; return sum; } 
  Int_t* getIntensityAll() { return fIntensity; }

  void setTurnID(Int_t turnID) { fTurnID = turnID; }
  void setRFID(Int_t rfID) { fRFID = rfID; }
  void setIntensity(const Int_t intensity[]) { for(Int_t i = 0; i < 33; ++i) fIntensity[i] = intensity[i]; }
  void setIntensity(Int_t i, Int_t val) { fIntensity[i] = val; }
  void setIntensity(Int_t val) { fIntensity[16] = val; }

  //Set the event info from another event
  void setEventInfo(SRawEvent* event);

  ///Clear the internal event structure
  void clear();

  ///Print for debugging purposes
  void print(); 

public:
  //Trigger type
  enum TriggerType 
    {
      NIM1 = triggerBit(1),
      NIM2 = triggerBit(2),
      NIM3 = triggerBit(3),
      NIM4 = triggerBit(4),
      NIM5 = triggerBit(5),
      MATRIX1 = triggerBit(6),
      MATRIX2 = triggerBit(7),
      MATRIX3 = triggerBit(8),
      MATRIX4 = triggerBit(9),
      MATRIX5 = triggerBit(10)
    };

private:
  //RunID, spillID, eventID
  Int_t fRunID;
  Int_t fEventID;
  Int_t fSpillID;

  //Trigger bit
  Int_t fTriggerBits;

  //Target pos
  Int_t fTargetPos;

  //Beam intensity information
  Int_t fTurnID;
  Int_t fRFID;
  Int_t fIntensity[33];   //16 before, one onset, and 16 after

  //Offline trigger simulation res
  Int_t fTriggerEmu;
  Int_t fNRoads[4];       //0, positive top; 1, positive bottom; 2, negative top; 3, negative bottom

  ///Hits of this event
  Int_t fNHits[nChamberPlanes+nHodoPlanes+nPropPlanes+1];  //0 for all hits, 1, 2, ..., 24 for number of hits in plane 1, 2, ..., 24
  std::vector<Hit> fAllHits;
  std::vector<Hit> fTriggerHits;

  ClassDef(SRawEvent, 7)
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
  Double_t costh;
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

  TVector3 p_stationH1[2];
  TVector3 v_stationH1[2];
  TVector3 p_stationH2[2];
  TVector3 v_stationH2[2];
  TVector3 p_stationH3[2];
  TVector3 v_stationH3[2];
  TVector3 p_stationH4[2];
  TVector3 v_stationH4[2];

  ClassDef(SRawMCEvent, 3) 
};

#endif
