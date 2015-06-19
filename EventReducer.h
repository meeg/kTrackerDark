/*
EventReducer.h

This class is intended to handle all the hit list manipulation/reduction,
The reduction methods could initialized once for all, and then it will update
the hit list stored in SRawEvent. It is declared as a friend class of SRawEvent.

Author: Kun Liu, liuk@fnal.gov
Created: 06-17-2015
*/

#ifndef _EVENTREDUCER_H
#define _EVENTREDUCER_H

#include "MODE_SWITCH.h"

#include <list>
#include <TString.h>
#include <TRandom.h>

#include "GeomSvc.h"
#include "SRawEvent.h"
#include "TriggerAnalyzer.h"

class EventReducer
{
public:
  EventReducer(TString options);
  ~EventReducer();

  //main external call
  int reduceEvent(SRawEvent* rawEvent);

  //sagitta ratio reducer
  void sagittaReducer();

  //hough transform reducer
  void houghReducer();

  //hit cluster remover
  void deClusterize();
  void processCluster(std::vector<std::list<Hit>::iterator>& cluster);

private:
  //pointer to geometry service, inited outside
  GeomSvc* p_geomSvc;

  //pointer to trigger analyzer, inited inside
  TriggerAnalyzer* p_triggerAna;

  //Random number
  TRandom rndm;

  //temporary container for the hit list
  std::list<Hit> hitlist;

  //flags of the hit manipulation method
  bool afterhit;            //after pulse removal
  bool hodomask;            //hodoscope masking
  bool outoftime;           //out of time hit removal
  bool decluster;           //remove hit clusters in chamber
  bool mergehodo;           //merge trigger hit with hit
  bool triggermask;         //use active trigger road for track masking
  bool sagitta;             //remove the hits which cannot form a sagitta triplet
  bool hough;               //remove the hits which cannot form a peak in hough space, will be implemented later
  bool externalpar;         //re-apply the alignment and calibration parameters
  bool realization;         //apply detector efficiency and resolution by dropping and smear
};

#endif
