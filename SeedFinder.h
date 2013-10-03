/*
SeedFinder.h

Definition of class Seed1D and SeedFinder, which is used to reconstruct
prop. tube seeds from sub-detectors as the seed input for kalman filter

Author: Kun Liu, liuk@fnal.gov
Created: 10-25-2011

*/

#ifndef _SEEDFINDER_H
#define _SEEDFINDER_H

#include "MODE_SWITCH.h"

#include <iostream>
#include <list>
#include <vector>
#include <utility>
#include <string>

#include "GeomSvc.h"
#include "SRawEvent.h"

class Seed1D
{
public:
  Seed1D();

  double ax, dax;
  double bx, dbx;
  double xchisq;

  std::list<int> xhits;

  int nHits; 
  double chisq;
  double quality;
  
  void update();
  
  double getExpPositionX(double z);
  double getExpPosErrorX(double z);

  bool hasHit(int hitID); 
  bool operator<(const Seed1D& elem) const;
  bool operator>(const Seed1D& elem) const;
  bool operator==(const Seed1D& elem) const;
  bool similarity(const Seed1D& elem) const;

  void print();
};

class SeedFinder
{
public:
  SeedFinder();

  ///Set event pointer
  void setEvent(SRawEvent *event_input);

  ///Process one event
  int processOneEvent(SRawEvent *event_input);

  ///Check the basic quality cuts
  bool acceptEvent();

  ///Requirement of a good seed
  bool acceptSeed(Seed1D& _seed);
  bool hodoMask(Seed1D& _seed);
  bool hodoMask(Seed1D& _seedx, Seed1D& _seedy);

  ///Initialize the seed candidates with random combinations
  void initializeSeedCandidates();

  ///Associate other hits
  void associateHits(Seed1D& _seed, double window, std::vector<int>& detectorIDs);

  ///Reduce the candidates list
  void reduceSeedList(std::list<Seed1D>& _seedlist_source, std::list<Seed1D>& _seedlist_target); 

  ///Add/remove a list of hits to exsiting seed 
  void addHitsToSeed(Seed1D& _seed, std::list<int>& _hitlist);
  int removeHitsFromSeed(Seed1D& _seed, std::list<int>& _hitlist_ref);

  ///Update the seed parameters
  void updateSeed(Seed1D& _seed);

  ///Direct linear fit
  int linearFit(double x[], double y[], double w[], int n, double *a, double *b, double *chisq, double *siga, double *sigb);

  ///Gets
  int getNSeeds() { return seed2d_final.size(); }
  std::list<Seed1D> getFinalSeeds() { return seed2d_final; }

  ///Print debug information
  void print();
  void printResults(std::string outputFileName);

  ///temp
  void setDetectorIDs(std::vector<int> start, std::vector<int> end, std::vector<int> all);
  void setMaskIDs(std::vector<int> mask) { detectors_mask.clear(); detectors_mask = mask; }
  void setSeedLimits(double ax, double bx) { ax_limit = ax; bx_limit = bx; }
   
private:
  SRawEvent *event;
  std::vector<Hit> hitAll;
  
  ///Seeds in x-z plane
  std::list<Seed1D> seed2d_candidates;
  std::list<Seed1D> seed2d_final;

  //Pattern recognization configurations
  std::vector<int> detectors_start_x;     //starting plane for initialization of 2D seeds
  std::vector<int> detectors_end_x;       //ending plane for initialization of 2D seeds
  std::vector<int> detectors_x;           //X/X' planes after the absorber
  std::vector<int> detectors_mask;

  //Seed limitations
  double ax_limit, bx_limit;

  ///Geometry service
  GeomSvc *geometrySvc;
};

#endif
