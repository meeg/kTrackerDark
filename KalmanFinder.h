/*
KalmanFinder.h

Use seed track from prop. tubes to find tracks
 
Author: Kun Liu, liuk@fnal.gov
Created: 10-14-2012
*/

#ifndef _KALMANFINDER_H
#define _KALMANFINDER_H

#include "MODE_SWITCH.h"

#include <list>
#include <vector>

#include <TFile.h>
#include <TTree.h>

#include "SRawEvent.h"
#include "KalmanUtil.h"
#include "KalmanTrack.h"

class KalmanFinder
{
public:
  KalmanFinder();
  ~KalmanFinder();

  ///Set event
  void setEvent(SRawEvent *ptr_evt);
  void setHitList(const std::vector<Hit>& hitAll);

  ///Event quality check
  bool acceptEvent();

  ///Main external call to process one seed
  int processOneSeed(Seed _seed);
  //int processOneSeed_single(Seed _seed);

  ///Initialize the kalman track finder
  void init();

  ///Search for hits in one detector Detector
  //void findHitsOnDetector(int detectorID, std::list<KalmanTrack>& _tracklist);
  void findHitsOnSuperDetector(int detectorID, std::list<KalmanTrack>& _tracklist); 
  double optimizedWindow(int detectorID, KalmanTrack& _track);

  bool addHitPairToTrack(SRawEvent::hit_pair _hitpair, int sign1, int sign2, KalmanTrack& _track);

  ///Refine track measurement with Kalman filter
  bool refineTrack(KalmanTrack& _track);

  ///Require adjacent wire only
  //bool reduceHitList(KalmanTrack& _track, std::list<int>& _hitlist);

  ///Reduce the track list by certain criteria
  void reduceTrackList(int detectorID, std::list<KalmanTrack>& _tracklist);
  void reduceTrackList(std::list<KalmanTrack>& _tracklist_source, std::list<KalmanTrack>& _tracklist_target);
  bool acceptTrack(int detectorID, KalmanTrack& _track);
  int getNHodoHits(int stationID, KalmanTrack& _track);

  ///Add single hits after entire search
  KalmanTrack associateSingles(KalmanTrack _track);

  ///Return the final tracks
  KalmanTrack getBestCandidate();
  std::list<KalmanTrack>& getAllCandidates() { return _track_candidates; }

  ///Print the output event display
  void printResults(std::string outputFileName, std::list<KalmanTrack>& _tracks);
  void printResults(std::string outputFileName, KalmanTrack& _track);
  
  ///Evaluation process
  void initEvaluation(std::string outputFileName);
  void finishEvaluation();
  void fillEvaluation(int trackID, KalmanTrack& _track);

private:
  ///Pointer to the raw event, and the hit list
  SRawEvent *event;
  std::vector<Hit> _hits;

  ///List of all the track candidates
  std::list<KalmanTrack> _track_candidates;

  ///Hodoscopes masks for each station
  std::vector<int> detectors_mask[3];
  std::list<int> hits_mask[3][4];

  ///Pointer to the geometry service
  GeomSvc *p_geomSvc;

  ///Evaluation of the track finding process
  TFile *_eval_file;
  TTree *_res_tree;

  ///residual tree
  int _runID;
  int _eventID;
  int _trackID;
  int _detectorID;
  double _mom_predicted;
  double _res_predicted;
  double _res_cov_predicted;
};

#endif
