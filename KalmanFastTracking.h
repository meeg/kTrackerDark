/*
KalmanFastTracking.h

Fast tracking utility of Kalman filter track, used to improve the tracking speed and also for online monitoring

Author: Kun Liu, liuk@fnal.gov
Created: 05-24-2013
*/

#ifndef _KALMANFASTTRACKING_H
#define _KALMANFASTTRACKING_H

#include "MODE_SWITCH.h"

#include <list>
#include <vector>

#include <Math/Factory.h>
#include <Math/Minimizer.h>
#include <Math/Functor.h>

#include "GeomSvc.h"
#include "SRawEvent.h"
#include "KalmanTrack.h"
#include "KalmanFitter.h"
#include "FastTracklet.h"

class KalmanFastTracking
{
public:
  KalmanFastTracking(bool flag = true);
  ~KalmanFastTracking();

  //Set the input event
  bool setRawEvent(SRawEvent* event_input);

  //Event quality cut
  bool acceptEvent(SRawEvent* rawEvent);

  ///Tracklet finding stuff
  //Build tracklets in a station
  void buildTrackletsInStation(int stationID, double* pos_exp = NULL, double* window = NULL);

  //Build back partial tracks using tracklets in station 2 & 3
  void buildBackPartialTracks();

  //Build global tracks by connecting station 23 tracklets and station 1 tracklets
  void buildGlobalTracks();

  //Fit tracklets
  int fitTracklet(Tracklet& tracklet);

  //Check the quality of tracklet, number of hits
  bool acceptTracklet(Tracklet& tracklet);
  bool muonID(Tracklet& tracklet);

  //Resolve left-right when possible
  void resolveLeftRight(SRawEvent::hit_pair hpair, int& LR1, int& LR2);
  void resolveLeftRight(Tracklet& tracklet, double threshold);
  void resolveSingleLeftRight(Tracklet& tracklet);

  //Remove bad hit if needed
  void removeBadHits(Tracklet& tracklet);

  //Reduce the list of tracklets, returns the number of elements reduced
  int reduceTrackletList(std::list<Tracklet>& tracklets);

  //Get exp postion and window using sagitta method in station 1
  void getSagittaWindowsInSt1(Tracklet& tracklet, double* pos_exp, double* window);
  void getExtrapoWindowsInSt1(Tracklet& tracklet, double* pos_exp, double* window);
  
  //Print the distribution of tracklets at detector back/front
  void printAtDetectorBack(int stationID, std::string outputFileName);

  ///Track fitting stuff
  //Convert Tracklet to KalmanTrack and solve left-right problem
  void processOneTracklet(Tracklet& tracklet);

  //Use Kalman fitter to fit a track
  bool fitTrack(KalmanTrack& kmtrk);

  //Resolve left right by Kalman fitting results
  void resolveLeftRight(KalmanTrack& kmtrk);

  ///Final output
  std::list<Tracklet>& getFinalTracklets() { return trackletsInSt[4]; }
  std::list<Tracklet>& getBackPartials() { return trackletsInSt[3]; }
  std::list<KalmanTrack>& getKalmanTracks() { return tracks; }

  ///Tool, a simple-minded chi square fit
  void chi2fit(int n, double x[], double y[], double& a, double& b);

private:
  //Raw event input
  SRawEvent* rawEvent;
  std::vector<Hit> hitAll;

  //Tracklets in one event, id = 0, 1, 2 for station 1, 2, 3+/-, id = 3 for station 2&3 combined, id = 4 for global tracks
  //Likewise for the next part
  std::list<Tracklet> trackletsInSt[5];

  //Final kalman tracks
  std::list<KalmanTrack> tracks;

  ///Configurations of tracklet finding
  //Hodo. IDs for masking
  std::vector<int> detectorIDs_mask[4];
  std::vector<int> detectorIDs_maskX[4];
  std::vector<int> detectorIDs_maskY[4];
  std::list<int> hitIDs_mask[4]; //hits in T/B, L/R are combined
  std::vector<int> stationIDs_mask[6];

  //prop. tube IDs for MUID -- 0 for x-z, 1 for y-z
  int detectorIDs_muid[2][4];
  std::list<int> hitIDs_muid[2][4];

  //Masking window sizes, index is the uniqueID defined by nElement*detectorID + elementID
  double z_mask[24];
  double x_mask_min[24][72];
  double x_mask_max[24][72];
  double y_mask_min[24][72];
  double y_mask_max[24][72];
  
  ///For following part, id = 0, 1, 2, 3 stand for station 1, 2, 3+, 3-
  //Super plane IDs for DCs
  std::vector<int> superIDs[4];
  
  //Window sizes for X-U combination
  double u_win[4];
  double u_costheta[4];
  double u_sintheta[4];

  //Plane angles for all planes
  double costheta_plane[25];
  double sintheta_plane[25];

  //Z positions
  double z_plane_x[4];
  double z_plane_u[4];
  double z_plane_v[4];
  double z_plane[25];

  //Maximum slope and intersection in each view
  double slope_max[25];
  double intersection_max[25];

  //Resolutions of all planes
  double resol_plane[25];

  //Cell width of all planes
  double spacing_plane[25];

  //Sagitta ratio in station 1 U/X/V
  double s_ratio[3];
  double s_sigma[3];
  int s_detectorID[3];

  //Current tracklets being processed
  Tracklet tracklet_curr;

  //Least chi square fitter and functor
  ROOT::Math::Minimizer* minimizer[2];
  ROOT::Math::Functor fcn;

  //Kalman fitter
  KalmanFitter* kmfitter;

  //Geometry service
  GeomSvc* p_geomSvc;

  //Flag for enable Kalman fitting
  bool enable_KF;
};

#endif
