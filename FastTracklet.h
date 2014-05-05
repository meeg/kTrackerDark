/*

Supporting structure of fast tracking, include SignedHit and Tracklet

Author: Kun Liu, liuk@fnal.gov
Created: 06-09-2013

*/

#ifndef _FASTTRACKLET_H
#define _FASTTRACKLET_H

#include "MODE_SWITCH.h"

#include <list>
#include <vector>

#include <TObject.h>
#include <TVector3.h>

#include "GeomSvc.h"
#include "SRawEvent.h"
#include "SRecEvent.h"

class SignedHit : public TObject
{
public:
  SignedHit();
  explicit SignedHit(int detectorID);
  SignedHit(Hit hit_input, int sign_input);
  
  //comparision operators for sorting
  bool operator<(const SignedHit elem) const { return hit.detectorID < elem.hit.detectorID; }
  bool operator==(const SignedHit elem) const { return hit.index == elem.hit.index; }

  //Get the real hit position
  double pos() { return hit.pos + sign*hit.driftDistance; }
  double pos(int sign_input) { return hit.pos + sign_input*hit.driftDistance; }

  //Data members
  Hit hit;
  int sign;

  ClassDef(SignedHit, 1)
};

class PropSegment : public TObject
{
public:
  PropSegment();

  //init -- temporary, only used for tests
  void init();

  //Quality cut
  bool isValid();

  //Number of hits
  int getNHits();

  //Fit the segment -- naive linear fit
  void fit();

  //track slope the interception
  double a;
  double b;
  double err_a;
  double err_b;

  //chisq of the segment
  double chisq;

  //Hit list -- only 4 hits at most
  SignedHit hits[4];

  ClassDef(PropSegment, 1)
};

class Tracklet : public TObject
{
public:
  Tracklet();

  //Basic quality cut
  bool isValid();

  //Sort hit list
  void sortHits() { hits.sort(); }

  //Get number of real hits
  int getNHits() const { return nXHits + nUHits + nVHits; }

  //Get the probabilities
  double getProb() const;

  //Get the chi square 
  double getChisq() const { return chisq; }

  //Get x and y positions at a given z
  double getExpPositionX(double z) const;
  double getExpPosErrorX(double z) const;
  double getExpPositionY(double z) const;
  double getExpPosErrorY(double z) const;
  double getExpPositionW(int detectorID);

  //Get momentum upstream/downstream
  TVector3 getMomentumSt1();
  TVector3 getMomentumSt3();

  //Get the i-th signed hit
  SignedHit getSignedHit(int index);

  //Kernal function to calculate chi square for minimizer
  double Eval(const double* par);
  double calcChisq();

  //Add dummy hits
  void addDummyHits();

  //Momentum estimation using back partial
  double getMomentum() const;

  //Decide charge by KMag bending direction
  int getCharge() const { return (tx*Z_KMAG_BEND + x0)/Z_KMAG_BEND > tx ? 1 : -1; }

  //Get the slope and intersection in station 1
  void getXZInfoInSt1(double& tx_st1, double& x0_st1);
  void getXZErrorInSt1(double& err_tx_st1, double& err_x0_st1);

  //For sorting tracklet list
  bool operator<(const Tracklet& elem) const;

  //For reducing similar tracklets
  bool similarity(const Tracklet& elem) const;

  //For adding two tracklets together to form a back partial track
  Tracklet operator+(const Tracklet& elem) const;

  //For adding two tracklets together to form a global track
  Tracklet operator*(const Tracklet& elem) const;

  //Debuggin output
  void print();

  //Convert to a SRecTrack
  SRecTrack getSRecTrack();

  //remember if jobopts tells us kmag is one
  bool m_enableKMag;

  //Station ID
  int stationID;

  //Number of hits
  int nXHits;
  int nUHits;
  int nVHits;

  //Chi square
  double chisq;

  //List of signed hits
  std::list<SignedHit> hits;

  //Corresponding prop. tube segments
  PropSegment seg_x;
  PropSegment seg_y;

  //Slope, intersection, momentum and their errors
  double tx;
  double ty;
  double x0;
  double y0;
  double invP;

  double err_tx;
  double err_ty;
  double err_x0;
  double err_y0;
  double err_invP;

  //Residuals of all pos
  double residual[24];

  ClassDef(Tracklet, 3)
};


#endif
