/*
FastTracklet.cxx

Implementation of class Tracklet

Author: Kun Liu, liuk@fnal.gov
Created: 05-28-2013
*/

#include <iostream>
#include <algorithm>
#include <cmath>

#include <TMath.h>
#include <TMatrixD.h>

#include "FastTracklet.h"

ClassImp(SignedHit)
ClassImp(Tracklet)

SignedHit::SignedHit()
{
  hit.index = -1;
  hit.detectorID = -1;
  hit.elementID = -1;
  sign = 0;
}

SignedHit::SignedHit(int detectorID)
{
  hit.index = -1;
  hit.detectorID = detectorID;
  hit.elementID = -1;
  sign = 0;
}

SignedHit::SignedHit(Hit hit_input, int sign_input)
{
  hit = hit_input;
  sign = sign_input;
}

Tracklet::Tracklet()
{
  nXHits = 0;
  nUHits = 0;
  nVHits = 0;
  chisq = 9999.;

  tx = 0.;
  ty = 0.;
  x0 = 0.;
  y0 = 0.;
  invP = 0.02;

  err_tx = -1.;
  err_ty = -1.;
  err_x0 = -1.;
  err_y0 = -1.;
  err_invP = -1.;

  stationID = -1;

  for(int i = 0; i < 24; i++) residual[i] = 999.;
}

bool Tracklet::isValid()
{
  if(stationID < 1 || stationID > 6) return false;

  if(fabs(tx) > TX_MAX || fabs(x0) > X0_MAX) return false;
  if(fabs(ty) > TY_MAX || fabs(y0) > Y0_MAX) return false;
  if(err_tx < 0 || err_ty < 0 || err_x0 < 0 || err_y0 < 0) return false;

  double prob = getProb();
  if(prob < PROB_LOOSE) return false;

  //Tracklets in each station
  int nHits = nXHits + nUHits + nVHits;
  if(stationID < 5)
    {
      if(nXHits < 1 || nUHits < 1 || nVHits < 1) return false;
      if(nHits < 4) return false;
      if(chisq > 15.) return false;
    }

  //Back partial
  if(stationID == 5)
    {
      if(nXHits < 2 || nUHits < 2 || nVHits < 2) return false;
      if(nHits < 8) return false; 
    }

  //Global tracks
  if(stationID == 6)
    {
      if(nXHits < 3 || nUHits < 3 || nVHits < 3) return false;
      if(nHits < 12) return false;
      if(prob < PROB_TIGHT) return false;
      
      if(KMAG_ON == 1)
	{
	  if(invP < INVP_MIN || invP > INVP_MAX) return false;
       	}
    }

  return true;
}

double Tracklet::getProb() const
{
  int ndf;
  if(stationID == 6 && KMAG_ON == 1)
    {
      ndf = getNHits() - 5;
    }
  else
    {
      ndf = getNHits() - 4;
    }

  return TMath::Prob(chisq, ndf);
}

double Tracklet::getExpPositionX(double z) const
{
  if(KMAG_ON == 1 && stationID >= 5 && z < Z_KMAG_BEND - 1.)
    {
      double tx_st1 = tx + PT_KICK_KMAG*invP*getCharge();
      double x0_st1 = tx*Z_KMAG_BEND + x0 - tx_st1*Z_KMAG_BEND;

      return x0_st1 + tx_st1*z;
    }
  else
    {
      return x0 + tx*z;
    }
}

double Tracklet::getExpPosErrorX(double z) const
{
  double err_x;
  if(KMAG_ON == 1 && stationID >= 5 && z < Z_KMAG_BEND - 1.)
    {
      double err_kick = err_invP*PT_KICK_KMAG;
      double err_tx_st1 = err_tx + err_kick;
      double err_x0_st1 = err_x0 + err_kick*Z_KMAG_BEND;

      err_x = err_x0_st1 + fabs(err_tx_st1*z);
    }
  else
    {
      err_x = fabs(err_tx*z) + err_x0;
    }

  if(z > Z_ABSORBER) err_x += 1.;
  return err_x;
}

double Tracklet::getExpPositionY(double z) const
{
  return y0 + ty*z;
}

double Tracklet::getExpPosErrorY(double z) const
{
  double err_y = fabs(err_ty*z) + err_y0;
  if(z > Z_ABSORBER) err_y += 1.;

  return err_y;
}

double Tracklet::getExpPositionW(int detectorID)
{
  GeomSvc* p_geomSvc = GeomSvc::instance();
  double z = p_geomSvc->getPlanePosition(detectorID);

  double x_exp = getExpPositionX(z);
  double y_exp = getExpPositionY(z);

  return p_geomSvc->getCostheta(detectorID)*x_exp + p_geomSvc->getSintheta(detectorID)*y_exp;
}

bool Tracklet::operator<(const Tracklet& elem) const
{
  //return nXHits + nUHits + nVHits - 0.4*chisq > elem.nXHits + elem.nUHits + elem.nVHits - 0.4*elem.chisq;
  if(getNHits() == elem.getNHits()) 
    {
      return chisq < elem.chisq;
    }
  else
    {
      return getProb() > elem.getProb();
    }
}

bool Tracklet::similarity(const Tracklet& elem) const
{
  int nCommonHits = 0;
  std::list<SignedHit>::const_iterator first = hits.begin();
  std::list<SignedHit>::const_iterator second = elem.hits.begin();

  while(first != hits.end() && second != elem.hits.end())
    {
      if((*first) < (*second))
	{
	  ++first;
	}
      else if((*second) < (*first))
	{
	  ++second;
	}
      else
	{
	  if((*first) == (*second)) nCommonHits++;
	  ++first;
	  ++second;
	}
    }

  if(nCommonHits/double(elem.getNHits()) > 0.33333) return true;
  return false;
}

double Tracklet::getMomentum() const
{
  //Ref. SEAQUEST-doc-453-v3 by Don. Geesaman
  //if(KMAG_ON == 0) return 1E8;

  double p = 50.;
  double charge = getCharge();

  double c1 = Z_FMAG_BEND*PT_KICK_FMAG*charge;
  double c2 = Z_KMAG_BEND*PT_KICK_KMAG*charge;
  double c3 = -x0;
  double c4 = ELOSS_KFMAG;
  double c5 = c4/2.;

  double b = c1/c3 + c2/c3 - c4 - c5;
  double c = c4*c5 - c1*c5/c3 - c2*c4/c3;

  double disc = b*b - 4*c;
  if(disc > 0.)
    {
      p = (-b + sqrt(disc))/2. - ELOSS_KFMAG;
    }

  if(p < 10. || p > 120. || disc < 0)
    {
      double k = fabs(getExpPositionX(Z_KFMAG_BEND)/Z_KFMAG_BEND - tx);
      p = 1./(0.00832161 + 0.184186*k - 0.104132*k*k) + ELOSS_ABSORBER;
    }

  return p;
}

void Tracklet::getXZInfoInSt1(double& tx_st1, double& x0_st1)
{
  if(KMAG_ON == 1)
    {
      tx_st1 = tx + PT_KICK_KMAG*invP*getCharge();
      x0_st1 = tx*Z_KMAG_BEND + x0 - tx_st1*Z_KMAG_BEND;
    }
  else
    {
      tx_st1 = tx;
      x0_st1 = x0;
    }  
}

void Tracklet::getXZErrorInSt1(double& err_tx_st1, double& err_x0_st1)
{
  if(KMAG_ON == 1)
    {
      double err_kick = err_invP*PT_KICK_KMAG;
      err_tx_st1 = err_tx + err_kick;
      err_x0_st1 = err_x0 + err_kick*Z_KMAG_BEND;
    }
  else
    {
      err_tx_st1 = err_tx;
      err_x0_st1 = err_x0;
    }
}

Tracklet Tracklet::operator+(const Tracklet& elem) const
{
  Tracklet tracklet;
  tracklet.stationID = 5;

  tracklet.nXHits = nXHits + elem.nXHits;
  tracklet.nUHits = nUHits + elem.nUHits;
  tracklet.nVHits = nVHits + elem.nVHits;

  tracklet.hits.assign(hits.begin(), hits.end());
  if(elem.stationID > stationID)
    {
      tracklet.hits.insert(tracklet.hits.end(), elem.hits.begin(), elem.hits.end());
    }
  else
    {
      tracklet.hits.insert(tracklet.hits.begin(), elem.hits.begin(), elem.hits.end());
    }

  tracklet.err_tx = 1./sqrt(1./err_tx/err_tx + 1./elem.err_tx/elem.err_tx);
  tracklet.err_ty = 1./sqrt(1./err_ty/err_ty + 1./elem.err_ty/elem.err_ty);
  tracklet.err_x0 = 1./sqrt(1./err_x0/err_x0 + 1./elem.err_x0/elem.err_x0);
  tracklet.err_y0 = 1./sqrt(1./err_y0/err_y0 + 1./elem.err_y0/elem.err_y0);

  tracklet.tx = (tx/err_tx/err_tx + elem.tx/elem.err_tx/elem.err_tx)*tracklet.err_tx*tracklet.err_tx;
  tracklet.ty = (ty/err_ty/err_ty + elem.ty/elem.err_ty/elem.err_ty)*tracklet.err_ty*tracklet.err_ty;
  tracklet.x0 = (x0/err_x0/err_x0 + elem.x0/elem.err_x0/elem.err_x0)*tracklet.err_x0*tracklet.err_x0;
  tracklet.y0 = (y0/err_y0/err_y0 + elem.y0/elem.err_y0/elem.err_y0)*tracklet.err_y0*tracklet.err_y0;

  tracklet.invP = 1./tracklet.getMomentum();
  tracklet.err_invP = 0.25*tracklet.invP;
  
  tracklet.calcChisq();
  return tracklet;
}

Tracklet Tracklet::operator*(const Tracklet& elem) const
{
  Tracklet tracklet;
  tracklet.stationID = 6;

  tracklet.nXHits = nXHits + elem.nXHits;
  tracklet.nUHits = nUHits + elem.nUHits;
  tracklet.nVHits = nVHits + elem.nVHits;

  tracklet.hits.assign(hits.begin(), hits.end());
  if(elem.stationID > stationID)
    {
      tracklet.hits.insert(tracklet.hits.end(), elem.hits.begin(), elem.hits.end());
    }
  else
    {
      tracklet.hits.insert(tracklet.hits.begin(), elem.hits.begin(), elem.hits.end());
    }

  if(elem.stationID == 5)
    {
      tracklet.tx = elem.tx;
      tracklet.ty = elem.ty;
      tracklet.x0 = elem.x0;
      tracklet.y0 = elem.y0;
      tracklet.invP = 1./elem.getMomentum();

      tracklet.err_tx = elem.err_tx;
      tracklet.err_ty = elem.err_ty;
      tracklet.err_x0 = elem.err_x0;
      tracklet.err_y0 = elem.err_y0;
      tracklet.err_invP = 0.25*tracklet.invP;
    }
  else
    {
      tracklet.tx = tx;
      tracklet.ty = ty;
      tracklet.x0 = x0;
      tracklet.y0 = y0;
      tracklet.invP = 1./getMomentum();

      tracklet.err_tx = err_tx;
      tracklet.err_ty = err_ty;
      tracklet.err_x0 = err_x0;
      tracklet.err_y0 = err_y0;
      tracklet.err_invP = 0.25*tracklet.invP;
    }

  tracklet.calcChisq();
  return tracklet;
}

void Tracklet::addDummyHits()
{
  std::vector<int> detectorIDs_all;
  for(int i = stationID*6 - 5; i <= stationID*6; i++) detectorIDs_all.push_back(i);

  std::vector<int> detectorIDs_now;
  for(std::list<SignedHit>::const_iterator iter = hits.begin(); iter != hits.end(); ++iter)
    {
      detectorIDs_now.push_back(iter->hit.detectorID);
    }

  std::vector<int> detectorIDs_miss(6);
  std::vector<int>::iterator iter = std::set_difference(detectorIDs_all.begin(), detectorIDs_all.end(), detectorIDs_now.begin(), detectorIDs_now.end(), detectorIDs_miss.begin());
  detectorIDs_miss.resize(iter - detectorIDs_miss.begin());

  for(std::vector<int>::iterator iter = detectorIDs_miss.begin(); iter != detectorIDs_miss.end(); ++iter)
    {
      SignedHit dummy;
      dummy.hit.detectorID = *iter;

      hits.push_back(dummy);
    }

  sortHits();
}

double Tracklet::calcChisq()
{
  GeomSvc* p_geomSvc = GeomSvc::instance();
  chisq = 0.;

  double tx_st1, x0_st1;
  if(stationID == 6 && KMAG_ON == 1)
    {
      getXZInfoInSt1(tx_st1, x0_st1);
    }

  for(std::list<SignedHit>::const_iterator iter = hits.begin(); iter != hits.end(); ++iter)
    {
      if(iter->hit.index < 0) continue;

      int detectorID = iter->hit.detectorID;
      int index = detectorID - 1;

      double sigma;
#ifdef COARSE_MODE
      if(iter->sign == 0) sigma = p_geomSvc->getPlaneSpacing(detectorID)/sqrt(12.);
#else
      //if(iter->sign == 0) sigma = fabs(iter->hit.driftDistance);
      if(iter->sign == 0) sigma = p_geomSvc->getPlaneSpacing(detectorID)/sqrt(12.);
#endif
      if(iter->sign != 0) sigma = p_geomSvc->getPlaneResolution(detectorID);

      double p = iter->hit.pos + iter->sign*fabs(iter->hit.driftDistance);
      if(KMAG_ON == 1 && stationID == 6 && detectorID <= 6)
	{
    	  residual[index] = p - p_geomSvc->getInterception(detectorID, tx_st1, ty, x0_st1, y0);
	}
      else
	{
	  residual[index] = p - p_geomSvc->getInterception(detectorID, tx, ty, x0, y0);
	}
     
      chisq += (residual[index]*residual[index]/sigma/sigma);
      //std::cout << iter->hit.detectorID << "  " << iter->hit.elementID << "  " << iter->sign << "  " << iter->hit.pos << "  " << iter->hit.driftDistance << "  " << costheta << "  " << sintheta << "  " << z << "  " << (x0_st1 + tx_st1*z) << "  " << (x0 + tx*z) << "  " << (y0 + ty*z) << "  " << sigma << std::endl;
    }

  //std::cout << chisq << std::endl;
  return chisq;
}

SignedHit Tracklet::getSignedHit(int index)
{
  int id = 0;
  for(std::list<SignedHit>::const_iterator iter = hits.begin(); iter != hits.end(); ++iter)
    {
      if(id == index) return *iter;
      id++;
    }

  SignedHit dummy;
  return dummy;
}

double Tracklet::Eval(const double* par)
{
  tx = par[0];
  ty = par[1];
  x0 = par[2];
  y0 = par[3];
  if(KMAG_ON == 1) invP = par[4];

  //std::cout << tx << "  " << ty << "  " << x0 << "  " << y0 << "  " << 1./invP << std::endl;
  return calcChisq();
}

SRecTrack Tracklet::getSRecTrack()
{
  GeomSvc* p_geomSvc = GeomSvc::instance();

  SRecTrack strack;
  strack.setChisq(chisq);
  for(std::list<SignedHit>::iterator iter = hits.begin(); iter != hits.end(); ++iter)
    {
      if(iter->hit.index < 0) continue;

      double z = p_geomSvc->getPlanePosition(iter->hit.detectorID);
      double tx_val, tx_err, x0_val, x0_err;
      if(iter->hit.detectorID <= 6)
	{
	  getXZInfoInSt1(tx_val, x0_val);
          getXZErrorInSt1(tx_err, x0_err);
	}
      else
	{
	  tx_val = tx;
	  x0_val = x0;
	  tx_err = err_tx;
	  x0_err = err_x0;
	}

      TMatrixD state(5, 1), covar(5, 5);
      state[0][0] = getCharge()*invP*sqrt((1. + tx_val*tx_val)/(1. + tx_val*tx_val + ty*ty));
      state[1][0] = tx_val;
      state[2][0] = ty;
      state[3][0] = getExpPositionX(z);
      state[4][0] = getExpPositionY(z);

      covar.Zero();
      covar[0][0] = err_invP*err_invP;
      covar[1][1] = tx_err*tx_err;
      covar[2][2] = err_ty*err_ty;
      covar[3][3] = getExpPosErrorX(z)*getExpPosErrorX(z);
      covar[4][4] = getExpPosErrorY(z)*getExpPosErrorY(z);
      
      strack.insertHitIndex(iter->hit.index*iter->sign);
      strack.insertStateVector(state);
      strack.insertCovariance(covar);
      strack.insertZ(z);
    }

  strack.swimToVertex();
  strack.setHodoHits();
  return strack;
}

TVector3 Tracklet::getMomentumSt1()
{
  double tx_st1, x0_st1;
  getXZInfoInSt1(tx_st1, x0_st1);

  double pz = 1./invP/sqrt(1. + tx_st1*tx_st1);
  return TVector3(pz*tx_st1, pz*ty, pz);
}

TVector3 Tracklet::getMomentumSt3()
{
  double pz = 1./invP/sqrt(1. + tx*tx);
  return TVector3(pz*tx, pz*ty, pz);
}

void Tracklet::print()
{
  using namespace std;

  calcChisq();

  cout << "Tracklet in station " << stationID << endl;
  cout << nXHits + nUHits + nVHits << " hits in this station with chisq = " << chisq << endl; 
  cout << "Momentum in z: " << 1./invP << " +/- " << err_invP/invP/invP << endl;
  cout << "Charge: " << getCharge() << endl;
  for(std::list<SignedHit>::iterator iter = hits.begin(); iter != hits.end(); ++iter)
    {
      if(iter->sign > 0) cout << "L: ";
      if(iter->sign < 0) cout << "R: ";
      if(iter->sign == 0) cout << "U: ";

      cout << iter->hit.index << " " << iter->hit.detectorID << "  " << iter->hit.elementID << "  " << residual[iter->hit.detectorID-1] << " === ";
    }
  cout << endl;

  cout << "X-Z: (" << tx << " +/- " << err_tx << ")*z + (" << x0 << " +/- " << err_x0 << ")" << endl;
  cout << "Y-Z: (" << ty << " +/- " << err_ty << ")*z + (" << y0 << " +/- " << err_y0 << ")" << endl;
  
  cout << "KMAG projection: X =  " << getExpPositionX(Z_KMAG_BEND) << " +/- " << getExpPosErrorX(Z_KMAG_BEND) << endl;  
  cout << "KMAG projection: Y =  " << getExpPositionY(Z_KMAG_BEND) << " +/- " << getExpPosErrorY(Z_KMAG_BEND) << endl;  
}
