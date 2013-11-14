/*
SeedFinder.cxx

Implimentation of class Seed1D and SeedFinder

Author: Kun Liu, liuk@fnal.gov
Created: 10-25-2011
*/


#include <iostream>
#include <list>
#include <vector>
#include <cmath>
#include <algorithm>

#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>

#include "GeomSvc.h"
#include "SeedFinder.h"

Seed1D::Seed1D()
{
  ax = -999.;
  bx = -999.;
  dax = 100.;
  dbx = 100.;
  xchisq = 0.;

  xhits.clear();

  update();
}

double Seed1D::getExpPositionX(double z)
{
  return ax*z+bx;
}

double Seed1D::getExpPosErrorX(double z)
{
  return dax*z+dbx;
}

void Seed1D::update()
{
  nHits = xhits.size();
  chisq = 0.;
  
  if(xchisq > 0) chisq += xchisq;
  quality = nHits - 0.2*chisq;
}

bool Seed1D::hasHit(int hitID)
{
  std::list<int>::iterator iter = find(xhits.begin(), xhits.end(), hitID);
  if(iter != xhits.end())
    {
      return true;
    }

  return false;
}

///Comparisons below are made first on number of hits, then chi square.
bool Seed1D::operator==(const Seed1D& elem) const
{
  return fabs(quality - elem.quality) < 0.001; 
}

bool Seed1D::operator<(const Seed1D& elem) const
{
  return quality > elem.quality;
}

bool Seed1D::operator>(const Seed1D& elem) const
{
  return quality < elem.quality;
}

bool Seed1D::similarity(const Seed1D& elem) const
{
  double delta_ax = fabs(ax - elem.ax);
  double delta_bx = fabs(bx - elem.bx);

  double sigma_ax = sqrt(dax*dax + elem.dax*elem.dax);
  double sigma_bx = sqrt(dbx*dbx + elem.dbx*elem.dbx);

  if(delta_ax < sigma_ax && delta_bx < sigma_bx) return true;
  return false;
}

void Seed1D::print()
{
  std::cout << "=========== Info of this Seed1D ============" << std::endl; 
  std::cout << "x-z info" << std::endl;
  std::cout << "In total has " << xhits.size() << " hits in X/X' layer." << std::endl;
  for(std::list<int>::iterator iter = xhits.begin(); iter != xhits.end(); ++iter)
    {
      std::cout << *iter << "  :  ";
    }
  std::cout << std::endl;
  std::cout << "ax = " << ax << ", bx = " << bx << ", chisq_x = " << xchisq << std::endl;
  std::cout << "dax = " << dax << ", dbx = " << dbx << ", chisq_x = " << xchisq << std::endl;
}

SeedFinder::SeedFinder()
{
  //Initialize the geometry service
  geometrySvc = GeomSvc::instance();
}

void SeedFinder::setEvent(SRawEvent *event_input)
{
  event = event_input;
  hitAll = event->getAllHits();

  seed2d_candidates.clear();
  seed2d_final.clear();
}

///Main external call
int SeedFinder::processOneEvent(SRawEvent *event_input)
{
  setEvent(event_input);
  if(!acceptEvent())
    {
      Log("Event quality is not good. Will continue to next event ... ");
      return 0;
    }

  ///Initialize x-z Seed candidates
  initializeSeedCandidates();
  
#ifdef _DEBUG_ON
  Log(seed2d_candidates.size() << " candidates in the seed init");
  for(std::list<Seed1D>::iterator iter = seed2d_candidates.begin(); iter != seed2d_candidates.end(); ++iter)
    {
      iter->print();
      Log("===================");
    }
#endif

  ///Associate the hits one other planes to the 2D Track candidates
  for(std::list<Seed1D>::iterator iter = seed2d_candidates.begin(); iter != seed2d_candidates.end(); ++iter)
    {
      associateHits(*iter, 0.9, detectors_x);
    }
  
  seed2d_candidates.sort();
  seed2d_candidates.unique();

#ifdef _DEBUG_ON
  Log(seed2d_candidates.size() << " candidates in the seed associate");
  Log(seed2d_candidates.size() << " 2D candidates at stage 1: ");
  for(std::list<Seed1D>::iterator iter = seed2d_candidates.begin(); iter != seed2d_candidates.end(); ++iter)
    {
      iter->print();
    }
#endif
 
  ///Remove the duplicate track candidates 
  if(seed2d_candidates.empty()) return 0;
  reduceSeedList(seed2d_candidates, seed2d_final);
 
#ifdef _DEBUG_ON
  Log(seed2d_final.size() << " reduced 2D candidates at stage 2: ");
  for(std::list<Seed1D>::iterator iter = seed2d_final.begin(); iter != seed2d_final.end(); ++iter)
    {
      iter->print();
    }
#endif

  ///Check the track quality of the final 2d tracks
  for(std::list<Seed1D>::iterator iter = seed2d_final.begin(); iter != seed2d_final.end(); )
    {
      if(!acceptSeed(*iter))
	{
	  iter = seed2d_final.erase(iter);
	}
      else
	{
	  iter++;
	}
    }

  return seed2d_final.size();
}

///Quality cuts for a event to be processed
bool SeedFinder::acceptEvent()
{
  int nHits_start_x = event->getNHitsInDetectors(detectors_start_x);
  int nHits_end_x = event->getNHitsInDetectors(detectors_end_x);
  //int nHits_x = event->getNHitsInDetectors(detectors_x);
 
#ifdef _DEBUG_ON
  Log(nHits_start_x << "  " << nHits_end_x);
#endif
  if(nHits_start_x < 2) return false;
  if(nHits_end_x < 2) return false;

  if(nHits_start_x > 20) return false;
  if(nHits_end_x > 20) return false;

  return true;  
}

//Quality cuts for a Track to be accepted
bool SeedFinder::acceptSeed(Seed1D& _seed)
{
  if(_seed.xhits.empty())
    { 
      //Log("No Hits.");
      return false;
    }

  if(_seed.xhits.size() < 4)
    {
      //Log("Les than 3 Hits.");
      return false;
    }
  if(_seed.xchisq > 40) return false;

  if(fabs(_seed.ax) > ax_limit) return false;
  if(fabs(_seed.bx) > bx_limit) return false;

  double z_mag = geometrySvc->getKMAGCenter();
  double x_mag = _seed.getExpPositionX(z_mag);
  double y_mag = 0.;

  if(!geometrySvc->isInKMAG(x_mag, y_mag))
    {
      //Log("Track projection is outside of KMAG!");
      return false;
    }

  if(!hodoMask(_seed)) 
    {
      //Log("Hodo check failed");
      return false;
    }

  return true; 
}

bool SeedFinder::hodoMask(Seed1D& _seedx, Seed1D& _seedy)
{
  int nMasks = 0;
  for(unsigned int i = 0; i < detectors_mask.size(); i++)
    {
      double z_exp = geometrySvc->getPlanePosition(detectors_mask[i]);
      double x_exp = _seedx.getExpPositionX(z_exp);
      double y_exp = _seedy.getExpPositionX(z_exp);
      double x_err = 0.1*geometrySvc->getPlaneSpacing(detectors_mask[i]);
      double y_err = 0.1*geometrySvc->getPlaneSpacing(detectors_mask[i]);;

      if(!geometrySvc->isInPlane(detectors_mask[i], x_exp, y_exp)) continue;

      std::list<int> hits_mask = event->getHitsIndexInDetector(detectors_mask[i]);
      if(hits_mask.empty()) continue;

      for(std::list<int>::iterator iter = hits_mask.begin(); iter != hits_mask.end(); ++iter)
	{
    	  double x_min, x_max, y_min, y_max;
	  geometrySvc->get2DBoxSize(hitAll[*iter].detectorID, hitAll[*iter].elementID, x_min, x_max, y_min, y_max);
          x_min -= x_err; x_max += x_err;
          y_min -= y_err; y_max += y_err;

	  if(x_exp > x_min && x_exp < x_max && y_exp > y_min && y_exp < y_max)
	    {
	      nMasks++;
	      if((i & 1) == 0) 
		{
		  //Log("Skpping next plane!");
		  i++;
		}

	      break;
	    }
	}
    }

  if(nMasks == int(detectors_mask.size()/2)) return true; 
  return false;
}

bool SeedFinder::hodoMask(Seed1D& _seed)
{
  std::string detectorName_first = geometrySvc->getDetectorName(hitAll[_seed.xhits.front()].detectorID);
  std::string detectorType;
  if(detectorName_first.find("X") != std::string::npos)
    {
      detectorType = "X";
    }
  else
    {
      detectorType = "Y";
    }

  int nMasks = 0;
  for(unsigned int i = 0; i < detectors_mask.size(); i++)
    {
      std::string detectorName = geometrySvc->getDetectorName(detectors_mask[i]);
      if(detectorName.find(detectorType.c_str()) == std::string::npos) continue;

      double z_exp = geometrySvc->getPlanePosition(detectors_mask[i]);
      double pos_exp = _seed.getExpPositionX(z_exp);

      std::list<int> hits_mask = event->getHitsIndexInDetector(detectors_mask[i]);
      if(hits_mask.empty()) continue;

      for(std::list<int>::iterator iter = hits_mask.begin(); iter != hits_mask.end(); ++iter)
	{
    	  double x_min, x_max, y_min, y_max;
	  geometrySvc->get2DBoxSize(hitAll[*iter].detectorID, hitAll[*iter].elementID, x_min, x_max, y_min, y_max);

	  double costheta = geometrySvc->getCostheta(hitAll[*iter].detectorID);
	  double sintheta = geometrySvc->getSintheta(hitAll[*iter].detectorID);
	  double pos_min = x_min*costheta + y_min*sintheta;
	  double pos_max = x_max*costheta + y_max*sintheta;

	  if(pos_exp > pos_min && pos_exp < pos_max)
	    {
	      nMasks++;
	      if((i & 1) == 0) 
		{
		  i++;
		}

	      break;
	    }
	}
    }

  if(nMasks == 2) return true;
  //if(nMasks > 0) return true;
  return false;
}

///Random combinations of two hits from the first and last plane avaible is used to 
///make Track candidates.
//
///After that all hits on the track is added to the Track, and the Track is updated
///Here double counting of hits is allowed
void SeedFinder::initializeSeedCandidates()
{
  std::list<int> hitID_start, hitID_end;
  hitID_start = event->getHitsIndexInDetectors(detectors_start_x);
  hitID_end = event->getHitsIndexInDetectors(detectors_end_x);

#ifdef _DEBUG_ON
  for(std::list<int>::iterator iter = hitID_start.begin(); iter != hitID_start.end(); ++iter)
    {
      std::cout << *iter << " === ";
    }
  std::cout << std::endl;
  
  for(std::list<int>::iterator iter = hitID_end.begin(); iter != hitID_end.end(); ++iter)
    {
      std::cout << *iter << " === ";
    }
  std::cout << std::endl;
#endif

  std::list<int>::iterator iter_i, iter_j;
  for(iter_i = hitID_start.begin(); iter_i != hitID_start.end(); ++iter_i)
    {
      for(iter_j = hitID_end.begin(); iter_j != hitID_end.end(); ++iter_j)
	{
	  //hitAll[*iter_i].print(); hitAll[*iter_j].print();
	  //std::cout << " === " << std::endl;

	  Seed1D _seed;

	  double x1 = hitAll[*iter_i].pos;
	  double x2 = hitAll[*iter_j].pos;
	  double z1 = geometrySvc->getPlanePosition(hitAll[*iter_i].detectorID);
	  double z2 = geometrySvc->getPlanePosition(hitAll[*iter_j].detectorID);

	  _seed.ax = (x1 - x2)/(z1 - z2);
	  _seed.bx = x1 - _seed.ax*z1;
	  _seed.xchisq = 0.;

	  _seed.xhits.push_back(*iter_i);
	  _seed.xhits.push_back(*iter_j);

          if(fabs(_seed.ax) < ax_limit && fabs(_seed.bx) < bx_limit)
	    {
	      seed2d_candidates.push_back(_seed);
	    }
	}
    }
}

void SeedFinder::associateHits(Seed1D& _seed, double window, std::vector<int>& detectorIDs)
{
  //_seed.print();
  unsigned int nDetectors = detectorIDs.size();
  for(unsigned int i = 0; i < nDetectors; i++)
    {	
      //LogDebug("Adding hits in detector: " << detectorIDs[i] << " named " << geometrySvc->getDetectorName(detectorIDs[i])); 
      double x_exp = _seed.getExpPositionX(geometrySvc->getPlanePosition(detectorIDs[i]));
  
      std::string detectorType = geometrySvc->getDetectorName(detectorIDs[i]);
      double range;
      if(detectorType.find("H") != std::string::npos)
	{
	  range = 0.49999999;
	}
      else
	{
	  range = 1.;
	}
      if(detectorIDs[i] < 35) range = range*1.5;

      std::list<int> _hitlist = event->getHitsIndexInDetector(detectorIDs[i], x_exp, range*geometrySvc->getPlaneSpacing(detectorIDs[i]));
	
      //LogDebug("Found " << _hitlist.size() << " hits in this detector!");
      if(_hitlist.size() > 0) addHitsToSeed(_seed, _hitlist);
    }	  
}

void SeedFinder::reduceSeedList(std::list<Seed1D>& _trklist_source, std::list<Seed1D>& _trklist_target)
{    
  _trklist_source.sort();
  /*
  for(std::list<Seed1D>::iterator iter = _trklist_source.begin(); iter != _trklist_source.end(); ++iter)
    {
      Log("The quality of this track is : " << iter->quality << ", nHits = " << iter->nHits << ", chisq = " << iter->chisq);
      iter->print();
    }
  */

  while(!_trklist_source.empty())
    {
      _trklist_target.push_back(_trklist_source.front());
      _trklist_source.pop_front();
      
      for(std::list<Seed1D>::iterator iter = _trklist_source.begin(); iter != _trklist_source.end(); )
	{
	  //Log("Working on this seed");
	  //iter->print();
	  if(removeHitsFromSeed(*iter, _trklist_target.back().xhits) <= 2)
	     {
	       //Log("Removed!");
	       iter = _trklist_source.erase(iter);
	       continue;
	     }

	  if(iter->similarity(_trklist_target.back()))
	    {
	      //Log("Removed!");
	      iter = _trklist_source.erase(iter);
	      continue;
	    }
  
	  ++iter;
	}
    }
}

void SeedFinder::addHitsToSeed(Seed1D& _seed, std::list<int>& _hitlist)
{
  //LogDebug("");
  if(_hitlist.empty()) return;

  unsigned int nHits_before = _seed.xhits.size();
  _seed.xhits.insert(_seed.xhits.end(), _hitlist.begin(), _hitlist.end());
  _seed.xhits.sort();
  _seed.xhits.unique();
  if(_seed.xhits.size() > nHits_before)
    {
      /*
      for(std::list<int>::iterator iter = hits_original->begin(); iter != hits_original->end(); ++iter)
	{
	  std::cout << *iter << " === ";
	}
      std::cout << std::endl;
      */

      updateSeed(_seed);
    }
}

int SeedFinder::removeHitsFromSeed(Seed1D& _seed, std::list<int>& _hitlist)
{
  std::list<int> commonHits; 
  commonHits.clear();	  
  set_intersection(_hitlist.begin(), _hitlist.end(), _seed.xhits.begin(), _seed.xhits.end(), back_inserter(commonHits));
	  
  if(_seed.xhits.size() - commonHits.size() > 2)	    
    {      
      for(std::list<int>::iterator iter = commonHits.begin(); iter != commonHits.end(); ++iter)
	{
	  _seed.xhits.remove(*iter);
	}
    }
  else
    {	
      return _seed.xhits.size() - commonHits.size();	    
    }

  if(!commonHits.empty())
    {
      updateSeed(_seed);
    }

  return _seed.xhits.size();
}

void SeedFinder::updateSeed(Seed1D& _seed)
{
  double xy[100], z[100], w[100];
  int nHits = 0;
  for(std::list<int>::iterator iter = _seed.xhits.begin(); iter != _seed.xhits.end(); ++iter)
    {
      double err = geometrySvc->getPlaneSpacing(hitAll[*iter].detectorID)/sqrt(12.);

      xy[nHits] = hitAll[*iter].pos;
      z[nHits] = geometrySvc->getPlanePosition(hitAll[*iter].detectorID);
      w[nHits] = 1./err/err;

      nHits++;
    }

  double _a, _b, _chisq, siga, sigb;
  linearFit(z, xy, w, nHits, &_a, &_b, &_chisq, &siga, &sigb);

  _seed.ax = _a;
  _seed.bx = _b;
  _seed.dax = siga;
  _seed.dbx = sigb;
  _seed.xchisq = _chisq;

  _seed.update();
}

///A least chisq fitter applied to a straight line
int SeedFinder::linearFit(double x[], double y[], double w[], int n, double *a, double *b, double *chisq, double *siga, double *sigb)
{
  double sum, sx, sy, sxx, sxy, syy, det;
  double chi;
  
  if(n < 2)
    {
      std::cout << "Should have at least two points!!!" << std::endl;
      return -1; 
    }
  
  sum = 0.;
  sx = 0.;
  sy = 0.;
  sxx = 0.;
  syy = 0.;
  sxy = 0.;
  
  for(int i = 0; i < n; i++)
    {
      sum += w[i];
      sx += w[i]*x[i];
      sy += w[i]*y[i];
      sxx += w[i]*x[i]*x[i];
      syy += w[i]*y[i]*y[i];
      sxy += w[i]*x[i]*y[i];
    }
  
  det = sum*sxx - sx*sx;
  if(fabs(det) < 1.0e-20)
    {
      *a = 1.0e20;
      *b = x[0];
      *chisq = 0.;
      *siga = 0.;
      *sigb = 0.;
      
      return 0;
    }
  
  *a = (sum*sxy - sx*sy)/det;
  *b = (sy*sxx - sxy*sx)/det;
  *siga = sqrt(fabs(sum/det));
  *sigb = sqrt(fabs(sxx/det));
  
  chi = 0.;
  for(int i = 0; i < n; i++)
    {
      chi += w[i]*(y[i] - (*a)*x[i] - (*b))*(y[i] - (*a)*x[i] - (*b));
    }
  
  *chisq = chi;
  
  return 1;  
}

void SeedFinder::print()
{
  std::cout << "The tracking info for Run " << event->getRunID() << ", Event " << event->getEventID() << std::endl; 
  std::cout << "This event has " << seed2d_final.size() << " successful Tracks. " << std::endl;
  int nTrack = 0;
  for(std::list<Seed1D>::iterator iter = seed2d_final.begin(); iter != seed2d_final.end(); ++iter)
    {
      std::cout << " == Track id. " << nTrack++ << ":" << std::endl;
      iter->print();
      std::cout << " == " << std::endl;
    }
  std::cout << "End of 3D final Track info. " << std::endl;
}

void SeedFinder::printResults(std::string outputFileName)
{
  double pos[3000], z[3000];
  double err[3000], dz[3000];

  std::list<int> hits = event->getHitsIndexInDetectors(detectors_x);
  std::string detectorName_first = geometrySvc->getDetectorName(hitAll[hits.front()].detectorID);
  std::string detectorType;
  if(detectorName_first.find("X") != std::string::npos)
    {
      detectorType = "X";
    }
  else
    {
      detectorType = "Y";
    }

  int nHits = 0;
  for(std::list<int>::iterator iter = hits.begin(); iter != hits.end(); ++iter)
    {
      pos[nHits] = hitAll[*iter].pos;
      z[nHits] = geometrySvc->getPlanePosition(hitAll[*iter].detectorID);
      //Log(z[nHits] << "  " << hitAll[*iter].detectorID);

      err[nHits] = 0.5*geometrySvc->getPlaneSpacing(hitAll[*iter].detectorID);
      dz[nHits] = 0.;

      nHits++;
    }

  /*
  std::list<int> hodo_hits = event->getHitsIndexInDetectors(detectors_mask);
  for(std::list<int>::iterator iter = hodo_hits.begin(); iter != hodo_hits.end(); ++iter)
    {
      std::string detectorName = geometrySvc->getDetectorName(hitAll[*iter].detectorID);
      if(detectorName.find(detectorType.c_str()) != std::string::npos)
	{
	  pos[nHits] = hitAll[*iter].pos;
	  err[nHits] = 0.5*geometrySvc->getPlaneSpacing(hitAll[*iter].detectorID);
	}
      else
	{
	  if(detectorType.find("X") != std::string::npos)
	    {
	      pos[nHits] = geometrySvc->getPlaneCenterX(hitAll[*iter].detectorID);
	      err[nHits] = 0.5*geometrySvc->getPlaneScaleX(hitAll[*iter].detectorID);
	    }
	  else
	    {
	      pos[nHits] = geometrySvc->getPlaneCenterY(hitAll[*iter].detectorID);
	      err[nHits] = 0.5*geometrySvc->getPlaneScaleY(hitAll[*iter].detectorID);
	    }
	}

      z[nHits] = geometrySvc->getPlanePosition(hitAll[*iter].detectorID);
      dz[nHits] = 0.;

      nHits++;
    }
  */
  TGraphErrors gr(nHits, z, pos, dz, err);

  TF1 f[3000];
  int nTracks = seed2d_final.size();
  for(std::list<Seed1D>::iterator iter = seed2d_final.begin(); iter != seed2d_final.end(); ++iter)
    {
      f[nTracks] = TF1("", "[0]+[1]*x", 1800, 2400);
      f[nTracks].SetParameter(0, iter->bx);
      f[nTracks].SetParameter(1, iter->ax);

      nTracks++;
    }

  TCanvas c1;
  c1.cd(); gr.Draw("AP");
  for(int i = 0; i < nTracks; i++)
    {
      f[i].Draw("same");
    }

  c1.Print(outputFileName.c_str());
}

void SeedFinder::setDetectorIDs(std::vector<int> start, std::vector<int> end, std::vector<int> all)
{
  detectors_start_x.clear();
  detectors_end_x.clear();
  detectors_x.clear();

  detectors_start_x = start;
  detectors_end_x = end;
  detectors_x = all;
}
