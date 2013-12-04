/*
KalmanFinder.cxx

Implementation of class Seed, KalmanTrack, and KalmanFinder

Author: Kun Liu, liuk@fnal.gov
Created: 10-14-2012
*/

#include <iostream>
#include <cmath>
#include <algorithm>

#include <TGraphErrors.h>
#include <TGraph.h>
#include <TCanvas.h>

#include "KalmanFinder.h"
#include "KalmanFitter.h"

KalmanFinder::KalmanFinder()
{
  init();

  p_geomSvc = GeomSvc::instance();

  detectors_mask[0] = p_geomSvc->getDetectorIDs("H1");
  detectors_mask[1] = p_geomSvc->getDetectorIDs("H2");
  detectors_mask[2] = p_geomSvc->getDetectorIDs("H3");

  _eval_file = NULL;
  _res_tree = NULL;
}

KalmanFinder::~KalmanFinder()
{
}

void KalmanFinder::setEvent(SRawEvent *ptr_evt)
{
  event = ptr_evt;
  setHitList(event->getAllHits());

  for(int i = 0; i < 3; i++)
    {
      for(unsigned int j = 0; j < detectors_mask[i].size(); j++)
	{
	  hits_mask[i][j].clear();
	  hits_mask[i][j] = event->getHitsIndexInDetector(detectors_mask[i][j]);
	}
    }
}

bool KalmanFinder::acceptEvent()
{
  ///Cuts on occupancy of DCs
  double occupancy_cut[4] = {0.15, 0.1, 0.1, 0.1};
  for(int i = 0; i < 4; i++)
    {
      for(int j = 1; j <= 6; j++)
	{
	  int detectorID = i*6 + j;
	  if(double(event->getNHitsInDetector(detectorID))/p_geomSvc->getPlaneNElements(detectorID) > occupancy_cut[i])
	    {
	      return false;
	    }
	}
    }

  ///Cuts on hodos
  for(int i = 0; i < 3; i++)
    {
      for(unsigned int j = 0; j < detectors_mask[i].size(); j++)
	{
	  if(hits_mask[i][j].empty()) return false;
	}
    }

  return true;
}

void KalmanFinder::setHitList(const std::vector<Hit>& hitAll)
{
  _hits.clear();
  _hits.assign(hitAll.begin(), hitAll.end());
 
#ifdef _DEBUG_ON
  for(std::vector<Hit>::iterator iter = _hits.begin(); iter != _hits.end(); ++iter)
    {
      iter->print();
    }
#endif
}

void KalmanFinder::init()
{
  _track_candidates.clear();
}

int KalmanFinder::processOneSeed(Seed _seed)
{
  init();

  ///Construct the first KalmanTrack from Seed
  KalmanTrack kmtrk_init(_seed);

  ///The search starts on ST3+ and ST3- separately, then combined to go for ST2&1
  std::list<KalmanTrack> _track_st3p;
  _track_st3p.push_back(kmtrk_init);
  for(int i = 9; i > 6; i--)
    {
      findHitsOnSuperDetector(i, _track_st3p);
      reduceTrackList(2*i-1, _track_st3p);
    }

  std::list<KalmanTrack> _track_st3m;
  _track_st3m.push_back(kmtrk_init);
  for(int i = 12; i > 9; i--)
    {
      findHitsOnSuperDetector(i, _track_st3m);
      reduceTrackList(2*i-1, _track_st3m);
    }

  ///Combine the track candidate list from St3+/-, and sort according to quality
  _track_candidates.insert(_track_candidates.end(), _track_st3p.begin(), _track_st3p.end());
  _track_candidates.insert(_track_candidates.end(), _track_st3m.begin(), _track_st3m.end());
  _track_candidates.sort();
  _track_candidates.unique();

  ///5 best candidates are kept after St3
  if(_track_candidates.empty()) return 0; 
  if(_track_candidates.size() > 5) _track_candidates.resize(5);

  ///Refine the track parameter before propagate to ST2
  for(std::list<KalmanTrack>::iterator iter = _track_candidates.begin(); iter != _track_candidates.end(); ++iter)
    {
      refineTrack(*iter);
    }
  _track_candidates.sort();

  ///Search on ST2
  for(int i = 6; i > 3; i--)
    {
      findHitsOnSuperDetector(i, _track_candidates);
      reduceTrackList(2*i-1, _track_candidates);
    }

  ///5 best candidates are kept for refinement
  if(_track_candidates.empty()) return 0;
  if(_track_candidates.size() > 5) _track_candidates.resize(5);
  for(std::list<KalmanTrack>::iterator iter = _track_candidates.begin(); iter != _track_candidates.end(); ++iter)
    { 
      refineTrack(*iter);
    }
  _track_candidates.sort();

  ///Only one candidate after searching on ST3 and ST2
  if(_track_candidates.size() > 1) _track_candidates.resize(1);

  ///Add a track with the same hit list but opposite charge
  KalmanTrack kmtrk_opposite = _track_candidates.front();
  kmtrk_opposite.flipCharge();
  refineTrack(kmtrk_opposite);
  _track_candidates.push_back(kmtrk_opposite);

  ///Update the momentum estimation using the precise position at ST2
  /// Because in ST3 and 2 the initial estimation won't be refined much
  for(std::list<KalmanTrack>::iterator iter = _track_candidates.begin(); iter != _track_candidates.end(); ++iter)
    {
      iter->updateMomentum();
    }

  ///Search on ST1
  for(int i = 3; i > 0; i--)
    {
      findHitsOnSuperDetector(i, _track_candidates);
      reduceTrackList(2*i-1, _track_candidates);

      if(_track_candidates.empty()) return 0;
    }
  _track_candidates.sort();
  
  return _track_candidates.size();
}

void KalmanFinder::findHitsOnSuperDetector(int detectorID, std::list<KalmanTrack>& _tracklist)
{ 
#ifdef _DEBUG_ON  
  Log("Looking for hits on " << 2*detectorID << " and " << 2*detectorID-1);   
  Log("Currently we have " << _tracklist.size() << " candidates");
#endif

  ///Retrieve the paired hit list on specific detector super plane
  std::list<SRawEvent::hit_pair> _pairlist_all = event->getHitPairsInSuperDetector(detectorID);
  
#ifdef _DEBUG_ON
  Log("Totally " << _pairlist_all.size() << " hit pairs in this Detector");
  for(std::list<SRawEvent::hit_pair>::iterator iter = _pairlist_all.begin(); iter != _pairlist_all.end(); ++iter)
    {
      Log("======================");
      _hits[iter->first].print();
      _hits[iter->second].print();
    }
#endif
  
  if(_pairlist_all.empty()) return;

  int iTrack = 0;
  std::list<KalmanTrack> _track_new;
  for(std::list<KalmanTrack>::iterator kmtrk = _tracklist.begin(); kmtrk != _tracklist.end(); )
    {
#ifdef _DEBUG_ON
      Log("Start working for a new candidate ... ");
      kmtrk->print();
#endif

      ///If the track is out of the world, remove the entire candidate
      if(!kmtrk->propagateTo(2*detectorID))
	{
    	  kmtrk = _tracklist.erase(kmtrk);
	  continue;
	}

      ///Open a approapriate window around the expected hit position
      double pos_exp, window;
      pos_exp = kmtrk->getExpPosition();
      window = optimizedWindow(detectorID, *kmtrk);
      std::list<SRawEvent::hit_pair> _pairlist = event->getHitPairsInSuperDetector(detectorID, pos_exp, window);

#ifdef _DEBUG_ON
      Log("Expected position is " << pos_exp << " +/- " << window);
      Log("Found " << _pairlist.size() << " paired hits!");
#endif

      if(_pairlist.empty())
	{
	  ++kmtrk;
	  continue;
	}
      
      ///Construct a new candidate track with each hit pair added
      for(std::list<SRawEvent::hit_pair>::iterator iter = _pairlist.begin(); iter != _pairlist.end(); ++iter)
	{
	  ///If the two tracks are off by more than one spacing then skip
	  if(fabs(_hits[iter->first].pos - _hits[iter->second].pos) > 1.1*p_geomSvc->getPlaneSpacing(_hits[iter->first].detectorID))
	    {
	      continue;
	    }

#ifdef COARSE_MODE
	  KalmanTrack kmtrk_new = *kmtrk;
	  if(addHitPairToTrack(*iter, 0, 0, kmtrk_new))
	    {
	      _track_new.push_back(kmtrk_new);
	    }
#else
	  ///Both from right side
	  KalmanTrack kmtrk_pp = *kmtrk;
	  if(addHitPairToTrack(*iter, 1, 1, kmtrk_pp))
	    {
	      _track_new.push_back(kmtrk_pp);
	    }

	  ///Need to apply a continue command if driftDistance is 0, so that hodo. hits could be included
	  ///Won't be needed by now	
	  
	  ///Both from left side
	  KalmanTrack kmtrk_mm = *kmtrk;
	  if(addHitPairToTrack(*iter, -1, -1, kmtrk_mm))
	    {
	      _track_new.push_back(kmtrk_mm);
	    }

	  ///First left, second right
	  KalmanTrack kmtrk_pm = *kmtrk;
	  if(addHitPairToTrack(*iter, 1, -1, kmtrk_pm))
	    {
	      _track_new.push_back(kmtrk_pm);
	    }

	  ///First right, second left
	  KalmanTrack kmtrk_mp = *kmtrk;
	  if(addHitPairToTrack(*iter, -1, 1, kmtrk_mp))
	    {
	      _track_new.push_back(kmtrk_mp);
	    }
#endif
	  /*
	  ///One left, one right
	  ///For ST1V, both assumptions needs to be considered
	  ///For other stations, either +/- or -/+ is considered according to their relative position
	  KalmanTrack kmtrk_pm = *kmtrk;
	  if(fabs(_hits[iter->first].pos - _hits[iter->second].pos) < 1E-4) //This is for station 1V 
	    {
	      if(addHitPairToTrack(*iter, 1, -1, kmtrk_pm))
		{
		  _track_new.push_back(kmtrk_pm);
		}

	      KalmanTrack kmtrk_mp = *kmtrk;
              if(addHitPairToTrack(*iter, -1, 1, kmtrk_mp))
		{
		  _track_new.push_back(kmtrk_mp);
		}
	    }
	  else if(_hits[iter->first].pos < _hits[iter->second].pos && addHitPairToTrack(*iter, 1, -1, kmtrk_pm))
	    {
	      _track_new.push_back(kmtrk_pm);
	    }
	  else if(_hits[iter->first].pos > _hits[iter->second].pos && addHitPairToTrack(*iter, -1, 1, kmtrk_pm))
	    {
	      _track_new.push_back(kmtrk_pm);
	    }
	  */

	  ///Fill the last inserted one for evaluation	 
     	  fillEvaluation(iTrack, _track_new.back());
	}
      
      ++iTrack;
      ++kmtrk;
    }

  ///Insert the new candidates into the candidate list and sort
  _tracklist.insert(_tracklist.end(), _track_new.begin(), _track_new.end());
  _tracklist.sort();
 
#ifdef _DEBUG_ON
  Log(_track_new.size() << " new candidates are added. ");
  Log("Candidates after super detector " << detectorID);
  for(std::list<KalmanTrack>::iterator kmtrk = _tracklist.begin(); kmtrk != _tracklist.end(); ++kmtrk)
    {
      kmtrk->print();
    }
#endif
}

double KalmanFinder::optimizedWindow(int detectorID, KalmanTrack& _track)
{
  ///By default, the window is decided by error matrix
  double window = _track.getExpPosError();
  if(_track.getNHits() == 0)
    {
      return window;
    }

  ///At the entrance of each station, the window is enlarged
  if(detectorID == 9 || detectorID == 24 || detectorID == 18)
    {
      return 3.*window;
    }

  if(detectorID == 6)
    {
      return 2.5*window;
    }

  ///At station 1, the window is also momentum dependent
  if(detectorID <= 3)
    {
      double p = _track.getNodeList().front().getFiltered().get_mom();
      if(p < 20.)
	{
	  return 5.*window + 3.;
	}
      else if(p < 40.)
	{
	  return 4.*window + 3.;
	}
      else
	{
	  return 3.*window + 3.;
	}
    }

  return 2.*window;
}

bool KalmanFinder::addHitPairToTrack(SRawEvent::hit_pair _hitpair, int sign1, int sign2, KalmanTrack& _track)
{
  Hit hit1 = _hits[_hitpair.first];
  hit1.index = hit1.index*sign1*(hit1.driftDistance > 0 ? 1 : -1);
  hit1.driftDistance = sign1*fabs(hit1.driftDistance);

  Hit hit2 = _hits[_hitpair.second];
  hit2.index = hit2.index*sign2*(hit2.driftDistance > 0 ? 1 : -1);
  hit2.driftDistance = sign2*fabs(hit2.driftDistance);

  ///Effective checksum, the local slope shouldn't be too much different with the expected value
  int detectorID = hit1.detectorID;
  double slop_exp = _track.getExpLocalSlop();
  double slop_err = _track.getExpLcSlopErr();
  double slop_local = (hit1.pos + hit1.driftDistance - hit2.pos - hit2.driftDistance)/(p_geomSvc->getPlanePosition(hit1.detectorID) - p_geomSvc->getPlanePosition(hit2.detectorID));

#ifdef _DEBUG_ON
  Log("Testing possibility: " << sign1 << " == " << sign2);
  Log("Slop_exp: " << slop_exp << ", slop_err: " << slop_err << ", slop_local: " << slop_local);
#endif

  if(detectorID > 6 && fabs(slop_local - slop_exp) > 4*slop_err) return false;
  if(fabs(slop_local) > 0.25) return false;

  if(!_track.addHit(hit1)) return false;
  if(!_track.propagateTo(hit2.detectorID)) return false;
  if(!_track.addHit(hit2)) return false;

#ifdef _DEBUG_ON
  Log("Added!");
#endif

  return true;
}

bool KalmanFinder::refineTrack(KalmanTrack& _track)
{
  ///Note, for sub tracks in station 2&3, refine is not necessarily successful.
  // if not, just skip the refine.
  // For track candidates with st1 hits, refine is required to be successful

  if(_track.getNodeList().empty()) return false;

  KalmanFitter *fitter = new KalmanFitter();
  fitter->setControlParameter(50, 0.001);
  
  if(fitter->processOneTrack(_track) == 0)
    {
      delete fitter;
      return false;
    }
  fitter->updateTrack(_track);
 
  delete fitter;
  return true;
}


void KalmanFinder::reduceTrackList(int detectorID, std::list<KalmanTrack>& _tracklist)
{
  if((detectorID/2)*2 == detectorID) return;

  int factor = (((24 - detectorID) % 6) + 1)/2;
  unsigned int max_size = 30/factor;

#ifdef _DEBUG_ON
  Log("Before: " << _tracklist.size()); 
#endif  
  
  for(std::list<KalmanTrack>::iterator iter = _tracklist.begin(); iter != _tracklist.end(); )
    {
      if(!acceptTrack(detectorID, *iter))
	{
	  //Log("Removed!");
	  iter = _tracklist.erase(iter);
	  continue;
	}
      else
	{
	  //Log("Kept!");
	  ++iter;
	}
    }

#ifdef _DEBUG_ON
  Log("After : " << _tracklist.size()); 
#endif

  _tracklist.sort();
  if(_tracklist.size() > max_size + 1)
    {
      _tracklist.resize(max_size);
    }
}

bool KalmanFinder::acceptTrack(int detectorID, KalmanTrack& _track)
{
#ifdef _DEBUG_ON
  Log(detectorID);
  _track.print();
#endif

  ///Basic quality checks
  double _chisq = _track.getChisq();
  if(_chisq < 0. || _chisq > 500.) return false;
 
  int nHits = _track.getNHits();
  if(nHits > 3)
    {
      double _p_up = _track.getMomentumUpstream();
      if(_p_up < 15. || _p_up > 100.) return false;
    }
 
  ///Check if the sign has been fliped in all the nodes
  if(_track.getNHits() > 2)
    {
      int charge_front = _track.getNodeList().front().getFiltered().get_charge();
      int charge_back = _track.getNodeList().back().getFiltered().get_charge();
      if(charge_front*charge_back < 0)
        {
#ifdef _DEBUG_ON
          Log(detectorID << ": Sign is flipped during the track finding!");
#endif
	  //_track.getNodeList().front().getFiltered().flip_charge();
	  return false;
        }
    }

  ///Test the hodoscope maskinig, not used for now
  int stationID = int((detectorID-1)/6) + 1;
  //if(!hodoMask(stationID, _track)) return false;

  ///Test the number of hits after each station
  if(detectorID % 6 == 1)
    {
      int nHits = _track.getNHitsInStation(stationID);
      if(nHits < 2) return false;
    }
  
  ///Test the difference before and after the KMag swim
  if(detectorID == 5)
    {
      double p_st2 = _track.getMomentumInStation(2);
      double p_st1 = _track.getMomentumInStation(1);
 
      int charge_assigned = _track.getAssignedCharge();
      int charge_kick = _track.getKickCharge();

      if(charge_assigned*charge_kick < 0) return false;
      if(fabs(p_st1 - p_st2) > 20.) return false;

      if(!refineTrack(_track)) return false;
    }

  return true;
}

int KalmanFinder::getNHodoHits(int stationID, KalmanTrack& _track)
{
  if(stationID == 3 || stationID == 4)
    {
      stationID = 2;
    }
  else
    {
      stationID = stationID - 1;
    }

  unsigned int nMask = 0;
  for(unsigned int i = 0; i < detectors_mask[stationID].size(); i++)
    {
      double z_exp = p_geomSvc->getPlanePosition(detectors_mask[stationID][i]);
      double x_exp, y_exp, x_err, y_err;
      Node *_node = _track.getNearestNodePtr(z_exp);
      _track.getExpPositionFast(z_exp, x_exp, y_exp, _node);
      _track.getExpPosErrorFast(z_exp, x_err, y_err, _node);

      for(std::list<int>::iterator iter = hits_mask[stationID][i].begin(); iter != hits_mask[stationID][i].end(); ++iter)
	{
	  double x_min, x_max, y_min, y_max;
	  p_geomSvc->get2DBoxSize(_hits[*iter].detectorID, _hits[*iter].elementID, x_min, x_max, y_min, y_max);
	  x_min -= x_err; x_max += x_err;
	  y_min -= y_err; y_max += y_err;

	  if(x_exp > x_min && x_exp < x_max && y_exp > y_min && y_exp < y_max)
	    {
	      nMask++;
	      if((i & 1) == 0)
		{
		  i++;
		}
	      break;
	    }
	}
    }

  return nMask;
}

KalmanTrack KalmanFinder::associateSingles(KalmanTrack _track)
{
  std::vector<int> detectorIDs = _track.getMissedDetectorIDs();
  std::list<Node> nodes_added;
  for(std::vector<int>::iterator iter = detectorIDs.begin(); iter != detectorIDs.end(); ++iter)
    {
      double x, y, dx, dy;
      double z = p_geomSvc->getPlanePosition(*iter);
      Node *node_nearest = _track.getNearestNodePtr(z);
      _track.getExpPositionFast(z, x, y, node_nearest);
      _track.getExpPosErrorFast(z, dx, dy, node_nearest);

      double pos_exp = p_geomSvc->getUinStereoPlane(*iter, x, y);
      double trk_err = fabs(p_geomSvc->getCostheta(*iter)*dx) + fabs(p_geomSvc->getSintheta(*iter)*dy);
      double spacing = p_geomSvc->getPlaneSpacing(*iter);
      double window = sqrt(trk_err*trk_err + spacing*spacing);

      std::list<int> hitlist = event->getHitsIndexInDetector(*iter, pos_exp, window);
      if(hitlist.empty()) continue;

      double delta_min = 1E6;
      int sign_min = 0;
      int index_min = -1;
      for(std::list<int>::iterator jter = hitlist.begin(); jter != hitlist.end(); ++jter)
	{
	  double delta_p = fabs(_hits[*jter].pos + fabs(_hits[*jter].driftDistance) - pos_exp);
	  double delta_m = fabs(_hits[*jter].pos - fabs(_hits[*jter].driftDistance) - pos_exp);
	  double delta = delta_p < delta_m ? delta_p : delta_m;
	  int sign = delta_p < delta_m ? 1 : -1;

	  if(delta_min > delta)
	    {
	      delta_min = delta;
	      sign_min = sign;
	      index_min = *jter;
	    }
	}

      if(index_min < 0) continue;

      double resol = p_geomSvc->getPlaneResolution(*iter);
      double delta_exp = 1.5*sqrt(trk_err*trk_err + resol*resol); //3 sigma window
      if(delta_min > delta_exp) continue;
      if(delta_min > 0.5*spacing) continue;

      Hit _hit = _hits[index_min];
      _hit.index = _hit.index*sign_min*(_hit.driftDistance > 0 ? 1 : -1);
      _hit.driftDistance = sign_min*fabs(_hit.driftDistance);

      Node _node(_hit);
      KalmanFilter *_kmfit = KalmanFilter::instance();
      _kmfit->setCurrTrkpar(*node_nearest);
      _kmfit->fit_node(_node);

      nodes_added.push_back(_node);
    }

  for(std::list<Node>::iterator iter = nodes_added.begin(); iter != nodes_added.end(); ++iter)
    {
      _track.getHitIndexList().push_back(iter->getHit().index);
      _track.getNodeList().push_back(*iter);
    }

  _track.getNodeList().sort();
  return _track;
}

KalmanTrack KalmanFinder::getBestCandidate()
{
  for(std::list<KalmanTrack>::iterator iter = _track_candidates.begin(); iter != _track_candidates.end(); ++iter)
    {
      if(!refineTrack(*iter)) continue;
 
      KalmanTrack track_final = associateSingles(*iter);
      if(track_final.getNHits() > iter->getNHits()) 
	{
	  if(refineTrack(track_final) && (track_final.getQuality() > iter->getQuality())) 
	    {
	      return track_final;
	    }
	}

      return *iter;
    }
	  
  KalmanTrack nll;
  nll.update();

  return nll;
}

void KalmanFinder::reduceTrackList(std::list<KalmanTrack>& _trklist_source, std::list<KalmanTrack>& _trklist_target)
{
  ///If two tracks share a lot of common hits, the one with better quality is kept.
  for(std::list<KalmanTrack>::iterator iter = _trklist_source.begin(); iter != _trklist_source.end(); )
    {
      if(!acceptTrack(1, *iter))
	{
	  iter = _trklist_source.erase(iter);
	  continue;
	}

      ++iter;
    }
  _trklist_source.sort();

  while(!_trklist_source.empty())
    {
      _trklist_target.push_back(_trklist_source.front());
      _trklist_source.pop_front();

      for(std::list<KalmanTrack>::iterator iter = _trklist_source.begin(); iter != _trklist_source.end(); )
	{
	  if(iter->similarity(_trklist_target.back()))
	    {
	      iter = _trklist_source.erase(iter);
	      continue;
	    }

	  ++iter;
	}
    }
}

void KalmanFinder::printResults(std::string outputFileName, std::list<KalmanTrack>& _tracks)
{
  double pos[2000], z[2000];
  double err[2000], dz[2000];

  std::vector<int> detectorIDs = p_geomSvc->getDetectorIDs("X");
  std::list<int> hitIDs = event->getHitsIndexInDetectors(detectorIDs);
  int nHits = 0;
  for(std::list<int>::iterator iter = hitIDs.begin(); iter != hitIDs.end(); ++iter)
    {
      int detectorID = _hits[*iter].detectorID;
      z[nHits] = p_geomSvc->getPlanePosition(detectorID);
      dz[nHits] = 0.;
      pos[nHits] = _hits[*iter].pos;
      err[nHits] = 0.5*p_geomSvc->getPlaneSpacing(detectorID);

      ++nHits;
    }

  TGraphErrors _hitvis(nHits, z, pos, dz, err);

  ///Tracks
  std::vector<TGraph> _trackvis;
  _trackvis.clear();
  for(std::list<KalmanTrack>::iterator iter = _tracks.begin(); iter != _tracks.end(); ++iter)
    {
      TGraph gr = iter->getXZProjection();
     _trackvis.push_back(gr);
    }

  TCanvas c1;
  c1.cd(); _hitvis.Draw("AP");
  for(std::vector<TGraph>::iterator iter = _trackvis.begin(); iter != _trackvis.end(); ++iter)
    {
      iter->SetMarkerColor(2);
      iter->SetMarkerStyle(2);
 
      iter->Draw("LP same");
    }

  std::stringstream filename(outputFileName);
  filename << event->getRunID() << "_" << event->getEventID() << ".eps";

  c1.SaveAs(filename.str().c_str());
}

void KalmanFinder::printResults(std::string outputFileName, KalmanTrack& _track)
{
  std::list<KalmanTrack> _tracks;
  _tracks.push_back(_track);

  printResults(outputFileName, _tracks);
}

void KalmanFinder::initEvaluation(std::string outputFileName)
{
  _eval_file = new TFile(outputFileName.c_str(), "recreate");
  
  _res_tree = new TTree("save", "save");

  _res_tree->Branch("runID", &_runID, "runID/I");
  _res_tree->Branch("eventID", &_eventID, "eventID/I");
  _res_tree->Branch("trackID", &_trackID, "trackID/I");
  _res_tree->Branch("detectorID", &_detectorID, "detectorID/I");
  _res_tree->Branch("res_predicted", &_res_predicted, "res_predicted/D");
  _res_tree->Branch("mom_predicted", &_mom_predicted, "mom_predicted/D");
  _res_tree->Branch("res_cov_predicted", &_res_cov_predicted, "res_cov_predicted/D");
}

void KalmanFinder::finishEvaluation()
{
  _eval_file->cd();
  _res_tree->Write();
  _eval_file->Close();
}

void KalmanFinder::fillEvaluation(int trackID, KalmanTrack& _track)
{
  if(_res_tree != NULL)
    {
      _runID = event->getRunID();
      _eventID = event->getEventID();
      _trackID = trackID;

      std::list<Node>::iterator _node = _track.getNodeList().begin();
      if(_node->isPredictionDone())
	{
	  _detectorID = _node->getHit().detectorID;
	  _mom_predicted = _node->getPredicted().get_mom();
          _res_predicted = _node->getPredictedResidual()[0][0]; 
	  _res_cov_predicted = sqrt(_node->getPredictedResidualCov()[0][0]);
          _res_tree->Fill();
	}

      ++_node;
      if(_node->isPredictionDone())
	{
          _detectorID = _node->getHit().detectorID;
	  _mom_predicted = _node->getPredicted().get_mom();
          _res_predicted = _node->getPredictedResidual()[0][0];
	  _res_cov_predicted = sqrt(_node->getPredictedResidualCov()[0][0]);
          _res_tree->Fill();
	}
    }
}
