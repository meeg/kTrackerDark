/*
SRawEvent.cxx

Implimention of the class SRawEvent

Author: Kun Liu, liuk@fnal.gov
Created: 10-24-2011
*/

#include <iostream>
#include <cmath>

#include <TRandom.h>
#include <TMath.h>
#include <TString.h>

#include "SRawEvent.h"

ClassImp(Hit)
ClassImp(SRawEvent)
ClassImp(SRawMCEvent)

bool Hit::operator<(const Hit& elem) const
{
  if(detectorID < elem.detectorID)
    {
      return true;
    }
  else if(detectorID > elem.detectorID)
    {
      return false;
    }
	  
  if(elementID < elem.elementID)
    {
      return true;
    }
  else if(elementID > elem.elementID)
    {
      return false;
    }

  if(tdcTime > elem.tdcTime)
    {
      return true;
    }
  else
    {
      return false;
    }
}

bool Hit::operator==(const Hit& elem) const
{
  if(detectorID == elem.detectorID && elementID == elem.elementID)
    {
      return true;
    }

  if(detectorID == elem.detectorID && fabs(pos - elem.pos) < 1E-3)
    {
      return true;
    }

  return false;
}

SRawEvent::SRawEvent()
{
  fAllHits.clear();
  for(Int_t i = 0; i < nChamberPlanes+nHodoPlanes+nPropPlanes+1; i++)
    {
      fNHits[i] = 0;
    }

  fRunID = -1;
  fEventID = -1;
  fSpillID = -1;

  fTriggerBits = 0;
  fTriggerHits.clear();
}

SRawEvent::~SRawEvent()
{
}

void SRawEvent::setEventInfo(Int_t runID, Int_t spillID, Int_t eventID)
{
  fRunID = runID;
  fEventID = eventID;
  fSpillID = spillID;
}

void SRawEvent::insertHit(Hit h)
{
  fAllHits.push_back(h);

  fNHits[0]++;
  fNHits[h.detectorID]++;
}

Int_t SRawEvent::findHit(Int_t detectorID, Int_t elementID)
{
  if(detectorID < 1 || detectorID > 48) return -1;
  if(elementID < 0) return -1;

  Hit h_dummy;
  h_dummy.detectorID = detectorID;
  h_dummy.elementID = elementID;

  /*
  This method produces problems in case of duplicate channels and thus people need to be cautious;
  It's okay here for two reasons:
     1. inTime is required when searching for trigger roads;
     2. hodoscope hit doesn't need tdcTime information as long as it's in-time;
  */
  Int_t idx_start;
  Int_t idx_end;
  if(detectorID <= 24)
    {
      idx_start = 0;
      idx_end = getNChamberHitsAll() - 1;
    }
  else if(detectorID <= 40)
    {
      idx_start = getNChamberHitsAll();
      idx_end = idx_start + getNHodoHitsAll();
    }
  else
    {
      idx_start = getNChamberHitsAll() + getNHodoHitsAll();
      idx_end = fNHits[0] - 1;
    }

  while(idx_start <= idx_end)
    {
      Int_t idx_mid = Int_t((idx_start + idx_end)/2);
      if(fAllHits[idx_mid] == h_dummy)
	{
	  return idx_mid;
	}
      else if(fAllHits[idx_mid] < h_dummy)
	{
	  idx_start = idx_mid + 1;
	}
      else
	{
	  idx_end = idx_mid - 1;
	}
    }

  return -1;
}

Hit SRawEvent::getHit(Int_t detectorID, Int_t elementID)
{
  Int_t hitID = findHit(detectorID, elementID);
  if(hitID >= 0) return getHit(hitID);

  Hit dummy;
  dummy.index = -1;
  dummy.detectorID = -1;
  dummy.elementID = -1;
  return dummy;
}

std::list<Int_t> SRawEvent::getHitsIndexInDetector(Int_t detectorID)
{
  std::list<Int_t> hit_list;
  hit_list.clear();

  for(Int_t i = 0; i < fNHits[0]; i++)
    {
      if(fAllHits[i].detectorID != detectorID) continue;

      hit_list.push_back(i);
    }

  return hit_list;
}


std::list<Int_t> SRawEvent::getHitsIndexInDetector(Int_t detectorID, Double_t x_exp, Double_t win)
{
  std::list<Int_t> hit_list;
  hit_list.clear();

  for(Int_t i = 0; i < fNHits[0]; i++)
    {
      if(fAllHits[i].detectorID != detectorID) continue;
      if(fabs(fAllHits[i].pos - x_exp) > win) continue;

      hit_list.push_back(i);
    }

  return hit_list;
}

std::list<Int_t> SRawEvent::getHitsIndexInSuperDetector(Int_t detectorID)
{
  std::list<Int_t> hit_list;
  hit_list.clear();

  for(Int_t i = 0; i < fNHits[0]; i++)
    {
      if((fAllHits[i].detectorID != 2*detectorID) && (fAllHits[i].detectorID != 2*detectorID-1)) continue;

      hit_list.push_back(i);
    }

  return hit_list;
}

std::list<Int_t> SRawEvent::getHitsIndexInDetectors(std::vector<Int_t>& detectorIDs)
{
  std::list<Int_t> hit_list;
  hit_list.clear();

  UInt_t nDetectors = detectorIDs.size();
  for(Int_t i = 0; i < fNHits[0]; i++)
    {
      for(UInt_t j = 0; j < nDetectors; j++)
	{
	  if(fAllHits[i].detectorID == detectorIDs[j])
	    {
	      hit_list.push_back(i);
	      break;
	    }
	}
    }

  return hit_list;
}

std::list<SRawEvent::hit_pair> SRawEvent::getHitPairsInSuperDetector(Int_t detectorID)
{
  std::list<SRawEvent::hit_pair> _hitpairs;
  std::list<int> _hitlist1 = getHitsIndexInDetector(2*detectorID);
  std::list<int> _hitlist2 = getHitsIndexInDetector(2*detectorID - 1);

  for(std::list<int>::iterator iter = _hitlist1.begin(); iter != _hitlist1.end(); ++iter)
    {
      for(std::list<int>::iterator jter = _hitlist2.begin(); jter != _hitlist2.end(); ++jter)
	{
	  if(abs(fAllHits[*iter].elementID - fAllHits[*jter].elementID) > 1) continue;
	  //if((fAllHits[*iter].elementID == fAllHits[*jter].elementID) || (fAllHits[*iter].elementID + 1 == fAllHits[*jter].elementID)) continue;
	  _hitpairs.push_back(std::make_pair(*iter, *jter));
	}
    }

  return _hitpairs;
}

std::list<SRawEvent::hit_pair> SRawEvent::getPartialHitPairsInSuperDetector(Int_t detectorID)
{
  std::list<SRawEvent::hit_pair> _hitpairs;
  std::list<int> _hitlist1 = getHitsIndexInDetector(2*detectorID);
  std::list<int> _hitlist2 = getHitsIndexInDetector(2*detectorID - 1);

  std::vector<int> _hitflag1(_hitlist1.size(), -1);
  std::vector<int> _hitflag2(_hitlist2.size(), -1);

  //Temp solutions here
  double spacing[13] = {0., 0.40, 0.40, 0.40, 1.3, 1.3, 1.3, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2};

  int index1 = -1;
  int index2 = -1;
  for(std::list<int>::iterator iter = _hitlist1.begin(); iter != _hitlist1.end(); ++iter)
    {
      index1++;
      index2 = -1;
      for(std::list<int>::iterator jter = _hitlist2.begin(); jter != _hitlist2.end(); ++jter)
	{
	  index2++;
	  if(fabs(fAllHits[*iter].pos - fAllHits[*jter].pos) > spacing[detectorID]) continue;
	  _hitpairs.push_back(std::make_pair(*iter, *jter));

	  _hitflag1[index1] = 1;
	  _hitflag2[index2] = 1;
	}
    }

  index1 = 0;
  for(std::list<int>::iterator iter = _hitlist1.begin(); iter != _hitlist1.end(); ++iter)
    {
      if(_hitflag1[index1] < 0) _hitpairs.push_back(std::make_pair(*iter, -1));
      ++index1;
    }

  index2 = 0;
  for(std::list<int>::iterator iter = _hitlist2.begin(); iter != _hitlist2.end(); ++iter)
    {
      if(_hitflag2[index2] < 0) _hitpairs.push_back(std::make_pair(*iter, -1));
      ++index2;
    }

  return _hitpairs;
}

std::list<SRawEvent::hit_pair> SRawEvent::getPartialHitPairsInSuperDetector(Int_t detectorID, Double_t x_exp, Double_t win)
{
  std::list<SRawEvent::hit_pair> _hitpairs;
  std::list<int> _hitlist1 = getHitsIndexInDetector(2*detectorID, x_exp, win);
  std::list<int> _hitlist2 = getHitsIndexInDetector(2*detectorID - 1, x_exp, win+3);

  std::vector<int> _hitflag1(_hitlist1.size(), -1);
  std::vector<int> _hitflag2(_hitlist2.size(), -1);

  //Temp solutions here
  double spacing[13] = {0., 0.40, 0.40, 0.40, 1.3, 1.3, 1.3, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2};

  int index1 = -1;
  int index2 = -1;
  for(std::list<int>::iterator iter = _hitlist1.begin(); iter != _hitlist1.end(); ++iter)
    {
      index1++;
      index2 = -1;
      for(std::list<int>::iterator jter = _hitlist2.begin(); jter != _hitlist2.end(); ++jter)
	{
	  index2++;
	  if(fabs(fAllHits[*iter].pos - fAllHits[*jter].pos) > spacing[detectorID]) continue;
	  _hitpairs.push_back(std::make_pair(*iter, *jter));

	  _hitflag1[index1] = 1;
	  _hitflag2[index2] = 1;
	}
    }

  index1 = 0;
  for(std::list<int>::iterator iter = _hitlist1.begin(); iter != _hitlist1.end(); ++iter)
    {
      if(_hitflag1[index1] < 0) _hitpairs.push_back(std::make_pair(*iter, -1));
      ++index1;
    }

  index2 = 0;
  for(std::list<int>::iterator iter = _hitlist2.begin(); iter != _hitlist2.end(); ++iter)
    {
      if(_hitflag2[index2] < 0) _hitpairs.push_back(std::make_pair(*iter, -1));
      ++index2;
    }

  return _hitpairs;
}

std::list<SRawEvent::hit_pair> SRawEvent::getHitPairsInSuperDetector(Int_t detectorID, Double_t x_exp, Double_t win)
{
  std::list<SRawEvent::hit_pair> _hitpairs;
  std::list<int> _hitlist1 = getHitsIndexInDetector(2*detectorID, x_exp, win);
  std::list<int> _hitlist2 = getHitsIndexInDetector(2*detectorID - 1, x_exp, win + 3.);

  for(std::list<int>::iterator iter = _hitlist1.begin(); iter != _hitlist1.end(); ++iter)
    {
      for(std::list<int>::iterator jter = _hitlist2.begin(); jter != _hitlist2.end(); ++jter)
	{
	  if(abs(fAllHits[*iter].elementID - fAllHits[*jter].elementID) > 1) continue;
	  //if((fAllHits[*iter].elementID == fAllHits[*jter].elementID) || (fAllHits[*iter].elementID + 1 == fAllHits[*jter].elementID)) continue;

	  _hitpairs.push_back(std::make_pair(*iter, *jter));
	}
    }

  return _hitpairs;
}


std::list<Int_t> SRawEvent::getAdjacentHitsIndex(Hit& _hit)
{
  std::list<Int_t> hit_list;
  hit_list.clear();

  Int_t detectorID = _hit.detectorID;
  Int_t detectorID_adj;
  if((detectorID/2)*2 == detectorID)
    {
      detectorID_adj = detectorID - 1;
    }
  else
    {
      detectorID_adj = detectorID + 1;
    }

  for(Int_t i = 0; i < fNHits[0]; i++)
    {
      if(fAllHits[i].detectorID == detectorID_adj && abs(fAllHits[i].elementID - _hit.elementID) <= 1)
	{
	  hit_list.push_back(i);
	}
    }

  return hit_list;
}



Int_t SRawEvent::getNChamberHitsAll()
{
  Int_t nHits = 0;
  for(Int_t i = 1; i <= nChamberPlanes; i++)
    {
      nHits += fNHits[i];
    }

  return nHits;
}

Int_t SRawEvent::getNHodoHitsAll()
{
  Int_t nHits = 0;
  for(Int_t i = nChamberPlanes+1; i <= nChamberPlanes+nHodoPlanes; i++)
    {
      nHits += fNHits[i];
    }

  return nHits;
}

Int_t SRawEvent::getNPropHitsAll()
{
  Int_t nHits = 0;
  for(Int_t i = nChamberPlanes+nHodoPlanes+1; i <= nChamberPlanes+nHodoPlanes+nPropPlanes; i++)
    {
      nHits += fNHits[i];
    }

  return nHits;
}

Int_t SRawEvent::getNHitsInDetectors(std::vector<Int_t>& detectorIDs)
{
  Int_t nHits = 0;
  UInt_t nDetectors = detectorIDs.size();
  for(UInt_t i = 0; i < nDetectors; i++)
    {
      for(Int_t j = 0; j <= nChamberPlanes+nHodoPlanes+nPropPlanes; j++)
	{
	  if(detectorIDs[i] == j)
	    {
	      nHits += fNHits[j];
	      break;
	    }
	}
    }

  return nHits;
}

Int_t SRawEvent::getNHitsInD1()
{
  Int_t nHits = 0;
  for(Int_t i = 1; i <= 6; i++)
    {
      nHits += fNHits[i];
    }

  return nHits;
}

Int_t SRawEvent::getNHitsInD2()
{
  Int_t nHits = 0;
  for(Int_t i = 7; i <= 12; i++)
    {
      nHits += fNHits[i];
    }

  return nHits;
}

Int_t SRawEvent::getNHitsInD3p()
{
  Int_t nHits = 0;
  for(Int_t i = 13; i <= 18; i++)
    {
      nHits += fNHits[i];
    }

  return nHits;
}

Int_t SRawEvent::getNHitsInD3m()
{
  Int_t nHits = 0;
  for(Int_t i = 19; i <= 24; i++)
    {
      nHits += fNHits[i];
    }

  return nHits;
}

void SRawEvent::reIndex(std::string option)
{
  bool _afterhit = false;
  bool _hodomask = false;
  bool _outoftime = false;
  bool _nonchamber = false;
  bool _decluster = false;
  bool _mergehodo = false;
  bool _triggermask = false;

  TString option_lower(option.c_str());
  option_lower.ToLower();
  if(option_lower.Contains("a")) _afterhit = true;
  if(option_lower.Contains("h")) _hodomask = true;
  if(option_lower.Contains("o")) _outoftime = true;
  if(option_lower.Contains("n")) _nonchamber = true;
  if(option_lower.Contains("c")) _decluster = true;
  if(option_lower.Contains("u")) _mergehodo = true;
  if(option_lower.Contains("t")) _triggermask = true;

  ///Dump the vector into a list and do the reduction
  std::list<Hit> hitlist_temp;
  hitlist_temp.clear();
  for(std::vector<Hit>::iterator iter = fAllHits.begin(); iter != fAllHits.end(); ++iter)
    {
      if(_outoftime && iter->inTime == 0) continue;
      if(_hodomask && iter->hodoMask == 0) continue;
      if(_nonchamber && iter->detectorID > 24) continue;
      if(_triggermask && iter->detectorID >= 24 && iter->detectorID <= 40 && iter->inTime != 2) continue;

      hitlist_temp.push_back(*iter);
    }

  if(_mergehodo)
    {
      for(std::vector<Hit>::iterator iter = fTriggerHits.begin(); iter != fTriggerHits.end(); ++iter) hitlist_temp.push_back(*iter);
    }

  ///Remove after hits
  hitlist_temp.sort();
  if(_afterhit)
    {
      hitlist_temp.unique();
    }
 
  if(_decluster)
    {
      deClusterize(hitlist_temp);
    }

  fAllHits.clear();
  fAllHits.assign(hitlist_temp.begin(), hitlist_temp.end());

  ///Reset the number of hits on each plane
  for(Int_t i = 0; i < nChamberPlanes+nHodoPlanes+nPropPlanes+1; i++)
    {
      fNHits[i] = 0;
    }

  for(UInt_t i = 0; i < fAllHits.size(); i++)
    {
      ++fNHits[fAllHits[i].detectorID];
    }

  fNHits[0] = fAllHits.size();
}

void SRawEvent::deClusterize(std::list<Hit>& hits)
{
  std::vector<std::list<Hit>::iterator> cluster;
  cluster.clear();
  for(std::list<Hit>::iterator hit = hits.begin(); hit != hits.end(); ++hit)
    {
      //if we already reached the hodo part, stop
      if(hit->detectorID > 24) break;

      if(cluster.size() == 0)
	{
	  cluster.push_back(hit);
	}
      else
	{
	  if(hit->detectorID != cluster.back()->detectorID)
	    {
	      processCluster(hits, cluster);
	      cluster.push_back(hit);
	    }
	  else if(hit->elementID - cluster.back()->elementID > 1)
	    {
	      processCluster(hits, cluster);
	      cluster.push_back(hit);      
	    }
	  else
	    {
	      cluster.push_back(hit);
	    }
	}
    }
}

void SRawEvent::processCluster(std::list<Hit>& hits, std::vector<std::list<Hit>::iterator>& cluster)
{
  //size-2 clusters, retain the hit with smaller driftDistance
  if(cluster.size() == 2)
    {
      double w_max = 0.9*0.5*(cluster.back()->pos - cluster.front()->pos);
      double w_min = 0.4*0.5*(cluster.back()->pos - cluster.front()->pos);
      if((cluster.front()->driftDistance > w_max && cluster.back()->driftDistance > w_min) || (cluster.front()->driftDistance > w_min && cluster.back()->driftDistance > w_max))
	{
	  cluster.front()->driftDistance > cluster.front()->driftDistance ? hits.erase(cluster.front()) : hits.erase(cluster.back());
	}
    }

  //size-larger-than-3, discard entirely
  if(cluster.size() >= 3)
    {
      double dt_mean = 0.;
      for(unsigned int i = 1; i < cluster.size(); ++i)
	{
	  dt_mean += fabs(cluster[i]->tdcTime - cluster[i-1]->tdcTime);
	}
      dt_mean = dt_mean/(cluster.size() - 1);

      if(dt_mean < 10.)
	{
	  //electric noise, discard them all
	  for(unsigned int i = 0; i < cluster.size(); ++i)
	    {
	      hits.erase(cluster[i]);
	    }
	}
      else
	{
	  //delta ray, keep the first and last
	  for(unsigned int i = 1; i < cluster.size() - 1; ++i)
	    {
	      hits.erase(cluster[i]);
	    }
	}
    }

  cluster.clear();
}

void SRawEvent::mixEvent(SRawEvent *event, int nBkgHits)
{
  event->reIndex("oah");

  for(Int_t i = 0; i < nChamberPlanes+nHodoPlanes+nPropPlanes+1; i++)
    {
      fNHits[i] = 0;
    }

  std::vector<Hit> hits_mix = event->getAllHits();
  for(std::vector<Hit>::iterator iter = fAllHits.begin(); iter != fAllHits.end(); ++iter)
    {
      for(std::vector<Hit>::iterator jter = hits_mix.begin(); jter != hits_mix.end(); ++jter)
	{
	  if(*iter == *jter)
	    {
	      hits_mix.erase(jter);
	      break;
	    }
	}
    }

  if(nBkgHits > 0)
    {
      TRandom rdn;
      double ratio = double(nBkgHits)/event->getNChamberHitsAll();
      for(std::vector<Hit>::iterator iter = hits_mix.begin(); iter != hits_mix.end(); )
	{
	  if(rdn.Rndm() > ratio)
	    {
	      iter = hits_mix.erase(iter);
	    }
	  else
	    {
	      ++iter;
	    }
	}
    }

  fAllHits.insert(fAllHits.end(), hits_mix.begin(), hits_mix.end());
  for(UInt_t i = 0; i < fAllHits.size(); i++)
    {
      fAllHits[i].index = i;
      fNHits[fAllHits[i].detectorID]++;
    }

  fNHits[0] = fAllHits.size();
}

void SRawEvent::setEventInfo(SRawEvent* event)
{
  //Set runID, eventID, spillID
  setEventInfo(event->getRunID(), event->getSpillID(), event->getEventID());

  //Set trigger bits
  setTriggerBits(event->getTriggerBits());

  //Set target position
  setTargetPos(event->getTargetPos());

  //Set bean info
  setTurnID(event->getTurnID());
  setRFID(event->getRFID());
  setIntensity(event->getIntensityAll());

  //Set the trigger emu info
  setTriggerEmu(event->isEmuTriggered());
  setNRoads(event->getNRoads());
}

void SRawEvent::clear()
{
  fAllHits.clear();
  for(Int_t i = 0; i < nChamberPlanes+1; i++)
    {
      fNHits[i] = 0;
    }

  fRunID = -1;
  fSpillID = -1;
  fEventID = -1;
  fTriggerBits = 0;
  fTriggerHits.clear();
}

void SRawEvent::setTriggerBits(Int_t triggers[])
{
  for(int i = 0; i < 10; ++i)
    {
      if(triggers[i] == 0) continue;
      fTriggerBits |= triggerBit(i+1);
    }
}

void SRawEvent::print()
{
  std::cout << "RunID: " << fRunID << ", EventID: " << fEventID << "===============" << std::endl;
  for(Int_t i = 1; i <= nChamberPlanes; i++)
    {
      std::cout << "Layer " << i << " has " << fNHits[i] << " hits." << std::endl;
    }
  std::cout << "===================================================================" << std::endl;

  return;
  for(std::vector<Hit>::iterator iter = fAllHits.begin(); iter != fAllHits.end(); ++iter)
    {
      iter->print();
    }
  std::cout << "===================================================================" << std::endl;
}
