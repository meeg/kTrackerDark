/*
SRawEvent.cxx

Implimention of the class SRawEvent

Author: Kun Liu, liuk@fnal.gov
Created: 10-24-2011
*/

#include <iostream>
#include <cmath>
#include <TRandom.h>

#include "SRawEvent.h"

ClassImp(Hit)
ClassImp(SRawEvent)
ClassImp(SRawMCEvent)

bool Hit::operator<(const Hit& elem) const
{
  if(inTime > elem.inTime)
    {
      return true;
    }
  else if(inTime < elem.inTime)
    {
      return false;
    }

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
  else
    {
      return false;
    }

  if(driftTime < elem.driftTime)
    {
      return true;
    }
  else
    {
      return false;
    }
}

bool Hit::sameChannel(const Hit& elem1, const Hit& elem2)
{
  if(elem1.detectorID == elem2.detectorID && elem1.elementID == elem2.elementID)
    {
      return true;
    }

  if(elem1.detectorID == elem2.detectorID && fabs(elem1.pos - elem2.pos) < 1E-3)
    {
      return true;
    }

  return false;
}

SRawEvent::SRawEvent()
{
  fAllHits.clear();
  for(Int_t i = 0; i < nChamberPlanes+1; i++)
    {
      fNHits[i] = 0;
    }

  fRunID = -1;
  fEventID = -1;
  fSpillID = -1;
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

  for(Int_t i = 0; i < fNHits[0]; i++)
    {
      if(fAllHits[i].detectorID == detectorID && fAllHits[i].elementID == elementID)
	{
	  return i;
	}
    }

  return -1;
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
  double spacing[13] = {0., 0.38, 0.38, 0.38, 1.22, 1.22, 1.22, 1.2, 1.2, 1.2, 0.6, 0.6, 0.6};

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
  double spacing[13] = {0., 0.38, 0.38, 0.38, 1.22, 1.22, 1.22, 1.2, 1.2, 1.2, 0.6, 0.6, 0.6};

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

  std::transform(option.begin(), option.end(), option.begin(), tolower);
  if(option.find("a") != std::string::npos) _afterhit = true;
  if(option.find("h") != std::string::npos) _hodomask = true;
  if(option.find("o") != std::string::npos) _outoftime = true;
  if(option.find("n") != std::string::npos) _nonchamber = true;

  ///Dump the vector into a list and do the reduction
  std::list<Hit> hitlist_temp;
  hitlist_temp.clear();
  for(UInt_t i = 0; i < fAllHits.size(); i++)
    {
      if(_outoftime && fAllHits[i].inTime == 0) continue;
      if(_hodomask && fAllHits[i].hodoMask == 0) continue;
      if(_nonchamber && fAllHits[i].detectorID > 24) continue;

      hitlist_temp.push_back(fAllHits[i]);
    }

  ///Remove after hits
  hitlist_temp.sort();
  if(_afterhit)
    {
      hitlist_temp.unique(Hit::sameChannel);
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
      fNHits[fAllHits[i].detectorID]++;
    }

  fNHits[0] = fAllHits.size();
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
	  if(Hit::sameChannel(*iter, *jter))
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

void SRawEvent::clear()
{
  fAllHits.clear();
  for(Int_t i = 0; i < nChamberPlanes+1; i++)
    {
      fNHits[i] = 0;
    }

  fRunID = -1;
  fEventID = -1;
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
