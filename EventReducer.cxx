#include "EventReducer.h"

EventReducer::EventReducer(TString options) : afterhit(false), hodomask(false), outoftime(false), decluster(false), mergehodo(false), triggermask(false), sagitta(false), hough(false), externalpar(false), realization(false)
{
  //parse the reducer setup
  options.ToLower();
  if(options.Contains("a")) afterhit = true;
  if(options.Contains("h")) hodomask = true;
  if(options.Contains("o")) outoftime = true;
  if(options.Contains("c")) decluster = true;
  if(options.Contains("m")) mergehodo = true;
  if(options.Contains("t")) triggermask = true;
  if(options.Contains("s")) sagitta = true;
  if(options.Contains("g")) hough = true;
  if(options.Contains("e")) externalpar = true;
  if(options.Contains("r")) realization = true;

  //initialize services
  p_geomSvc = GeomSvc::instance();
  if(triggermask)
    {
      p_triggerAna = new TriggerAnalyzer();
      p_triggerAna->init();
      p_triggerAna->buildTriggerTree();
    }

  //set random seed
  rndm.SetSeed(0);
}

EventReducer::~EventReducer()
{
  if(triggermask)
    {
      delete p_triggerAna;
    }
}

int EventReducer::reduceEvent(SRawEvent* rawEvent)
{
  int nHits_before = rawEvent->getNHitsAll();

  //Label the hits which are not on an active trigger road as intime and trigger masked
  if(triggermask)
    {
      p_triggerAna->trimEvent(rawEvent);
    }

  //dump the vector of hits from SRawEvent to a list first
  hitlist.clear();
  for(std::vector<Hit>::iterator iter = rawEvent->fAllHits.begin(); iter != rawEvent->fAllHits.end(); ++iter)
    {
      if(realization && iter->detectorID <= 24)
	{
	  if(rndm.Rndm() > 0.94) continue;
	}

      if(outoftime && (!iter->isInTime())) continue;
      if(hodomask && iter->detectorID <= 24 && (!iter->isHodoMask())) continue;
      if(triggermask && iter->detectorID > 24 && iter->detectorID <= 40 && (!iter->isTriggerMask())) continue;

      //only temporary before the mapping is fixed
      if((iter->detectorID == 17 || iter->detectorID == 18) && iter->elementID >= 97 && iter->elementID <= 104)
	{
	  iter->detectorID = iter->detectorID == 17 ? 18 : 17;
	  iter->pos = p_geomSvc->getMeasurement(iter->detectorID, iter->elementID);
	  //iter->driftDistance = p_geomSvc->getDriftDistance(iter->detectorID, iter->tdcTime);
	  //iter->setInTime(p_geomSvc->isInTime(iter->detectorID, iter->tdcTime));
	}

      if(externalpar)
	{
	  iter->pos = p_geomSvc->getMeasurement(iter->detectorID, iter->elementID);
	  iter->driftDistance = p_geomSvc->getDriftDistance(iter->detectorID, iter->tdcTime);
	  //iter->setInTime(p_geomSvc->isInTime(iter->detectorID, iter->tdcTime));
	}

      if(realization && iter->detectorID <= 24)
	{
	  iter->driftDistance += rndm.Gaus(0., 0.04);
	}

      hitlist.push_back(*iter);
    }

  //Remove after hits
  hitlist.sort();
  if(afterhit) hitlist.unique();

  //Remove hit clusters
  if(decluster) deClusterize();

  //Remove the hits by sagitta ratio
  if(sagitta) sagittaReducer();

  //Push everything back to SRawEvent
  rawEvent->fAllHits.clear();
  rawEvent->fAllHits.assign(hitlist.begin(), hitlist.end());

  rawEvent->reIndex();
  return nHits_before - rawEvent->fNHits[0];
}

void EventReducer::sagittaReducer()
{
  //find index for D1, D2, and D3
  int nHits_D1 = 0;
  int nHits_D2 = 0;
  int nHits_D3 = 0;
  for(std::list<Hit>::iterator iter = hitlist.begin(); iter != hitlist.end(); ++iter)
    {
      if(iter->detectorID > 24) break;
      if(iter->detectorID <= 6)
	{
	  ++nHits_D1;
	}
      else if(iter->detectorID <= 12)
	{
	  ++nHits_D2;
	}
      else
	{
	  ++nHits_D3;
	}
    }
  int idx_D1 = nHits_D1;
  int idx_D2 = nHits_D1 + nHits_D2;
  int idx_D3 = nHits_D1 + nHits_D2 + nHits_D3;

  //Loop over all hits
  std::vector<Hit> hitTemp;
  hitTemp.assign(hitlist.begin(), hitlist.end());

  std::vector<int> flag(hitTemp.size(), -1);
  for(int i = idx_D2; i < idx_D3; ++i)
    {
      double z3 = p_geomSvc->getPlanePosition(hitTemp[i].detectorID);
      double slope = hitTemp[i].pos/z3;
      for(int j = idx_D1; j < idx_D2; ++j)
	{
	  if(p_geomSvc->getPlaneType(hitTemp[i].detectorID) != p_geomSvc->getPlaneType(hitTemp[j].detectorID)) continue;

	  double z2 = p_geomSvc->getPlanePosition(hitTemp[j].detectorID);
	  if(fabs((hitTemp[i].pos - hitTemp[j].pos)/(z2 - z3)) > TX_MAX) continue;
	  double s2 = hitTemp[j].pos - z2*slope;

	  for(int k = 0; k < idx_D1; ++k)
	    {
	      if(p_geomSvc->getPlaneType(hitTemp[i].detectorID) != p_geomSvc->getPlaneType(hitTemp[k].detectorID)) continue;
	      double s1 = hitTemp[k].pos - slope*p_geomSvc->getPlanePosition(hitTemp[k].detectorID);
	
	      if(fabs(s1/s2 - 1.77) < 0.25)
		{
		  flag[i] = 1;
		  flag[j] = 1;
		  flag[k] = 1;
		}
	    }
	}
    }

  int idx = 0;
  for(std::list<Hit>::iterator iter = hitlist.begin(); iter != hitlist.end(); )
    {
      if(flag[idx] < 0)
	{
	  iter = hitlist.erase(iter);
	}
      else
	{
	  ++iter;
	}

      ++idx;
      if(idx >= idx_D3) break;
    }
}

void EventReducer::deClusterize()
{
  std::vector<std::list<Hit>::iterator> cluster;
  cluster.clear();
  for(std::list<Hit>::iterator hit = hitlist.begin(); hit != hitlist.end(); ++hit)
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
	      processCluster(cluster);
	      cluster.push_back(hit);
	    }
	  else if(hit->elementID - cluster.back()->elementID > 1)
	    {
	      processCluster(cluster);
	      cluster.push_back(hit);
	    }
	  else
	    {
	      cluster.push_back(hit);
	    }
	}
    }
}

void EventReducer::processCluster(std::vector<std::list<Hit>::iterator>& cluster)
{
  unsigned int clusterSize = cluster.size();

  //size-2 clusters, retain the hit with smaller driftDistance
  if(clusterSize == 2)
    {
      double w_max = 0.9*0.5*(cluster.back()->pos - cluster.front()->pos);
      double w_min = w_max/9.*4.; //double w_min = 0.6*0.5*(cluster.back()->pos - cluster.front()->pos);

      if((cluster.front()->driftDistance > w_max && cluster.back()->driftDistance > w_min) || (cluster.front()->driftDistance > w_min && cluster.back()->driftDistance > w_max))
	{
	  cluster.front()->driftDistance > cluster.back()->driftDistance ? hitlist.erase(cluster.front()) : hitlist.erase(cluster.back());
	}
      else if(fabs(cluster.front()->tdcTime - cluster.back()->tdcTime) < 8. && cluster.front()->detectorID >= 13 && cluster.front()->detectorID <= 18)
	{
	  hitlist.erase(cluster.front());
	  hitlist.erase(cluster.back());
	}
    }

  //size-larger-than-3, discard entirely
  if(clusterSize >= 3)
    {
      double dt_mean = 0.;
      for(unsigned int i = 1; i < clusterSize; ++i)
	{
	  dt_mean += fabs(cluster[i]->tdcTime - cluster[i-1]->tdcTime);
	}
      dt_mean = dt_mean/(clusterSize - 1);

      if(dt_mean < 10.)
	{
	  //electric noise, discard them all
	  for(unsigned int i = 0; i < clusterSize; ++i)
	    {
	      hitlist.erase(cluster[i]);
	    }
	}
      else
	{
	  /*
	  double dt_rms = 0.;
       	  for(unsigned int i = 1; i < clusterSize; ++i)
 	    {
	       double dt = fabs(cluster[i]->tdcTime - cluster[i-1]->tdcTime);
	       dt_rms += ((dt - dt_mean)*(dt - dt_mean));
	    }
	  dt_rms = sqrt(dt_rms/(clusterSize - 1));
	
	  //delta ray, keep the first and last
	  if(dt_rms < 5.)*/
	    {
	      for(unsigned int i = 1; i < clusterSize - 1; ++i)
		{
		  hitlist.erase(cluster[i]);
		}
	    }
	}
    }

  cluster.clear();
}
