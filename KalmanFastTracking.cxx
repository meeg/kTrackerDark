/*
KalmanFastTracking.cxx

Implementation of class Tracklet, KalmanFastTracking

Author: Kun Liu, liuk@fnal.gov
Created: 05-28-2013
*/

#include <iostream>
#include <algorithm>
#include <cmath>

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TBox.h>
#include <TMatrixD.h>

#include "KalmanFitter.h"
#include "KalmanFastTracking.h"

KalmanFastTracking::KalmanFastTracking(bool flag)
{
  using namespace std;
#ifdef _DEBUG_ON
  cout << "Initialization of KalmanFastTracking ..." << endl;
  cout << "========================================" << endl;
#endif

  enable_KF = flag;

  //Initialize Kalman fitter
  if(enable_KF)
    {
      kmfitter = new KalmanFitter();
      kmfitter->setControlParameter(50, 0.001);
    }

  //Initialize minuit minimizer
  minimizer[0] = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex");
  minimizer[1] = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Combined");
  if(KMAG_ON == 1)
    {
      fcn = ROOT::Math::Functor(&tracklet_curr, &Tracklet::Eval, 5);
    }
  else
    {
      fcn = ROOT::Math::Functor(&tracklet_curr, &Tracklet::Eval, 4);
    }

  for(int i = 0; i < 2; ++i)
    {
      minimizer[i]->SetMaxFunctionCalls(1000000);
      minimizer[i]->SetMaxIterations(100);
      minimizer[i]->SetTolerance(1E-2);
      minimizer[i]->SetFunction(fcn);
      minimizer[i]->SetPrintLevel(0);
    }

  //Minimize ROOT output
  extern Int_t gErrorIgnoreLevel;
  gErrorIgnoreLevel = 9999;

  //Initialize geometry service
  p_geomSvc = GeomSvc::instance();

  //Initialize plane angles for all planes
  for(int i = 1; i <= 24; i++)
    {
      costheta_plane[i] = p_geomSvc->getCostheta(i);
      sintheta_plane[i] = p_geomSvc->getSintheta(i);
    }

  //Initialize hodoscope IDs
  detectorIDs_mask[0] = p_geomSvc->getDetectorIDs("H1"); 
  detectorIDs_mask[1] = p_geomSvc->getDetectorIDs("H2"); 
  detectorIDs_mask[2] = p_geomSvc->getDetectorIDs("H3"); 
  detectorIDs_mask[3] = p_geomSvc->getDetectorIDs("H4"); 
  detectorIDs_maskX[0] = p_geomSvc->getDetectorIDs("H1[TB]"); 
  detectorIDs_maskX[1] = p_geomSvc->getDetectorIDs("H2[TB]"); 
  detectorIDs_maskX[2] = p_geomSvc->getDetectorIDs("H3[TB]"); 
  detectorIDs_maskX[3] = p_geomSvc->getDetectorIDs("H4[TB]"); 
  detectorIDs_maskY[0] = p_geomSvc->getDetectorIDs("H1[LR]"); 
  detectorIDs_maskY[1] = p_geomSvc->getDetectorIDs("H2[LR]"); 
  detectorIDs_maskY[2] = p_geomSvc->getDetectorIDs("H4Y1[LR]"); 
  detectorIDs_maskY[3] = p_geomSvc->getDetectorIDs("H4Y2[LR]"); 

  //Masking stations for tracklets in station-1, 2, 3+/-
  stationIDs_mask[0].push_back(1);
  stationIDs_mask[1].push_back(2);
  stationIDs_mask[2].push_back(3);
  stationIDs_mask[3].push_back(3);

  //Masking stations for back partial
  stationIDs_mask[4].push_back(2);
  stationIDs_mask[4].push_back(3);
  stationIDs_mask[4].push_back(4);

  //Masking stations for global track
  stationIDs_mask[5].push_back(1);
  stationIDs_mask[5].push_back(2);
  stationIDs_mask[5].push_back(3);
  stationIDs_mask[5].push_back(4);

  //prop. tube IDs for mu id
  detectorIDs_muid[0][0] = 43;
  detectorIDs_muid[0][1] = 44;
  detectorIDs_muid[0][2] = 45;
  detectorIDs_muid[0][3] = 46;
  detectorIDs_muid[1][0] = 41;
  detectorIDs_muid[1][1] = 42;
  detectorIDs_muid[1][2] = 47;
  detectorIDs_muid[1][3] = 48;

  //Initialize masking window sizes, with 15% contingency, for station-2, increase that to 20%
  for(int i = 25; i <= 48; i++)
    {
      double factor = 0.15;
      if(i > 28 && i < 33) factor = 0.2; //for station-2

      z_mask[i-25] = p_geomSvc->getPlanePosition(i);
      for(int j = 1; j <= p_geomSvc->getPlaneNElements(i); j++)
	{
	  double x_min, x_max, y_min, y_max;
	  p_geomSvc->get2DBoxSize(i, j, x_min, x_max, y_min, y_max);
	  
	  x_min -= (factor*(x_max - x_min));
	  x_max += (factor*(x_max - x_min));
	  y_min -= (factor*(y_max - y_min));
	  y_max += (factor*(y_max - y_min));

	  x_mask_min[i-25][j-1] = x_min;
	  x_mask_max[i-25][j-1] = x_max;
	  y_mask_min[i-25][j-1] = y_min;
	  y_mask_max[i-25][j-1] = y_max;
	}      
    }

#ifdef _DEBUG_ON
  cout << "========================" << endl;
  cout << "Hodo. masking settings: " << endl;
  for(int i = 0; i < 5; i++)
    {
      cout << "For station " << i+1 << endl;
      for(std::vector<int>::iterator iter = detectorIDs_mask[i].begin(); iter != detectorIDs_mask[i].end(); ++iter) cout << "All: " << *iter << endl;
      for(std::vector<int>::iterator iter = detectorIDs_maskX[i].begin(); iter != detectorIDs_maskX[i].end(); ++iter) cout << "X: " << *iter << endl;
      for(std::vector<int>::iterator iter = detectorIDs_maskY[i].begin(); iter != detectorIDs_maskY[i].end(); ++iter) cout << "Y: " << *iter << endl;
    }
  
  for(int i = 0; i < 6; ++i)
    {
      std::cout << "Masking stations for tracklets with stationID = " << i + 1 << ": " << std::endl; 
      for(std::vector<int>::iterator iter = stationIDs_mask[i].begin(); iter != stationIDs_mask[i].end(); ++iter)
	{
	  std::cout << *iter << "  ";
	}
      std::cout << std::endl;
    }
 
#endif

  //Initialize super stationIDs
  for(int i = 0; i < 4; i++) superIDs[i].clear();
  superIDs[0].push_back((p_geomSvc->getDetectorIDs("D1X")[0] + 1)/2);
  superIDs[0].push_back((p_geomSvc->getDetectorIDs("D1U")[0] + 1)/2);
  superIDs[0].push_back((p_geomSvc->getDetectorIDs("D1V")[0] + 1)/2);
  superIDs[1].push_back((p_geomSvc->getDetectorIDs("D2X")[0] + 1)/2);
  superIDs[1].push_back((p_geomSvc->getDetectorIDs("D2U")[0] + 1)/2);
  superIDs[1].push_back((p_geomSvc->getDetectorIDs("D2V")[0] + 1)/2);
  superIDs[2].push_back((p_geomSvc->getDetectorIDs("D3pX")[0] + 1)/2);
  superIDs[2].push_back((p_geomSvc->getDetectorIDs("D3pU")[0] + 1)/2);
  superIDs[2].push_back((p_geomSvc->getDetectorIDs("D3pV")[0] + 1)/2);
  superIDs[3].push_back((p_geomSvc->getDetectorIDs("D3mX")[0] + 1)/2);
  superIDs[3].push_back((p_geomSvc->getDetectorIDs("D3mU")[0] + 1)/2);
  superIDs[3].push_back((p_geomSvc->getDetectorIDs("D3mV")[0] + 1)/2);

#ifdef _DEBUG_ON
  cout << "=============" << endl;
  cout << "Chamber IDs: " << endl;
  for(int i = 0; i < 4; i++)
    {
      for(int j = 0; j < 3; j++) cout << i << "  " << j << ": " << superIDs[i][j] << endl;
    }

  //Initialize widow sizes for X-U matching and z positions of all chambers
  cout << "======================" << endl;
  cout << "U plane window sizes: " << endl;
#endif
  double u_factor[] = {5., 5., 15., 15.};
  for(int i = 0; i < 4; i++)
    {
      int xID = 2*superIDs[i][0] - 1;
      int uID = 2*superIDs[i][1] - 1;
      int vID = 2*superIDs[i][2] - 1;
      double spacing = p_geomSvc->getPlaneSpacing(uID);
      double x_span = p_geomSvc->getPlaneScaleY(uID);

      z_plane_x[i] = 0.5*(p_geomSvc->getPlanePosition(xID) + p_geomSvc->getPlanePosition(xID+1));
      z_plane_u[i] = 0.5*(p_geomSvc->getPlanePosition(uID) + p_geomSvc->getPlanePosition(uID+1));
      z_plane_v[i] = 0.5*(p_geomSvc->getPlanePosition(vID) + p_geomSvc->getPlanePosition(vID+1));

      u_costheta[i] = costheta_plane[uID];
      u_sintheta[i] = sintheta_plane[uID];
      
      //u_win[i] = fabs(0.5*x_span/(spacing/sintheta_plane[uID])) + 2.*spacing + u_factor[i];
      u_win[i] = fabs(0.5*x_span*sintheta_plane[uID]) + TX_MAX*fabs((z_plane_u[i] - z_plane_x[i])*u_costheta[i]) + TY_MAX*fabs((z_plane_u[i] - z_plane_x[i])*u_sintheta[i]) + 2.*spacing + u_factor[i];

#ifdef _DEBUG_ON
      cout << "Station " << i << ": " << u_win[i] << endl; 
#endif
    }

  //Initialize Z positions and maximum parameters of all planes
  for(int i = 1; i <= 24; i++) 
    {
      z_plane[i] = p_geomSvc->getPlanePosition(i);
      slope_max[i] = costheta_plane[i]*TX_MAX + sintheta_plane[i]*TY_MAX;
      intersection_max[i] = costheta_plane[i]*X0_MAX + sintheta_plane[i]*Y0_MAX;

#ifdef COARSE_MODE
      resol_plane[i] = p_geomSvc->getPlaneSpacing(i)/sqrt(12.);
#else 
      resol_plane[i] = p_geomSvc->getPlaneResolution(i);
#endif

      spacing_plane[i] = p_geomSvc->getPlaneSpacing(i);
    }

#ifdef _DEBUG_ON
  cout << "======================================" << endl;
  cout << "Maximum local slope and intersection: " << endl;
#endif
  for(int i = 1; i <= 12; i++)
    {
      double d_slope = (p_geomSvc->getPlaneResolution(2*i - 1) + p_geomSvc->getPlaneResolution(2*i))/(z_plane[2*i] - z_plane[2*i-1]);
      double d_intersection = d_slope*z_plane[2*i];

      slope_max[2*i-1] += d_slope;
      intersection_max[2*i-1] += d_intersection;
      slope_max[2*i] += d_slope;
      intersection_max[2*i] += d_intersection;

#ifdef _DEBUG_ON
      cout << "Super plane " << i << ": " << slope_max[2*i-1] << "  " << intersection_max[2*i-1] << endl; 
#endif
    }

  //Initialize sagitta ratios, index 0, 1, 2 are for U, X, V
  s_ratio[0] = 1.77; s_sigma[0] = 0.2; s_detectorID[0] = 12;
  s_ratio[1] = 1.77; s_sigma[1] = 0.2; s_detectorID[1] = 10;
  s_ratio[2] = 1.77; s_sigma[2] = 0.2; s_detectorID[2] = 8;
}

KalmanFastTracking::~KalmanFastTracking()
{
  if(enable_KF) delete kmfitter;
  delete minimizer[0];
  delete minimizer[1];
}

bool KalmanFastTracking::setRawEvent(SRawEvent* event_input)
{
  rawEvent = event_input;
  if(!acceptEvent(rawEvent)) return false;
  hitAll = event_input->getAllHits();
#ifdef _DEBUG_ON
  for(std::vector<Hit>::iterator iter = hitAll.begin(); iter != hitAll.end(); ++iter) iter->print();
#endif

  //Initialize hodo and masking IDs
  for(int i = 0; i < 4; i++)
    {
      //std::cout << "For station " << i << std::endl;
      hitIDs_mask[i].clear();
      hitIDs_mask[i] = rawEvent->getHitsIndexInDetectors(detectorIDs_maskX[i]);

      //for(std::list<int>::iterator iter = hitIDs_mask[i].begin(); iter != hitIDs_mask[i].end(); ++iter) std::cout << *iter << " " << hitAll[*iter].detectorID << " === ";
      //std::cout << std::endl;
    }

  //Initialize prop. tube IDs
  for(int i = 0; i < 2; ++i)
    {
      for(int j = 0; j < 4; ++j)
	{
	  hitIDs_muid[i][j].clear();
	  hitIDs_muid[i][j] = rawEvent->getHitsIndexInDetector(detectorIDs_muid[i][j]);
	}
    }

  //Initialize tracklet lists
  for(int i = 0; i < 5; i++) trackletsInSt[i].clear();
  tracks.clear();

  //Build tracklets in station 2, 3+, 3-
  //When i = 3, works for st3+, for i = 4, works for st3-
  buildTrackletsInStation(2); 
  if(trackletsInSt[1].empty()) 
    {
#ifdef _DEBUG_ON
      LogInfo("Failed in tracklet build at station 2");
#endif
      return false;
    }

  buildTrackletsInStation(3); buildTrackletsInStation(4);
  if(trackletsInSt[2].empty()) 
    {
#ifdef _DEBUG_ON
      LogInfo("Failed in tracklet build at station 3");
#endif
      return false;
    }

  //Build back partial tracks in station 2, 3+ and 3-
  buildBackPartialTracks();
 
  //Connect tracklets in station 2/3 and station 1 to form global tracks
  buildGlobalTracks();
 
#ifdef _DEBUG_ON 
  for(int i = 0; i <= 4; i++)
    {
      std::cout << "=======================================================================================" << std::endl;
      LogInfo("Final tracklets in station: " << i+1 << " is " << trackletsInSt[i].size()); 
      for(std::list<Tracklet>::iterator tracklet = trackletsInSt[i].begin(); tracklet != trackletsInSt[i].end(); ++tracklet)
	{
	  tracklet->print();
	}
      std::cout << "=======================================================================================" << std::endl;
   }
#endif

  if(trackletsInSt[4].empty()) return false;
  if(!enable_KF) return true;

  //If there is no possibility of a dimuon, return 
  if(DIMUON_MODE == 1)
    {
      int nPlus = 0;
      int nMinus = 0;
      for(std::list<Tracklet>::iterator tracklet = trackletsInSt[4].begin(); tracklet != trackletsInSt[4].end(); ++tracklet)
	{
	  if(tracklet->getCharge() > 0) 
	    {
	      ++nPlus;
	    }
	  else
	    {
	      ++nMinus;
	    }
	}

      if(nPlus < 1 || nMinus < 1) return false;
    }
 
  //Build kalman tracks
  for(std::list<Tracklet>::iterator tracklet = trackletsInSt[4].begin(); tracklet != trackletsInSt[4].end(); ++tracklet)
    {
      processOneTracklet(*tracklet);
    }

#ifdef _DEBUG_ON
  LogInfo(tracks.size() << " final tracks:");
  for(std::list<KalmanTrack>::iterator kmtrk = tracks.begin(); kmtrk != tracks.end(); ++kmtrk)
    {
      kmtrk->print();
    }
#endif

  return true;
}

bool KalmanFastTracking::acceptEvent(SRawEvent* rawEvent)
{
#ifdef _DEBUG_ON
  LogInfo("D1: " << rawEvent->getNHitsInD1());
  LogInfo("D2: " << rawEvent->getNHitsInD2());
  LogInfo("D3p: " << rawEvent->getNHitsInD3p());
  LogInfo("D3m: " << rawEvent->getNHitsInD3m());
  LogInfo("H1: " << rawEvent->getNHitsInDetectors(detectorIDs_maskX[0]));
  LogInfo("H2: " << rawEvent->getNHitsInDetectors(detectorIDs_maskX[1]));
  LogInfo("H3: " << rawEvent->getNHitsInDetectors(detectorIDs_maskX[2]));
  LogInfo("H4: " << rawEvent->getNHitsInDetectors(detectorIDs_maskX[3]));
#endif

  if(rawEvent->getNHitsInD1() > 200) return false;
  if(rawEvent->getNHitsInD2() > 100) return false;
  if(rawEvent->getNHitsInD3p() > 100) return false;
  if(rawEvent->getNHitsInD3m() > 100) return false;
  if(rawEvent->getNHitsInDetectors(detectorIDs_maskX[0]) > 15) return false;
  if(rawEvent->getNHitsInDetectors(detectorIDs_maskX[1]) > 10) return false;
  if(rawEvent->getNHitsInDetectors(detectorIDs_maskX[2]) > 10) return false;
  if(rawEvent->getNHitsInDetectors(detectorIDs_maskX[3]) > 10) return false;
  if(rawEvent->getNPropHitsAll() > 250) return false;  
  return true;
}

void KalmanFastTracking::buildBackPartialTracks()
{
#ifndef ALIGNMENT_MODE
  //Temporary container for a simple chisq fit
  int nHitsX2, nHitsX3;
  double z_fit[4], x_fit[4];
  double a, b;
#endif

  for(std::list<Tracklet>::iterator tracklet3 = trackletsInSt[2].begin(); tracklet3 != trackletsInSt[2].end(); ++tracklet3)
    {
#ifndef ALIGNMENT_MODE
      //Extract the X hits only from station-3 tracks
      nHitsX3 = 0;
      for(std::list<SignedHit>::iterator ptr_hit = tracklet3->hits.begin(); ptr_hit != tracklet3->hits.end(); ++ptr_hit)
	{
	  if(ptr_hit->hit.index < 0) continue;
	  if(p_geomSvc->getPlaneType(ptr_hit->hit.detectorID) == 1) 
	    {
	      z_fit[nHitsX3] = z_plane[ptr_hit->hit.detectorID];
	      x_fit[nHitsX3] = ptr_hit->hit.pos;
	      ++nHitsX3;
	    }
	}
#endif
      Tracklet tracklet_best;
      for(std::list<Tracklet>::iterator tracklet2 = trackletsInSt[1].begin(); tracklet2 != trackletsInSt[1].end(); ++tracklet2)
	{
#ifndef ALIGNMENT_MODE
	  //Extract the X hits from station-2 tracke
	  nHitsX2 = nHitsX3;
	  for(std::list<SignedHit>::iterator ptr_hit = tracklet2->hits.begin(); ptr_hit != tracklet2->hits.end(); ++ptr_hit)
	    {
	      if(ptr_hit->hit.index < 0) continue;
	      if(p_geomSvc->getPlaneType(ptr_hit->hit.detectorID) == 1) 
		{
		  z_fit[nHitsX2] = z_plane[ptr_hit->hit.detectorID];
		  x_fit[nHitsX2] = ptr_hit->hit.pos;
		  ++nHitsX2;
		}
	    }
	  
	  //Apply a simple linear fit to get rough estimation of X-Z slope and intersection
	  chi2fit(nHitsX2, z_fit, x_fit, a, b);
	  if(fabs(a) > TX_MAX || fabs(b) > X0_MAX) continue;

	  //Project to proportional tubes to see if there is enough
	  int nPropHits = 0;
	  for(int i = 0; i < 4; ++i)
	    {
	      double x_exp = a*z_mask[detectorIDs_muid[0][i] - 25] + b;
	      for(std::list<int>::iterator iter = hitIDs_muid[0][i].begin(); iter != hitIDs_muid[0][i].end(); ++iter)
		{
		  if(fabs(hitAll[*iter].pos - x_exp) < 2.54)
		    {
		      ++nPropHits;
		      break;
		    }
		}
	      if(nPropHits > 0) break;
	    }
	  if(nPropHits == 0) continue;
#endif

	  Tracklet tracklet_23 = (*tracklet2) + (*tracklet3);
#ifdef _DEBUG_ON
	  LogInfo("Using following two tracklets:");
	  tracklet2->print();
	  tracklet3->print();
	  LogInfo("Yield this combination:");
	  tracklet_23.print();
#endif
	  fitTracklet(tracklet_23);
	  if(tracklet_23.chisq > 3000.)
	    {
#ifdef _DEBUG_ON
	      tracklet_23.print();
	      LogInfo("Impossible combination!");
#endif
	      continue;
	    }

#ifndef COARSE_MODE
	  resolveLeftRight(tracklet_23, 25.);
	  resolveLeftRight(tracklet_23, 100.);
#endif
	  ///Remove bad hits if needed
	  //removeBadHits(tracklet_23);


#ifdef _DEBUG_ON
	  LogInfo("New tracklet: ");
	  tracklet_23.print(); 

	  LogInfo("Current best:");
	  tracklet_best.print();

	  LogInfo("Comparison: " << (tracklet_23 < tracklet_best));
	  LogInfo("Quality: " << acceptTracklet(tracklet_23));
#endif

	  //If current tracklet is better than the best tracklet up-to-now
	  if(acceptTracklet(tracklet_23) && tracklet_23 < tracklet_best)
	    {
	      tracklet_best = tracklet_23;
	    }
#ifdef _DEBUG_ON
	  else
	    {
	      LogInfo("Rejected!!");
	    }
#endif
	}

      if(tracklet_best.isValid()) trackletsInSt[3].push_back(tracklet_best);
    }

  reduceTrackletList(trackletsInSt[3]);
  trackletsInSt[3].sort();
}

void KalmanFastTracking::buildGlobalTracks()
{
  double pos_exp[3], window[3];
  for(std::list<Tracklet>::iterator tracklet23 = trackletsInSt[3].begin(); tracklet23 != trackletsInSt[3].end(); ++tracklet23)
    {
      //Calculate the window in station 1
      if(KMAG_ON)
	{
	  getSagittaWindowsInSt1(*tracklet23, pos_exp, window);
	}
      else
	{
	  getExtrapoWindowsInSt1(*tracklet23, pos_exp, window);
	}

#ifdef _DEBUG_ON
      LogInfo("Using this back partial: ");
      tracklet23->print();
      for(int i = 0; i < 3; i++) LogInfo("Extrapo: " << pos_exp[i] << "  " << window[i]);
#endif

      trackletsInSt[0].clear();
      buildTrackletsInStation(1, pos_exp, window);
      
      Tracklet tracklet_best;
      for(std::list<Tracklet>::iterator tracklet1 = trackletsInSt[0].begin(); tracklet1 != trackletsInSt[0].end(); ++tracklet1)
	{
#ifdef _DEBUG_ON
	  LogInfo("With this station 1 track:");
	  tracklet1->print();
#endif

	  Tracklet tracklet_global = (*tracklet23) * (*tracklet1);
	  fitTracklet(tracklet_global);
  
#ifndef COARSE_MODE
	  ///Resolve the left-right with a tight pull cut, then a loose one, then resolve by single projections
	  resolveLeftRight(tracklet_global, 100.);
	  resolveLeftRight(tracklet_global, 1000.);
          resolveSingleLeftRight(tracklet_global);
#endif
	  ///Remove bad hits if needed
	  removeBadHits(tracklet_global);

#ifdef _DEBUG_ON
	  LogInfo("New tracklet: ");
	  tracklet_global.print(); 

	  LogInfo("Current best:");
	  tracklet_best.print();

	  LogInfo("Comparison: " << (tracklet_global < tracklet_best));
	  LogInfo("Quality   : " << acceptTracklet(tracklet_global));
#endif
	  if(acceptTracklet(tracklet_global) && tracklet_global < tracklet_best)
	    {
#ifdef _DEBUG_ON
	      LogInfo("Accepted!!!");
#endif
	      tracklet_best = tracklet_global;
    	    }
#ifdef _DEBUG_ON
	    {
	      LogInfo("Rejected!!!");
	    }
#endif
	}

      if(tracklet_best.isValid()) trackletsInSt[4].push_back(tracklet_best);
    }

  trackletsInSt[4].sort();
}

void KalmanFastTracking::resolveLeftRight(Tracklet& tracklet, double threshold)
{
#ifdef _DEBUG_ON
  LogInfo("Left right for this track..");
  tracklet.print();
#endif

  //Check if the track has been updated
  bool isUpdated = false;

  //Four possibilities
  int possibility[4][2] = {{1, 1}, {1, -1}, {-1, 1}, {-1, -1}};

  //Total number of hit pairs in this tracklet
  int nPairs = tracklet.hits.size()/2;
  
  int nResolved = 0;
  std::list<SignedHit>::iterator hit1 = tracklet.hits.begin();
  std::list<SignedHit>::iterator hit2 = tracklet.hits.begin(); ++hit2;
  while(true)
    {
#ifdef _DEBUG_ON
      LogInfo(hit1->hit.index << "  " << hit2->sign << " === " << hit2->hit.index << "  " << hit2->sign);
      int detectorID1 = hit1->hit.detectorID;
      int detectorID2 = hit2->hit.detectorID;
      LogInfo("Hit1: " << tracklet.getExpPositionX(z_plane[detectorID1])*costheta_plane[detectorID1] + tracklet.getExpPositionY(z_plane[detectorID1])*sintheta_plane[detectorID1] << "  " << hit1->hit.pos + hit1->hit.driftDistance << "  " << hit1->hit.pos - hit1->hit.driftDistance);
      LogInfo("Hit2: " << tracklet.getExpPositionX(z_plane[detectorID2])*costheta_plane[detectorID2] + tracklet.getExpPositionY(z_plane[detectorID2])*sintheta_plane[detectorID2] << "  " << hit2->hit.pos + hit2->hit.driftDistance << "  " << hit2->hit.pos - hit2->hit.driftDistance);
#endif

      if(hit1->hit.index > 0 && hit2->hit.index > 0 && hit1->sign*hit2->sign == 0)
	{
	  int index_min = -1;
	  double pull_min = 1E6;
	  for(int i = 0; i < 4; i++)
	    {
	      double slope_local = (hit1->pos(possibility[i][0]) - hit2->pos(possibility[i][1]))/(z_plane[hit1->hit.detectorID] - z_plane[hit2->hit.detectorID]);
	      double inter_local = hit1->pos(possibility[i][0]) - slope_local*z_plane[hit1->hit.detectorID];
	   
	      if(fabs(slope_local) > slope_max[hit1->hit.detectorID] || fabs(inter_local) > intersection_max[hit1->hit.detectorID]) continue;

	      double tx, ty, x0, y0;
	      double err_tx, err_ty, err_x0, err_y0;
	      if(tracklet.stationID == 6 && hit1->hit.detectorID <= 6)
		{
		  tracklet.getXZInfoInSt1(tx, x0);
		  tracklet.getXZErrorInSt1(err_tx, err_x0);
		}
	      else
		{
		  tx = tracklet.tx;
		  x0 = tracklet.x0;
		  err_tx = tracklet.err_tx;
		  err_x0 = tracklet.err_x0;
		}
	      ty = tracklet.ty;
	      y0 = tracklet.y0;
	      err_ty = tracklet.err_ty;
	      err_y0 = tracklet.err_y0;

	      double slope_exp = costheta_plane[hit1->hit.detectorID]*tx + sintheta_plane[hit1->hit.detectorID]*ty;
	      double err_slope = fabs(costheta_plane[hit1->hit.detectorID]*err_tx) + fabs(sintheta_plane[hit2->hit.detectorID]*err_ty);
	      double inter_exp = costheta_plane[hit1->hit.detectorID]*x0 + sintheta_plane[hit1->hit.detectorID]*y0;
	      double err_inter = fabs(costheta_plane[hit1->hit.detectorID]*err_x0) + fabs(sintheta_plane[hit2->hit.detectorID]*err_y0);

	      double pull = sqrt((slope_exp - slope_local)*(slope_exp - slope_local)/err_slope/err_slope + (inter_exp - inter_local)*(inter_exp - inter_local)/err_inter/err_inter);
	      if(pull < pull_min)
		{
		  index_min = i;
		  pull_min = pull;
		}

#ifdef _DEBUG_ON
	      LogInfo(hit1->hit.detectorID << ": " << i << "  " << possibility[i][0] << "  " << possibility[i][1]);
	      LogInfo(tx << "  " << x0 << "  " << ty << "  " << y0);
	      LogInfo("Slope: " << slope_local << "  " << slope_exp << "  " << err_slope);
	      LogInfo("Intersection: " << inter_local << "  " << inter_exp << "  " << err_inter);
	      LogInfo("Current: " << pull << "  " << index_min << "  " << pull_min);
#endif
	    }

	  //LogInfo("Final: " << index_min << "  " << pull_min);
	  if(index_min >= 0 && pull_min < threshold)//((tracklet.stationID == 5 && pull_min < 25.) || (tracklet.stationID == 6 && pull_min < 100.)))
	    {
	      hit1->sign = possibility[index_min][0];
	      hit2->sign = possibility[index_min][1];

	      isUpdated = true;
	    }
	}

      ++nResolved;
      if(nResolved >= nPairs) break;

      ++hit1; ++hit1;
      ++hit2; ++hit2;
    }

  if(isUpdated) fitTracklet(tracklet);
}

void KalmanFastTracking::resolveSingleLeftRight(Tracklet& tracklet)
{
#ifdef _DEBUG_ON
  LogInfo("Single left right for this track..");
  tracklet.print();
#endif

  //Check if the track has been updated
  bool isUpdated = false;

  for(std::list<SignedHit>::iterator hit_sign = tracklet.hits.begin(); hit_sign != tracklet.hits.end(); ++hit_sign)
    {
      if(hit_sign->hit.index < 0 || hit_sign->sign != 0) continue;

      int detectorID = hit_sign->hit.detectorID;
      double pos_exp = tracklet.getExpPositionX(z_plane[detectorID])*costheta_plane[detectorID] + tracklet.getExpPositionY(z_plane[detectorID])*sintheta_plane[detectorID];
      hit_sign->sign = pos_exp > hit_sign->hit.pos ? 1 : -1;

      isUpdated = true;
    }

  if(isUpdated) fitTracklet(tracklet);
}

void KalmanFastTracking::removeBadHits(Tracklet& tracklet)
{
#ifdef _DEBUG_ON
  LogInfo("Removing hits for this track..");
  tracklet.calcChisq();
  tracklet.print();
#endif

  //Check if the track has beed updated
  bool isUpdated = true;
  while(isUpdated)
    {
      isUpdated = false;
      tracklet.calcChisq();

      SignedHit* hit_remove = NULL;
      double res_remove = -1.;
      for(std::list<SignedHit>::iterator hit_sign = tracklet.hits.begin(); hit_sign != tracklet.hits.end(); ++hit_sign)
	{
	  if(hit_sign->hit.index < 0) continue;

	  int detectorID = hit_sign->hit.detectorID; 
	  double res_curr = fabs(tracklet.residual[detectorID-1]);
	  if(res_remove < res_curr)
	    {
	      res_remove = res_curr;
	      hit_remove = &(*hit_sign);
	    }
	}

      if(hit_remove != NULL && res_remove > HIT_REJECT*resol_plane[hit_remove->hit.detectorID])
	{
#ifdef _DEBUG_ON
	  LogInfo("Dropping this hit: " << res_remove << "  " << HIT_REJECT*resol_plane[hit_remove->hit.detectorID]);
	  hit_remove->hit.print();
#endif

	  hit_remove->hit.index = -1;
	  int planeType = p_geomSvc->getPlaneType(hit_remove->hit.detectorID);
	  if(planeType == 1)
	    {
	      --tracklet.nXHits;
	    } 
	  else if(planeType == 2)
	    {
	      --tracklet.nUHits;
	    }
	  else
	    {
	      --tracklet.nVHits;
	    }

	  isUpdated = true;
	}

      if(isUpdated) fitTracklet(tracklet);
    }
}

void KalmanFastTracking::resolveLeftRight(SRawEvent::hit_pair hpair, int& LR1, int& LR2)
{
  LR1 = 0;
  LR2 = 0;

  //If either hit is missing, no left-right can be assigned
  if(hpair.first < 0 || hpair.second < 0)
    {
      return;
    }

  int possibility[4][2] = {{1, 1}, {1, -1}, {-1, 1}, {-1, -1}};
  int nResolved = 0;
  for(int i = 0; i < 4; i++)
    {
      if(nResolved > 1) break;

      int hitID1 = hpair.first;
      int hitID2 = hpair.second;
      double slope_local = (hitAll[hitID1].pos + possibility[i][0]*hitAll[hitID1].driftDistance - hitAll[hitID2].pos - possibility[i][1]*hitAll[hitID2].driftDistance)/(z_plane[hitAll[hitID1].detectorID] - z_plane[hitAll[hitID2].detectorID]);
      double intersection_local = hitAll[hitID1].pos + possibility[i][0]*hitAll[hitID1].driftDistance - slope_local*z_plane[hitAll[hitID1].detectorID];
      
      //LogInfo(i << "  " << nResolved << "  " << slope_local << "  " << intersection_local);
      if(fabs(slope_local) < slope_max[hitAll[hitID1].detectorID] && fabs(intersection_local) < intersection_max[hitAll[hitID1].detectorID])
	{
  	  nResolved++;
	  LR1 = possibility[i][0];
	  LR2 = possibility[i][1];
	}
    }

  if(nResolved > 1)
    {
      LR1 = 0;
      LR2 = 0;
    }

  //LogInfo("Final: " << LR1 << "  " << LR2);
}

void KalmanFastTracking::buildTrackletsInStation(int stationID, double* pos_exp, double* window)
{
#ifdef _DEBUG_ON
  LogInfo("Building tracklets in station " << stationID);
#endif

  //actuall ID of the tracklet lists
  int sID = stationID - 1;
  int listID = sID;
  if(listID == 3) listID = 2;

  //Extract the X, U, V hit pairs
  std::list<SRawEvent::hit_pair> pairs_X, pairs_U, pairs_V;
  if(pos_exp == NULL)
    {
      pairs_X = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][0]);
      pairs_U = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][1]);
      pairs_V = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][2]);
    }
  else
    {
      //Note that in pos_exp[], index 0 stands for U, index 1 stands for X, index 2 stands for V
      pairs_X = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][0], pos_exp[1], window[1]);
      pairs_U = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][1], pos_exp[0], window[0]);
      pairs_V = rawEvent->getPartialHitPairsInSuperDetector(superIDs[sID][2], pos_exp[2], window[2]);
    }

#ifdef _DEBUG_ON
  LogInfo("Hit pairs in this event: ");
  for(std::list<SRawEvent::hit_pair>::iterator iter = pairs_X.begin(); iter != pairs_X.end(); ++iter) LogInfo("X :" << iter->first << "  " << iter->second << "  " << hitAll[iter->first].index << " " << (iter->second < 0 ? -1 : hitAll[iter->second].index));
  for(std::list<SRawEvent::hit_pair>::iterator iter = pairs_U.begin(); iter != pairs_U.end(); ++iter) LogInfo("U :" << iter->first << "  " << iter->second << "  " << hitAll[iter->first].index << " " << (iter->second < 0 ? -1 : hitAll[iter->second].index));
  for(std::list<SRawEvent::hit_pair>::iterator iter = pairs_V.begin(); iter != pairs_V.end(); ++iter) LogInfo("V :" << iter->first << "  " << iter->second << "  " << hitAll[iter->first].index << " " << (iter->second < 0 ? -1 : hitAll[iter->second].index));
#endif

  if(pairs_X.empty() || pairs_U.empty() || pairs_V.empty())
    {
#ifdef _DEBUG_ON
      LogInfo("Not all view has hits in station " << stationID);
#endif
      return;
    }

  //X-U combination first, then add V pairs
  for(std::list<SRawEvent::hit_pair>::iterator xiter = pairs_X.begin(); xiter != pairs_X.end(); ++xiter)
    {
      //U projections from X plane
      double x_pos = xiter->second >= 0 ? 0.5*(hitAll[xiter->first].pos + hitAll[xiter->second].pos) : hitAll[xiter->first].pos;
      double u_min = x_pos*u_costheta[sID] - u_win[sID];
      double u_max = u_min + 2.*u_win[sID];

#ifdef _DEBUG_ON
      LogInfo("Trying X hits " << xiter->first << "  " << xiter->second << "  " << hitAll[xiter->first].elementID << " at " << x_pos);
      LogInfo("U plane window:" << u_min << "  " << u_max);
#endif
      for(std::list<SRawEvent::hit_pair>::iterator uiter = pairs_U.begin(); uiter != pairs_U.end(); ++uiter)
	{
	  double u_pos = uiter->second >= 0 ? 0.5*(hitAll[uiter->first].pos + hitAll[uiter->second].pos) : hitAll[uiter->first].pos;
#ifdef _DEBUG_ON
	  LogInfo("Trying U hits " << uiter->first << "  " << uiter->second << "  " << hitAll[uiter->first].elementID << " at " << u_pos);
#endif
	  if(u_pos < u_min || u_pos > u_max) continue;

	  //V projections from X and U plane
          double z_x = xiter->second >= 0 ? z_plane_x[sID] : z_plane[hitAll[xiter->first].detectorID];
          double z_u = uiter->second >= 0 ? z_plane_u[sID] : z_plane[hitAll[uiter->first].detectorID];
          double z_v = z_plane_v[sID];
	  double v_win1 = spacing_plane[hitAll[uiter->first].detectorID]*2.*u_costheta[sID];
	  double v_win2 = fabs((z_u + z_v - 2.*z_x)*u_costheta[sID]*TX_MAX);
	  double v_win3 = fabs((z_v - z_u)*u_sintheta[sID]*TY_MAX);
	  double v_win = v_win1 + v_win2 + v_win3 + 2.*spacing_plane[hitAll[uiter->first].detectorID];
	  double v_min = 2*x_pos*u_costheta[sID] - u_pos - v_win;
	  double v_max = v_min + 2.*v_win;

#ifdef _DEBUG_ON	  
      	  LogInfo("V plane window:" << v_min << "  " << v_max);
#endif
	  for(std::list<SRawEvent::hit_pair>::iterator viter = pairs_V.begin(); viter != pairs_V.end(); ++viter)
	    {
	      double v_pos = viter->second >= 0 ? 0.5*(hitAll[viter->first].pos + hitAll[viter->second].pos) : hitAll[viter->first].pos;
#ifdef _DEBUG_ON
    	      LogInfo("Trying V hits " << viter->first << "  " << viter->second << "  " << hitAll[viter->first].elementID << " at " << v_pos);
#endif
	      if(v_pos < v_min || v_pos > v_max) continue;

	      //Now add the tracklet
	      int LR1 = 0;
	      int LR2 = 0;
	      Tracklet tracklet_new;
              tracklet_new.stationID = stationID;

	      //resolveLeftRight(*xiter, LR1, LR2);
	      if(xiter->first >= 0) { tracklet_new.hits.push_back(SignedHit(hitAll[xiter->first], LR1)); tracklet_new.nXHits++; }
	      if(xiter->second >= 0) { tracklet_new.hits.push_back(SignedHit(hitAll[xiter->second], LR2)); tracklet_new.nXHits++; }

	      //resolveLeftRight(*uiter, LR1, LR2);
	      if(uiter->first >= 0) { tracklet_new.hits.push_back(SignedHit(hitAll[uiter->first], LR1)); tracklet_new.nUHits++; }
	      if(uiter->second >= 0) { tracklet_new.hits.push_back(SignedHit(hitAll[uiter->second], LR2)); tracklet_new.nUHits++; }

	      //resolveLeftRight(*viter, LR1, LR2);
	      if(viter->first >= 0) { tracklet_new.hits.push_back(SignedHit(hitAll[viter->first], LR1)); tracklet_new.nVHits++; }
	      if(viter->second >= 0) { tracklet_new.hits.push_back(SignedHit(hitAll[viter->second], LR2)); tracklet_new.nVHits++; }
	   
	      tracklet_new.sortHits();
	      if(!tracklet_new.isValid())
		{
		  fitTracklet(tracklet_new);
		}
	      else
		{
		  continue;
		}
	   
#ifdef _DEBUG_ON
	      tracklet_new.print(); 
#endif
	      if(acceptTracklet(tracklet_new)) 
		{
		  trackletsInSt[listID].push_back(tracklet_new);
		}
#ifdef _DEBUG_ON
	      else
		{
		  LogInfo("Rejected!!!");
		}
#endif
	    }
	}
    } 

  //Reduce the tracklet list and add dummy hits
  //reduceTrackletList(trackletsInSt[listID]);
  for(std::list<Tracklet>::iterator iter = trackletsInSt[listID].begin(); iter != trackletsInSt[listID].end(); ++iter)
    {
      iter->addDummyHits();
    }
}

bool KalmanFastTracking::acceptTracklet(Tracklet& tracklet)
{
  //Tracklet itself is okay with enough hits (4-out-of-6) and small chi square
  if(!tracklet.isValid()) 
    {
#ifdef _DEBUG_ON
      LogInfo("Failed in quality check!");
#endif
      return false;
    }

  //LogInfo(tracklet.stationID);
  int nHodoHits = 0;
  for(std::vector<int>::iterator stationID = stationIDs_mask[tracklet.stationID-1].begin(); stationID != stationIDs_mask[tracklet.stationID-1].end(); ++stationID)
    {
      bool masked = false;
      for(std::list<int>::iterator iter = hitIDs_mask[*stationID-1].begin(); iter != hitIDs_mask[*stationID-1].end(); ++iter)
       	{
	  int detectorID = hitAll[*iter].detectorID;
	  int elementID = hitAll[*iter].elementID;
		  
	  int idx1 = detectorID - 25;
	  int idx2 = elementID - 1;
  
	  double z_hodo = z_mask[idx1];
    	  double x_hodo = tracklet.getExpPositionX(z_hodo);
	  double y_hodo = tracklet.getExpPositionY(z_hodo);
    	  double err_x = 3.*tracklet.getExpPosErrorX(z_hodo);
	  double err_y = 3.*tracklet.getExpPosErrorY(z_hodo);

	  double x_min = x_mask_min[idx1][idx2] - err_x;
	  double x_max = x_mask_max[idx1][idx2] + err_x;
	  double y_min = y_mask_min[idx1][idx2] - err_y;
	  double y_max = y_mask_max[idx1][idx2] + err_y;

#ifdef _DEBUG_ON
	  LogInfo(*iter);
	  hitAll[*iter].print();
	  LogInfo(nHodoHits << "/" << nMinimum << ":  " << z_hodo << "  " << x_hodo << " +/- " << err_x << "  " << y_hodo << " +/-" << err_y << " : " << x_min << "  " << x_max << "  " << y_min << "  " << y_max);
#endif
	  if(x_hodo > x_min && x_hodo < x_max && y_hodo > y_min && y_hodo < y_max)
    	    {
    	      nHodoHits++;
	      masked = true;

	      break;
    	    }
	}

      if(!masked) return false;
    }

#ifdef _DEBUG_ON
  LogInfo(tracklet.stationID << "  " << nHodoHits << "  " << stationIDs_mask[tracklet.stationID-1].size());
#endif

  //For back partials, require projection inside KMAG, and muon id in prop. tubes
  if(tracklet.stationID > 4)
    {
      if(!p_geomSvc->isInKMAG(tracklet.getExpPositionX(Z_KMAG_BEND), tracklet.getExpPositionY(Z_KMAG_BEND))) return false;
      if(!muonID(tracklet)) return false;
    }

  //If everything is fine ...
  return true;
}

bool KalmanFastTracking::muonID(Tracklet& tracklet)
{
  //Set the cut value on multiple scattering
  double cut = tracklet.stationID == 6 ? MUID_REJECT*(MUID_P0 + MUID_P1/tracklet.invP + MUID_P2/tracklet.invP/tracklet.invP) : 0.03; 

  double slope[2] = {tracklet.tx, tracklet.ty};
  PropSegment* segs[2] = {&(tracklet.seg_x), &(tracklet.seg_y)};
  for(int i = 0; i < 2; ++i)
    {
      segs[i]->init();
      for(int j = 0; j < 4; ++j)
	{
	  int index = detectorIDs_muid[i][j] - 25;
	  double x_exp = tracklet.getExpPositionX(z_mask[index]);
	  double y_exp = tracklet.getExpPositionY(z_mask[index]);
	  double pos_exp = p_geomSvc->getInterceptionFast(detectorIDs_muid[i][j], x_exp, y_exp);

	  if(!p_geomSvc->isInPlane(detectorIDs_muid[i][j], x_exp, y_exp)) continue;

	  double dist_min = 1E6;
	  for(std::list<int>::iterator iter = hitIDs_muid[i][j].begin(); iter != hitIDs_muid[i][j].end(); ++iter)
	    {
	      double pos = hitAll[*iter].pos;
	      if(fabs(pos - pos_exp) > 5.08) continue;

	      double dist_l = fabs(pos - hitAll[*iter].driftDistance - pos_exp);
	      double dist_r = fabs(pos + hitAll[*iter].driftDistance - pos_exp);
	      double dist = dist_l < dist_r ? dist_l : dist_r;
	      if(dist < dist_min)
		{
		  dist_min = dist;
		  if(dist < 2.54)
		    {
		      segs[i]->hits[j].hit = hitAll[*iter];
		      segs[i]->hits[j].sign = fabs(pos - hitAll[*iter].driftDistance - pos_exp) < fabs(pos + hitAll[*iter].driftDistance - pos_exp) ? -1 : 1;
		    }
		}
	    }
	}
      
      segs[i]->fit();
      if(!(segs[i]->isValid() && fabs(slope[i] - segs[i]->a) < cut)) return false;
    }

  return true;
}

int KalmanFastTracking::fitTracklet(Tracklet& tracklet)
{
  tracklet_curr = tracklet;

  //idx = 0, using simplex; idx = 1 using migrad
  int idx = 1;
#ifdef _ENABLE_MULTI_MINI
  if(tracklet.stationID < 5) idx = 0;
#endif

  minimizer[idx]->SetLimitedVariable(0, "tx", tracklet.tx, 0.001, -TX_MAX, TX_MAX);
  minimizer[idx]->SetLimitedVariable(1, "ty", tracklet.ty, 0.001, -TY_MAX, TY_MAX);
  minimizer[idx]->SetLimitedVariable(2, "x0", tracklet.x0, 0.1, -X0_MAX, X0_MAX);
  minimizer[idx]->SetLimitedVariable(3, "y0", tracklet.y0, 0.1, -Y0_MAX, Y0_MAX);
  if(KMAG_ON == 1)
    {
      minimizer[idx]->SetLimitedVariable(4, "invP", tracklet.invP, 0.001*tracklet.invP, INVP_MIN, INVP_MAX);
    }
  minimizer[idx]->Minimize();

  tracklet.tx = minimizer[idx]->X()[0];
  tracklet.ty = minimizer[idx]->X()[1];
  tracklet.x0 = minimizer[idx]->X()[2];
  tracklet.y0 = minimizer[idx]->X()[3];
 
  tracklet.err_tx = minimizer[idx]->Errors()[0];
  tracklet.err_ty = minimizer[idx]->Errors()[1];
  tracklet.err_x0 = minimizer[idx]->Errors()[2];
  tracklet.err_y0 = minimizer[idx]->Errors()[3];

  if(KMAG_ON == 1 && tracklet.stationID == 6)
    {
      tracklet.invP = minimizer[idx]->X()[4];
      tracklet.err_invP = minimizer[idx]->Errors()[4];
    }

  tracklet.chisq = minimizer[idx]->MinValue();
 
  int status = minimizer[idx]->Status();
  return status;
}

int KalmanFastTracking::reduceTrackletList(std::list<Tracklet>& tracklets)
{
  std::list<Tracklet> targetList;

  tracklets.sort();
  while(!tracklets.empty())
    {
      targetList.push_back(tracklets.front());
      tracklets.pop_front();

#ifdef _DEBUG_ON_LEVEL_2
      LogInfo("Current best tracklet in reduce");
      targetList.back().print();
#endif

      for(std::list<Tracklet>::iterator iter = tracklets.begin(); iter != tracklets.end(); )
	{
	  if(iter->similarity(targetList.back()))
	    {
#ifdef _DEBUG_ON_LEVEL_2
	      LogInfo("Removing this tracklet: ");
	      iter->print();
#endif
	      iter = tracklets.erase(iter);
	      continue;
	    }
	  else
	    {
	      ++iter;
	    }
	}
    }

  tracklets.assign(targetList.begin(), targetList.end());
  return 0;
}

void KalmanFastTracking::getExtrapoWindowsInSt1(Tracklet& tracklet, double* pos_exp, double* window)
{
  if(tracklet.stationID != 5)
    {
      for(int i = 0; i < 3; i++)
	{
	  pos_exp[i] = 9999.;
	  window[i] = 0.;
	}

      return;
    }

  for(int i = 0; i < 3; i++)
    {
      int detectorID = 2*i+2;
      double z_st1 = z_plane[detectorID];
      double x_st1 = tracklet.getExpPositionX(z_st1);
      double y_st1 = tracklet.getExpPositionY(z_st1);
      double err_x = tracklet.getExpPosErrorX(z_st1);
      double err_y = tracklet.getExpPosErrorY(z_st1);

      pos_exp[i] = p_geomSvc->getUinStereoPlane(detectorID, x_st1, y_st1);
      window[i] = 5.*(fabs(costheta_plane[detectorID]*err_x) + fabs(sintheta_plane[detectorID]*err_y));
    }
}

void KalmanFastTracking::getSagittaWindowsInSt1(Tracklet& tracklet, double* pos_exp, double* window)
{
  if(tracklet.stationID != 5)
    {
      for(int i = 0; i < 3; i++) 
	{
	  pos_exp[i] = 9999.;
	  window[i] = 0.;
	}

      return;
    }

  double z_st3 = z_plane[tracklet.hits.back().hit.detectorID];
  double x_st3 = tracklet.getExpPositionX(z_st3);
  double y_st3 = tracklet.getExpPositionY(z_st3);

  //For U, X, and V planes
  for(int i = 0; i < 3; i++)
    {
      double pos_st3 = p_geomSvc->getUinStereoPlane(s_detectorID[i], x_st3, y_st3);

      double z_st2 = z_plane[s_detectorID[i]];
      double x_st2 = tracklet.getExpPositionX(z_st2);
      double y_st2 = tracklet.getExpPositionY(z_st2);
      double pos_st2 = p_geomSvc->getUinStereoPlane(s_detectorID[i], x_st2, y_st2);
      double sagitta_st2 = pos_st2 - pos_st3*z_st2/z_st3;

      double z_st1 = z_plane[2*i+2];
      pos_exp[i] = sagitta_st2*s_ratio[i] + pos_st3*z_st1/z_st3;
      window[i] = fabs(5.*sagitta_st2*s_sigma[i]);
    }
}

void KalmanFastTracking::printAtDetectorBack(int stationID, std::string outputFileName)
{
  TCanvas c1;

  std::vector<double> x, y, dx, dy;
  for(std::list<Tracklet>::iterator iter = trackletsInSt[stationID].begin(); iter != trackletsInSt[stationID].end(); ++iter)
    {
      double z = p_geomSvc->getPlanePosition(iter->stationID*6);
      x.push_back(iter->getExpPositionX(z));
      y.push_back(iter->getExpPositionY(z));
      dx.push_back(iter->getExpPosErrorX(z));
      dy.push_back(iter->getExpPosErrorY(z));
    }

  TGraphErrors gr(x.size(), &x[0], &y[0], &dx[0], &dy[0]);
  gr.SetMarkerStyle(8);

  //Add detector frames
  std::vector<double> x_f, y_f, dx_f, dy_f;
  x_f.push_back(p_geomSvc->getPlaneCenterX(stationID*6 + 6));
  y_f.push_back(p_geomSvc->getPlaneCenterY(stationID*6 + 6));
  dx_f.push_back(p_geomSvc->getPlaneScaleX(stationID*6 + 6)*0.5);
  dy_f.push_back(p_geomSvc->getPlaneScaleY(stationID*6 + 6)*0.5);

  if(stationID == 2)
    {
      x_f.push_back(p_geomSvc->getPlaneCenterX(stationID*6 + 12));
      y_f.push_back(p_geomSvc->getPlaneCenterY(stationID*6 + 12));
      dx_f.push_back(p_geomSvc->getPlaneScaleX(stationID*6 + 12)*0.5);
      dy_f.push_back(p_geomSvc->getPlaneScaleY(stationID*6 + 12)*0.5);
    }
 
  TGraphErrors gr_frame(x_f.size(), &x_f[0], &y_f[0], &dx_f[0], &dy_f[0]); 
  gr_frame.SetLineColor(kRed);
  gr_frame.SetLineWidth(2);
  gr_frame.SetFillColor(15);
  
  c1.cd();
  gr_frame.Draw("A2[]");
  gr.Draw("Psame");

  c1.SaveAs(outputFileName.c_str());
}

void KalmanFastTracking::processOneTracklet(Tracklet& tracklet)
{
  //tracklet.print();
  KalmanTrack kmtrk;

  //Set the whole hit and node list
  for(std::list<SignedHit>::iterator iter = tracklet.hits.begin(); iter != tracklet.hits.end(); ++iter)
    {
      if(iter->hit.index < 0) continue;

      Node node_add(*iter);
      kmtrk.getNodeList().push_back(node_add);
      kmtrk.getHitIndexList().push_back(iter->sign*iter->hit.index);
    }

  //Set initial state
  TrkPar trkpar_curr;
  trkpar_curr._z = p_geomSvc->getPlanePosition(kmtrk.getNodeList().back().getHit().detectorID);
  trkpar_curr._state_kf[0][0] = tracklet.getCharge()*tracklet.invP/sqrt(1. + tracklet.tx*tracklet.tx + tracklet.ty*tracklet.ty);
  trkpar_curr._state_kf[1][0] = tracklet.tx;
  trkpar_curr._state_kf[2][0] = tracklet.ty;
  trkpar_curr._state_kf[3][0] = tracklet.getExpPositionX(trkpar_curr._z);
  trkpar_curr._state_kf[4][0] = tracklet.getExpPositionY(trkpar_curr._z);

  trkpar_curr._covar_kf.Zero();
  trkpar_curr._covar_kf[0][0] = 0.001;//1E6*tracklet.err_invP*tracklet.err_invP;
  trkpar_curr._covar_kf[1][1] = 0.01;//1E6*tracklet.err_tx*tracklet.err_tx;
  trkpar_curr._covar_kf[2][2] = 0.01;//1E6*tracklet.err_ty*tracklet.err_ty;
  trkpar_curr._covar_kf[3][3] = 100;//1E6*tracklet.getExpPosErrorX(trkpar_curr._z)*tracklet.getExpPosErrorX(trkpar_curr._z);
  trkpar_curr._covar_kf[4][4] = 100;//1E6*tracklet.getExpPosErrorY(trkpar_curr._z)*tracklet.getExpPosErrorY(trkpar_curr._z);

  kmtrk.setCurrTrkpar(trkpar_curr);
  kmtrk.getNodeList().back().getPredicted() = trkpar_curr;

  //Fit the track first with possibily a few nodes unresolved
  if(!fitTrack(kmtrk)) return;

  //Resolve left-right based on the current solution, re-fit if anything changed
  //resolveLeftRight(kmtrk);
  if(kmtrk.isValid()) tracks.push_back(kmtrk);
}

bool KalmanFastTracking::fitTrack(KalmanTrack& kmtrk)
{
  if(kmtrk.getNodeList().empty()) return false;

  if(kmfitter->processOneTrack(kmtrk) == 0)
    {
      return false;
    }
  kmfitter->updateTrack(kmtrk);

  return true;
}

void KalmanFastTracking::resolveLeftRight(KalmanTrack& kmtrk)
{
  bool isUpdated = false;

  std::list<int>::iterator hitID = kmtrk.getHitIndexList().begin();
  for(std::list<Node>::iterator node = kmtrk.getNodeList().begin(); node != kmtrk.getNodeList().end(); )
    {
      if(*hitID == 0)
	{
	  double x_hit = node->getSmoothed().get_x();
	  double y_hit = node->getSmoothed().get_y();
	  double pos_hit = p_geomSvc->getUinStereoPlane(node->getHit().detectorID, x_hit, y_hit);
	  
	  int sign = 0;
	  if(pos_hit > node->getHit().pos)
	    {
	      sign = 1;
	    }
	  else
	    {
	      sign = -1;
	    }

	  //update the node list
	  TMatrixD m(1, 1), dm(1, 1);
	  m[0][0] = node->getHit().pos + sign*node->getHit().driftDistance;
	  dm[0][0] = p_geomSvc->getPlaneResolution(node->getHit().detectorID)*p_geomSvc->getPlaneResolution(node->getHit().detectorID); 
	  node->setMeasurement(m, dm);
	  *hitID = sign*node->getHit().index;

	  isUpdated = true;
	}

      ++node;
      ++hitID;
    }

  if(isUpdated) fitTrack(kmtrk);
}

void KalmanFastTracking::chi2fit(int n, double x[], double y[], double& a, double& b)
{
  double sum = 0.;
  double sx = 0.;
  double sy = 0.;
  double sxx = 0.;
  double syy = 0.;
  double sxy = 0.;

  for(int i = 0; i < n; ++i)
    {
      ++sum;
      sx += x[i];
      sy += y[i];
      sxx += (x[i]*x[i]);
      syy += (y[i]*y[i]);
      sxy += (x[i]*y[i]);
    }

  double det = sum*sxx - sx*sx;
  if(fabs(det) < 1E-20) 
    {
      a = 0.;
      b = 0.;

      return;
    }

  a = (sum*sxy - sx*sy)/det;
  b = (sy*sxx - sxy*sx)/det;
}
