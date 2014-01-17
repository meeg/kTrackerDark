/*
VertexFit.cxx

Implementation of the VertexFit, in the beginning only a old-fashioned primary vertex finder
is implemented

Author: Kun Liu, liuk@fnal.gov
Created: 2-8-2012
*/

#include <iostream>
#include <cmath>
#include <TMatrixD.h>

#include "VertexFit.h"

VertexFit::VertexFit()
{
  ///In construction, initialize the projector for the vertex node
  TMatrixD proj(2, 5);
  proj.Zero();

  proj[0][3] = 1.;
  proj[1][4] = 1.;

  _node_vertex.getMeasurement().ResizeTo(2, 1);
  _node_vertex.getMeasurementCov().ResizeTo(2, 2);
  _node_vertex.getProjector().ResizeTo(2, 5);

  _node_vertex.getProjector() = proj;

  _node_vertex.getMeasurementCov().Zero();
  _node_vertex.getMeasurementCov()[0][0] = 1.;
  _node_vertex.getMeasurementCov()[1][1] = 1.;

  _max_iteration = 100;
  _tolerance = 1E-3;

  _kmfit = KalmanFilter::instance();
  _extrapolator.init(GEOMETRY_VERSION, true);
  
  ///Single track finding doesn't require a propagation matrix
  _extrapolator.setPropCalc(false);
}

VertexFit::~VertexFit()
{
}

bool VertexFit::setRecEvent(SRecEvent* recEvent)
{
  //if the single vertex is not set, set it first
  int nTracks = recEvent->getNTracks();
  for(int i = 0; i < nTracks; ++i)
    {
      SRecTrack& recTrack = recEvent->getTrack(i);
      if(recTrack.getChisqVertex() < 0.)
	{
	  recTrack.setZVertex(findSingleMuonVertex(recTrack));
	}
    }

  std::vector<int> idx_pos = recEvent->getChargedTrackIDs(+1);
  std::vector<int> idx_neg = recEvent->getChargedTrackIDs(-1);

  int nPos = idx_pos.size();
  int nNeg = idx_neg.size();
  if(nPos*nNeg == 0) return false;

  for(int i = 0; i < nPos; ++i)
    {
      SRecTrack track_pos = recEvent->getTrack(idx_pos[i]);
      if(!track_pos.isValid()) continue;
      for(int j = 0; j < nNeg; ++j)
	{
	  SRecTrack track_neg = recEvent->getTrack(idx_neg[j]);
	  if(!track_neg.isValid()) continue;

	  SRecDimuon dimuon;
	  dimuon.trackID_pos = idx_pos[i];
	  dimuon.trackID_neg = idx_neg[j];

	  dimuon.p_pos_single = track_pos.getMomentumVertex();
	  dimuon.p_neg_single = track_neg.getMomentumVertex();
	  dimuon.vtx_pos = track_pos.getVertex();
	  dimuon.vtx_neg = track_neg.getVertex();

	  //Start prepare the vertex fit
	  init();
	  addTrack(0, track_pos);
	  addTrack(1, track_neg);
	  addHypothesis(0.5*(dimuon.vtx_pos[2] + dimuon.vtx_neg[2]), 50.);
	  processOnePair();

	  //Retrieve the results
	  track_pos.setZVertex(getVertexZ0());
	  track_neg.setZVertex(getVertexZ0());
	  dimuon.p_pos = track_pos.getMomentumVertex();
	  dimuon.p_neg = track_neg.getMomentumVertex();
	  dimuon.chisq_kf = getKFChisq();
	  dimuon.chisq_vx = getVXChisq();
	  dimuon.vtx.SetXYZ(_vtxpar_curr._r[0][0], _vtxpar_curr._r[1][0], _vtxpar_curr._r[2][0]);
	  dimuon.calcVariables();

	  recEvent->insertDimuon(dimuon);
	}
    }

  if(recEvent->getNDimuons() > 0) return true;
  return false; 
}

void VertexFit::init()
{
  _trkpar_curr.clear();

  z_vertex.clear();
  chisq_km.clear();
  chisq_vx.clear();

  z_start.clear();
  sig_z_start.clear();

  ///Two default starting points
  addHypothesis(100., 50.);
  addHypothesis(-130., 50.);
}

void VertexFit::setStartingVertex(double z_start, double sigz_start)
{
  ///Initialize the starting vertex with a guess and large error 
  _vtxpar_curr._r.Zero();
  _vtxpar_curr._r[2][0] = z_start;

  _vtxpar_curr._cov.Zero();
  _vtxpar_curr._cov[0][0] = 1.5;
  _vtxpar_curr._cov[1][1] = 1.5;
  _vtxpar_curr._cov[2][2] = sigz_start*sigz_start;
  
  _chisq_vertex = 0.;
  _chisq_kalman = 0.;
}

int VertexFit::processOnePair()
{
  double chisq_min = 1E6;
  int index_min = -1; 
  for(unsigned int i = 0; i < z_start.size(); i++)
    {
#ifdef _DEBUG_ON
      LogInfo("Testing starting point: " << z_start[i]);
#endif

      setStartingVertex(z_start[i], sig_z_start[i]);
      findVertex();

      chisq_km.push_back(_chisq_kalman);
      chisq_vx.push_back(_chisq_vertex);
      z_vertex.push_back(_vtxpar_curr._r[2][0]);

      if(_chisq_kalman < chisq_min)
	{
	  index_min = i;
	  chisq_min = _chisq_kalman;
	}
    }

  if(index_min < 0) return 0;

  _chisq_kalman = chisq_km[index_min];
  _chisq_vertex = chisq_vx[index_min];
  _vtxpar_curr._r[2][0] = z_vertex[index_min];
  if(z_vertex[index_min] > 1E4) _vtxpar_curr._r[2][0] = z_start.back();

  return index_min+1;
}

void VertexFit::addTrack(int index, KalmanTrack& _track)
{
  if(_track.getNodeList().front().isSmoothDone())
    {
      _trkpar_curr.push_back(_track.getNodeList().front().getSmoothed());
    }
  else if(_track.getNodeList().front().isFilterDone())
    {
      _trkpar_curr.push_back(_track.getNodeList().front().getFiltered());
    }
  else
    {
      _trkpar_curr.push_back(_track.getNodeList().front().getPredicted());
    }
}

void VertexFit::addTrack(int index, SRecTrack& _track)
{
  TrkPar _trkpar;
  _trkpar._state_kf = _track.getStateVector(0);
  _trkpar._covar_kf = _track.getCovariance(0);
  _trkpar._z = _track.getZ(0);

  _trkpar_curr.push_back(_trkpar);
}

void VertexFit::addTrack(int index, TrkPar& _trkpar)
{
  _trkpar_curr.push_back(_trkpar);  
}

int VertexFit::findVertex()
{
  int nIter = 0;
  for(; nIter < _max_iteration; nIter++)
    {
#ifdef _DEBUG_ON
      LogInfo("Iteration: " << nIter);
#endif

      _chisq_vertex = 0.;
      _chisq_kalman = 0.;
      for(unsigned int j = 0; j < _trkpar_curr.size(); j++)
	{
	  _node_vertex.resetFlags();
	  //_node_vertex.getMeasurement() = _vtxpar_curr._r.GetSub(0, 1, 0, 0);
	  //_node_vertex.getMeasurementCov() = _vtxpar_curr._cov.GetSub(0, 1, 0, 1);
	  _node_vertex.setZ(_vtxpar_curr._r[2][0]);

	  _kmfit->setCurrTrkpar(_trkpar_curr[j]);
	  if(!_kmfit->fit_node(_node_vertex))
	    {
#ifdef _DEBUG_ON
	      LogInfo("Vertex fit for this track failed!");
#endif
	      _chisq_kalman = 1E5;	    
	      break;
	    }
	  else
	    {
	      _chisq_kalman += _node_vertex.getChisq();
	    }

	  updateVertex();
  	}

      ///break the iteration if the z0 converges
#ifdef _DEBUG_ON
      LogInfo("At this iteration: ");
      LogInfo(_vtxpar_curr._r[2][0] << " ===  " << _node_vertex.getZ() << " : " << _chisq_kalman << " === " << _chisq_vertex << " : " << _vtxpar_curr._r[0][0] << " === " << _vtxpar_curr._r[1][0]);
#endif
      if(_vtxpar_curr._r[2][0] < -260. || _vtxpar_curr._r[2][0] > 700.)
	{
	  _chisq_kalman = 1E5;
	  _chisq_vertex = 1E5;
	  break;
	} 

      if(nIter > 0 && fabs(_vtxpar_curr._r[2][0] - _node_vertex.getZ()) < _tolerance) break;
    } 

  return nIter+1;
}

void VertexFit::updateVertex()
{
  double p = fabs(1./_node_vertex.getFiltered()._state_kf[0][0]);
  double px = _node_vertex.getFiltered()._state_kf[1][0];
  double py = _node_vertex.getFiltered()._state_kf[2][0];
  double pz = sqrt(p*p - px*px - py*py);

  ///Set the projector matrix from track state vector to the coordinate
  TMatrixD H(2, 3);
  H.Zero();

  H[0][0] = 1.;
  H[1][1] = 1.;
  H[0][2] = -px/pz;
  H[1][2] = -py/pz;

  TMatrixD vertex_dummy(3, 1);
  vertex_dummy.Zero();
  vertex_dummy[2][0] = _vtxpar_curr._r[2][0];
  
  TMatrixD mxy = _node_vertex.getFiltered()._state_kf.GetSub(3, 4, 0, 0);
  TMatrixD Vxy = _node_vertex.getFiltered()._covar_kf.GetSub(3, 4, 3, 4);
  TMatrixD S = SMatrix::invertMatrix(Vxy + SMatrix::getABCt(H, _vtxpar_curr._cov, H));
  TMatrixD K = SMatrix::getABtC(_vtxpar_curr._cov, H, S);
  TMatrixD zeta = mxy - H*(_vtxpar_curr._r - vertex_dummy);
  TMatrixD _r_filtered = _vtxpar_curr._r + K*zeta;
  TMatrixD _cov_filtered = _vtxpar_curr._cov - K*H*(_vtxpar_curr._cov);

  _chisq_vertex += SMatrix::getAtBC(zeta, S, zeta)[0][0];

  _vtxpar_curr._r = _r_filtered;
  _vtxpar_curr._cov = _cov_filtered;
}

double VertexFit::findSingleMuonVertex(SRecTrack& _track)
{
  TrkPar _trkpar_start;
  _trkpar_start._state_kf = _track.getStateVector(0);
  _trkpar_start._covar_kf = _track.getCovariance(0);
  _trkpar_start._z = _track.getZ(0);

  return findSingleMuonVertex(_trkpar_start);
}

double VertexFit::findSingleMuonVertex(Node& _node_start)
{
  TrkPar _trkpar_start;
  if(_node_start.isSmoothDone())
    {
      _trkpar_start = _node_start.getSmoothed();
    }
  else if(_node_start.isFilterDone())
    {
      _trkpar_start = _node_start.getFiltered();
    }
  else
    {
      _trkpar_start = _node_start.getPredicted();
    }

  return findSingleMuonVertex(_trkpar_start);
}

double VertexFit::findSingleMuonVertex(TrkPar& _trkpar_start)
{
  _extrapolator.setInitialStateWithCov(_trkpar_start._z, _trkpar_start._state_kf, _trkpar_start._covar_kf);
  double z_vertex = _extrapolator.extrapolateToIP(-250., 5.);
  
  return z_vertex;
}

void VertexFit::print()
{
  using namespace std;

  for(unsigned int i = 0; i < z_start.size(); i++)
    {
      cout << "============= Hypothesis " << i << " ============" << endl;
      cout << "z_start = " << z_start[i] << ", sigz_start = " << sig_z_start[i] << endl;
      cout << "Found z_vertex = " << z_vertex[i] << endl;
      cout << "With chisq_km = " << chisq_km[i] << " and chisq_vx = " << chisq_vx[i] << endl; 
    }
}
