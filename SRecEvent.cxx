/*
SRecEvent.cxx

Implimentation of the calss SRecTrack and SRecEvent

Author: Kun Liu, liuk@fnal.gov
Created: 01-21-2013
*/

#include <cmath>

#include <TLorentzVector.h>
#include <TMatrixD.h>

#include "SRecEvent.h"
#include "KalmanUtil.h"
#include "KalmanFilter.h"

ClassImp(SRecTrack)
ClassImp(SRecDimuon)
ClassImp(SRecEvent)

SRecTrack::SRecTrack()
{
  fChisq = -99.;

  fHitIndex.clear();
  fState.clear();
  fCovar.clear();
  fZ.clear();
  fChisqAtNode.clear();

  fChisqVertex = -99.;
  fVtxPar[0] = 999.;
  fVtxPar[1] = 999.;
  fVtxPar[2] = 999.;
  fStateVertex.ResizeTo(5, 1);
  fCovarVertex.ResizeTo(5, 5);
}

bool SRecTrack::operator<(const SRecTrack& elem) const
{
  if((Double_t)fHitIndex.size() - 0.4*fChisq > (Double_t)elem.fHitIndex.size() - 0.4*elem.fChisq)
    {
      return true;
    }
  else
    {
      return false;
    }
}

void SRecTrack::setZVertex(Double_t z)
{
  Node _node_vertex;
  _node_vertex.setZ(z);

  TMatrixD m(2, 1), cov(2, 2),  proj(2, 5);
  m[0][0] = 0.;
  m[1][0] = 0.;

  cov.Zero();
  cov[0][0] = BEAM_SPOT_X*BEAM_SPOT_X;
  cov[1][1] = BEAM_SPOT_Y*BEAM_SPOT_Y;

  proj.Zero();
  proj[0][3] = 1.;
  proj[1][4] = 1.;

  _node_vertex.getMeasurement().ResizeTo(2, 1);
  _node_vertex.getMeasurementCov().ResizeTo(2, 2);
  _node_vertex.getProjector().ResizeTo(2, 5);

  _node_vertex.getMeasurement() = m;
  _node_vertex.getMeasurementCov() = cov;
  _node_vertex.getProjector() = proj;

  TrkPar _trkpar_curr;
  _trkpar_curr._state_kf = fState[0];
  _trkpar_curr._covar_kf = fCovar[0];
  _trkpar_curr._z = fZ[0];

  KalmanFilter *kmfit = KalmanFilter::instance();
  kmfit->setCurrTrkpar(_trkpar_curr);
  kmfit->fit_node(_node_vertex);

  fChisqVertex = _node_vertex.getChisq();
  fVtxPar[0] = _node_vertex.getFiltered().get_x();
  fVtxPar[1] = _node_vertex.getFiltered().get_y();
  fVtxPar[2] = z;

  fStateVertex = _node_vertex.getFiltered()._state_kf;
  fCovarVertex = _node_vertex.getFiltered()._covar_kf;
}

void SRecTrack::setVertexFast(TVector3 mom, TVector3 pos)
{
  fVtxPar[0] = pos[0];
  fVtxPar[1] = pos[1];
  fVtxPar[2] = pos[2];

  fStateVertex[0][0] = getCharge()/mom.Mag();
  fStateVertex[1][0] = mom[0]/mom[2];
  fStateVertex[2][0] = mom[1]/mom[2];
  fStateVertex[3][0] = pos[0];
  fStateVertex[4][0] = pos[1];

  fCovarVertex.UnitMatrix();
}

bool SRecTrack::isVertexValid()
{
  if(fChisqVertex > 50.) return false;

  return true;
}

bool SRecTrack::isHodoMasked()
{
  for(Int_t i = 0; i < 3; i++)
    {
      if(fNHodoHits[i] < 1) return false;
    }

  return true;
}

Int_t SRecTrack::getNearestNode(Double_t z)
{
  Int_t nNodes = getNHits();
  Double_t deltaZ_curr = 0;
  Double_t deltaZ_prev = 1E6;
  for(Int_t i = 0; i < nNodes; i++)
    {
      deltaZ_curr = fabs(z - fZ[i]);
      if(deltaZ_curr > deltaZ_prev)
	{
	  return i - 1;
	}

      deltaZ_prev = deltaZ_curr;
    }

  return nNodes - 1;
}

void SRecTrack::getExpPositionFast(Double_t z, Double_t& x, Double_t& y, Int_t iNode)
{
  if(iNode < 0 || iNode >= getNHits())
    {
      iNode = getNearestNode(z);
    }

  TrkPar _trkpar;
  _trkpar._state_kf = fState[iNode];
  _trkpar._covar_kf = fCovar[iNode];

  Double_t z_ref = fZ[iNode];
  Double_t x_ref = _trkpar.get_x();
  Double_t y_ref = _trkpar.get_y();
  Double_t axz = _trkpar.get_dxdz();
  Double_t ayz = _trkpar.get_dydz();

  x = x_ref + axz*(z - z_ref);
  y = y_ref + ayz*(z - z_ref);
}

void SRecTrack::getExpPosErrorFast(Double_t z, Double_t& dx, Double_t& dy, Int_t iNode)
{
  if(iNode < 0 || iNode >= getNHits())
    {
      iNode = getNearestNode(z);
    }

  Double_t z_ref = fZ[iNode];
  Double_t dx_ref = sqrt((fCovar[iNode])[3][3]);
  Double_t dy_ref = sqrt((fCovar[iNode])[4][4]);
  Double_t daxz = sqrt((fCovar[iNode])[1][1]);
  Double_t dayz = sqrt((fCovar[iNode])[2][2]);

  dx = 2.*(dx_ref + fabs(daxz*(z - z_ref)));
  dy = 2.*(dy_ref + fabs(dayz*(z - z_ref)));
}

double SRecTrack::getExpMomentumFast(Double_t z, Int_t iNode)
{
  Double_t px, py, pz;
  return getExpMomentumFast(z, px, py, pz, iNode);
}

Double_t SRecTrack::getExpMomentumFast(Double_t z, Double_t& px, Double_t& py, Double_t& pz, Int_t iNode)
{
  if(iNode < 0 || iNode >= getNHits())
    {
      iNode = getNearestNode(z);
    }

   return getMomentum(fState[iNode], px, py, pz);

}

Double_t SRecTrack::getMomentum(TMatrixD& state, Double_t& px, Double_t& py, Double_t& pz)
{
  Double_t p = 1./fabs(state[0][0]);
  pz = p/sqrt(1. + state[1][0]*state[1][0] + state[2][0]*state[2][0]);
  px = p*state[1][0];
  py = p*state[2][0];

  return p;
}

Double_t SRecTrack::getPosition(TMatrixD& state, Double_t& x, Double_t& y)
{
  x = state[3][0];
  y = state[4][0];

  return sqrt(x*x + y*y);
}

TLorentzVector SRecTrack::getMomentumVertex()
{
  Double_t mmu = 0.10566;
  Double_t px, py, pz, E;

  getMomentumVertex(px, py, pz);
  E = sqrt(px*px + py*py + pz*pz + mmu*mmu);

  return TLorentzVector(px, py, pz, E);
}

bool SRecTrack::isValid()
{
  //Vertex valid
  if(!isVertexValid()) return false;

  //Number of hits cut
  if(getNHits() < 12) return false;

  //Total chisq, may change to cut on prob
  if(getChisq() > 20.) return false;

  //hodo and prop. tube masking
  if(!isHodoMasked()) return false;

  return true;
}

void SRecTrack::print()
{
  std::cout << "=============== Reconstructed track ==================" << std::endl;
  std::cout << "This candidate has " << fHitIndex.size() << " hits!" << std::endl;
  std::cout << "Most upstream momentum is: " << 1./fabs((fState[0])[0][0]) << std::endl;
  std::cout << "Chi square of the track is: " << fChisq << std::endl;

  std::cout << "Hodoscope hits: " << std::endl;
  for(Int_t i = 0; i < 3; i++) std::cout << "Station " << i+1 << ": " << fNHodoHits[i] << std::endl;

  std::cout << "Current vertex position: " << std::endl;
  for(Int_t i = 0; i < 3; i++) std::cout << fVtxPar[i] << "  ";
  std::cout << std::endl;

  std::cout << "Momentum at vertex: " << 1./fabs(fStateVertex[0][0]) << std::endl; 
  std::cout << "Chi square at vertex: " << fChisqVertex << std::endl;
}

void SRecDimuon::calcVariables()
{
  Double_t mp = 0.938;
  Double_t ebeam = 120.;

  TLorentzVector p_beam(0., 0., sqrt(ebeam*ebeam - mp*mp), ebeam);
  TLorentzVector p_target(0., 0., 0., mp);

  TLorentzVector p_cms = p_beam + p_target;
  TVector3 bv_cms = p_cms.BoostVector();
  Double_t s = p_cms.M2();

  TLorentzVector p_sum = p_pos + p_neg;
  mass = p_sum.M();
  pT = p_sum.Perp();

  p_sum.Boost(-bv_cms);
  xF = 2.*p_sum.Pz()/TMath::Sqrt(s);
  costh = p_sum.CosTheta();
  Double_t tau = p_sum.M2()/s;
  Double_t y = 0.5*TMath::Log((p_sum.E() + p_sum.Pz())/(p_sum.E() - p_sum.Pz()));
  x1 = TMath::Sqrt(tau)*TMath::Exp(y);
  x2 = TMath::Sqrt(tau)*TMath::Exp(-y);
}

SRecEvent::SRecEvent()
{
  fRunID = -1;
  fSpillID = -1;
  fEventID = -1;

  clear();
}

void SRecEvent::setRawEvent(SRawEvent *rawEvent)
{
  setEventInfo(rawEvent->getRunID(), rawEvent->getSpillID(), rawEvent->getEventID());
  fRunID = rawEvent->getRunID();
  fSpillID = rawEvent->getSpillID();
  fEventID = rawEvent->getEventID();

  for(Int_t i = 0; i < rawEvent->getNHitsAll(); i++)
    {
      fLocalID.insert(std::map<Int_t, Int_t>::value_type(rawEvent->getHit(i).index, i));
    }

  fAllTracks.clear();
}

void SRecEvent::setEventInfo(Int_t runID, Int_t spillID, Int_t eventID)
{
  fRunID = runID;
  fSpillID = spillID;
  fEventID = eventID;
}

std::vector<Int_t> SRecEvent::getChargedTrackIDs(Int_t charge)
{
  std::vector<Int_t> trkIDs;
  trkIDs.clear();

  Int_t nTracks = getNTracks();
  for(Int_t i = 0; i < nTracks; i++)
    {
      if(fAllTracks[i].getCharge() == charge)
	{
	  trkIDs.push_back(i);
	}
    }

  return trkIDs;
}

void SRecEvent::clear()
{
  fAllTracks.clear();
  fLocalID.clear();
  fDimuons.clear();
}
