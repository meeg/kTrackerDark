#include <iostream>
#include <cmath>
#include <algorithm>
#include <string>
#include <time.h>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TRandom.h>
#include <TMatrixD.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <TRandom.h>

#include "GeomSvc.h"
#include "SRawEvent.h"
#include "KalmanUtil.h"
#include "KalmanTrack.h"
#include "KalmanFilter.h"
#include "KalmanFinder.h"
#include "KalmanFitter.h"
#include "VertexFit.h"
#include "SRecEvent.h"

TLorentzVector getMom(double px, double py, double pz)
{
  double mass_mu = 0.10566;
  double E = sqrt(px*px + py*py + pz*pz + mass_mu*mass_mu);
  
  TLorentzVector p;
  p.SetPxPyPzE(px, py, pz, E);
      
  return p;
}

void calc_variables(TLorentzVector& p1, TLorentzVector& p2, double& mass, double& pT, double& xF, double& x1, double& x2)
{
  double mp = 0.938;
  double ebeam = 120.;

  TLorentzVector p_beam(0., 0., sqrt(ebeam*ebeam - mp*mp), ebeam);
  TLorentzVector p_target(0., 0., 0., mp);

  TLorentzVector p_cms = p_beam + p_target;
  TVector3 bv_cms = p_cms.BoostVector();
  double s = p_cms.M2();

  TLorentzVector p_sum = p1 + p2;
  mass = p_sum.M();
  pT = p_sum.Perp();

  p_sum.Boost(-bv_cms);
  xF = 2*p_sum.Pz()/TMath::Sqrt(s);
  double tau = p_sum.M2()/s;
  double y = 0.5*std::log((p_sum.E() + p_sum.Pz())/(p_sum.E() - p_sum.Pz()));

  x1 = TMath::Sqrt(tau)*TMath::Exp(y);
  x2 = TMath::Sqrt(tau)*TMath::Exp(-y);
}

bool trackOK(SRecTrack& _track)
{
  if(!_track.isVertexValid()) return false;
  if(_track.getChisq() > 5.) return false;
  if(_track.getChisqVertex() > 15.) return false;
  if(_track.getNHits() < 12) return false;
  if(!_track.isHodoMasked()) return false;
  return true;
}

int main(int argc, char *argv[])
{
  //Initialize geometry service
  Log("Initializing geometry service ... ");
  GeomSvc* geometrySvc = GeomSvc::instance();
  geometrySvc->init(GEOMETRY_VERSION);

  //Retrieve the raw event
  Log("Retrieving the event stored in ROOT file ... ");
#ifdef MC_MODE
  SRawMCEvent* rawEvent = new SRawMCEvent();
#else
  SRawEvent* rawEvent = new SRawEvent();
#endif
  SRecEvent* recEvent = new SRecEvent();

  TFile *dataFile = new TFile(argv[1], "READ");
  TTree *dataTree = (TTree *)dataFile->Get("save");

  dataTree->SetBranchAddress("rawEvent", &rawEvent);
  dataTree->SetBranchAddress("recEvent", &recEvent);

  double chisq_pair;
  double chisq_single;
  double z_pair;
  double z_single;
  double dz_single;

  double p1_pair, p1x_pair, p1y_pair, p1z_pair;
  double p2_pair, p2x_pair, p2y_pair, p2z_pair;

  double mass_pair;
  double pT, xF, x1, x2;

  TFile *saveFile = new TFile(argv[2], "recreate");
  TTree *saveTree = new TTree("save", "save");

  saveTree->Branch("z_pair", &z_pair, "z_pair/D");
  saveTree->Branch("z_single", &z_single, "z_single/D");
  saveTree->Branch("dz_single", &dz_single, "dz_single/D");

  saveTree->Branch("chisq_pair", &chisq_pair, "chisq_pair/D");
  saveTree->Branch("chisq_single", &chisq_single, "chisq_single/D");
  saveTree->Branch("p1_pair", &p1_pair, "p1_pair/D");
  saveTree->Branch("p1x_pair", &p1x_pair, "p1x_pair/D");
  saveTree->Branch("p1y_pair", &p1y_pair, "p1y_pair/D");
  saveTree->Branch("p1z_pair", &p1z_pair, "p1z_pair/D");
  saveTree->Branch("p2_pair", &p2_pair, "p2_pair/D");
  saveTree->Branch("p2x_pair", &p2x_pair, "p2x_pair/D");
  saveTree->Branch("p2y_pair", &p2y_pair, "p2y_pair/D");
  saveTree->Branch("p2z_pair", &p2z_pair, "p2z_pair/D");
  
  saveTree->Branch("mass_pair", &mass_pair, "mass_pair/D");
  saveTree->Branch("pT", &pT, "pT/D");
  saveTree->Branch("xF", &xF, "xF/D");
  saveTree->Branch("x1", &x1, "x1/D");
  saveTree->Branch("x2", &x2, "x2/D");
  
  //Initialize track finder
  Log("Initializing the track finder and kalman filter ... ");
  VertexFit *vtxfit = new VertexFit();

  vector<SRecTrack> ptracks, mtracks;
  vector<int> pflags, mflags;
  ptracks.clear(); mtracks.clear();
  pflags.clear(); mflags.clear();
  int nEvtMax = dataTree->GetEntries();
  for(int i = 0; i < nEvtMax; i++)
    {
      dataTree->GetEntry(i);
      rawEvent->reIndex("oah");
      
      int nTracks = recEvent->getNTracks();
      for(int j = 0; j < nTracks; j++)
	{
	  SRecTrack& _track = recEvent->getTrack(j);
	  _track.setZVertex(vtxfit->findSingleMuonVertex(_track));

	  if(trackOK(_track)) 
	    {
	      if(_track.getCharge() > 0)
		{
		  ptracks.push_back(_track);
		  pflags.push_back(1);
		}
	      else
		{
		  mtracks.push_back(_track);
		  mflags.push_back(1);
		}
	    }
	}

      rawEvent->clear();
      recEvent->clear();
    }

  TRandom rnd;
  int nPlus = ptracks.size();
  int nMinus = mtracks.size();
  int nEntries = atoi(argv[3]);
  cout << nPlus << "  " << nMinus << endl;
  while(saveTree->GetEntries() < nEntries)
    {
      cout << saveTree->GetEntries() << endl;

      int id1 = int(rnd.Rndm()*nPlus);
      int id2 = int(rnd.Rndm()*nMinus);

      //Neither of the tracks should be used
      if(pflags[id1] < 0 || mflags[id2] < 0) continue;
      
      //Total momentum constrain
      double px1, py1, pz1;
      double px2, py2, pz2;
      ptracks[id1].getMomentumVertex(px1, py1, pz1);
      mtracks[id1].getMomentumVertex(px2, py2, pz2);
      if(pz1 + pz1 > 120.) continue;

      //Z separation constrain
      dz_single = ptracks[id1].getZVertex() - mtracks[id2].getZVertex();
      z_single = (ptracks[id1].getZVertex() + mtracks[id2].getZVertex())/2.;
      if(fabs(dz_single) > 50.) continue;

      vtxfit->init();
      vtxfit->addTrack(0, ptracks[id1]);
      vtxfit->addTrack(1, mtracks[id2]);

      vtxfit->addHypothesis(z_single, 50);
      vtxfit->processOneEvent();

      z_pair = vtxfit->getVertexZ0();
      chisq_pair = vtxfit->getKFChisq();

      ptracks[id1].setZVertex(z_pair);
      mtracks[id2].setZVertex(z_pair);

      p1_pair = ptracks[id1].getMomentumVertex(p1x_pair, p1y_pair, p1z_pair);
      p2_pair = mtracks[id2].getMomentumVertex(p2x_pair, p2y_pair, p2z_pair);

      TLorentzVector p1_4mom = getMom(p1x_pair, p1y_pair, p1z_pair);
      TLorentzVector p2_4mom = getMom(p2x_pair, p2y_pair, p2z_pair);
      calc_variables(p1_4mom, p2_4mom, mass_pair, pT, xF, x1, x2);

      cout << id1 << "  " << id2 << endl;
      pflags[id1] = -1;
      mflags[id2] = -1;
      saveTree->Fill();
    }

  saveFile->cd();
  saveTree->Write();
  saveFile->Close();

  delete vtxfit;

  return 1;
}
