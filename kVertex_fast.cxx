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

#include "GeomSvc.h"
#include "SRawEvent.h"
#include "KalmanUtil.h"
#include "KalmanTrack.h"
#include "KalmanFilter.h"
#include "KalmanFinder.h"
#include "KalmanFitter.h"
#include "VertexFit.h"
#include "SRecEvent.h"
#include "FastTracklet.h"

TLorentzVector getMom(double px, double py, double pz)
{
  double mass_mu = 0.10566;
  double E = sqrt(px*px + py*py + pz*pz + mass_mu*mass_mu);
  
  TLorentzVector p;
  p.SetPxPyPzE(px, py, pz, E);
      
  return p;
}

void calc_variables(TLorentzVector p1, TLorentzVector p2, double& mass, double& pT, double& xF, double& x1, double& x2)
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

int main(int argc, char *argv[])
{
  //Initialize geometry service
  Log("Initializing geometry service ... ");
  GeomSvc* geometrySvc = GeomSvc::instance();
  geometrySvc->init(GEOMETRY_VERSION);

  //Retrieve the raw event
  Log("Retrieving the event stored in ROOT file ... ");

  SRawEvent *rawEvent = new SRawEvent();
  TClonesArray* tracklets = new TClonesArray("Tracklet");

  TFile *dataFile = new TFile(argv[1], "READ");
  TTree *dataTree = (TTree *)dataFile->Get("save");

  dataTree->SetBranchAddress("rawEvent", &rawEvent);
  dataTree->SetBranchAddress("tracklets", &tracklets);

  int nDimuons;
  int pID[30], mID[30];
  double p1_pair[30], p1x_pair[30], p1y_pair[30], p1z_pair[30];
  double p1_single[30], p1x_single[30], p1y_single[30], p1z_single[30];
  double p2_pair[30], p2x_pair[30], p2y_pair[30], p2z_pair[30];
  double p2_single[30], p2x_single[30], p2y_single[30], p2z_single[30];

  double r1_single[30], r2_single[30];
  double r1_pair[30], r2_pair[30];
  double angle_pair[30];

  double chisq_vx[30];
  double chisq_pair[30];
  double chisq_single[30];
  double z_pair[30];
  double z_single[30];
  double dz_single[30];

  int hypoID[30];
  double mass_pair[30];
  double mass_single[30];
  double pT_pair[30], xF_pair[30], x1_pair[30], x2_pair[30];
  double pT_single[30], xF_single[30], x1_single[30], x2_single[30];

  TFile *saveFile = new TFile(argv[2], "recreate");
  TTree *saveTree = new TTree("save", "save");

  saveTree->Branch("rawEvent", &rawEvent, 256000, 99);
  saveTree->Branch("tracklets", &tracklets, 256000, 99);
  
  saveTree->Branch("nDimuons", &nDimuons, "nDimuons/I");
  saveTree->Branch("hypoID", hypoID, "hypoID[nDimuons]/I");
  saveTree->Branch("pID", pID, "pID[nDimuons]/I");
  saveTree->Branch("mID", mID, "mID[nDimuons]/I");
 
  saveTree->Branch("z_pair", z_pair, "z_pair[nDimuons]/D");
  saveTree->Branch("z_single", z_single, "z_single[nDimuons]/D");
  saveTree->Branch("dz_single", dz_single, "dz_single[nDimuons]/D");

  saveTree->Branch("chisq_vx", chisq_vx, "chisq_vx[nDimuons]/D");
  saveTree->Branch("chisq_pair", chisq_pair, "chisq_pair[nDimuons]/D");
  saveTree->Branch("chisq_single", chisq_single, "chisq_single[nDimuons]/D");
  saveTree->Branch("p1_pair", p1_pair, "p1_pair[nDimuons]/D");
  saveTree->Branch("p1x_pair", p1x_pair, "p1x_pair[nDimuons]/D");
  saveTree->Branch("p1y_pair", p1y_pair, "p1y_pair[nDimuons]/D");
  saveTree->Branch("p1z_pair", p1z_pair, "p1z_pair[nDimuons]/D");
  saveTree->Branch("p2_pair", p2_pair, "p2_pair[nDimuons]/D");
  saveTree->Branch("p2x_pair", p2x_pair, "p2x_pair[nDimuons]/D");
  saveTree->Branch("p2y_pair", p2y_pair, "p2y_pair[nDimuons]/D");
  saveTree->Branch("p2z_pair", p2z_pair, "p2z_pair[nDimuons]/D");
  
  saveTree->Branch("p1_single", p1_single, "p1_single[nDimuons]/D");
  saveTree->Branch("p1x_single", p1x_single, "p1x_single[nDimuons]/D");
  saveTree->Branch("p1y_single", p1y_single, "p1y_single[nDimuons]/D");
  saveTree->Branch("p1z_single", p1z_single, "p1z_single[nDimuons]/D");
  saveTree->Branch("p2_single", p2_single, "p2_single[nDimuons]/D");
  saveTree->Branch("p2x_single", p2x_single, "p2x_single[nDimuons]/D");
  saveTree->Branch("p2y_single", p2y_single, "p2y_single[nDimuons]/D");
  saveTree->Branch("p2z_single", p2z_single, "p2z_single[nDimuons]/D");
 
  saveTree->Branch("r1_single", r1_single, "r1_single[nDimuons]/D");
  saveTree->Branch("r2_single", r2_single, "r2_single[nDimuons]/D");
  saveTree->Branch("r1_pair", r1_pair, "r1_pair[nDimuons]/D");
  saveTree->Branch("r2_pair", r2_pair, "r2_pair[nDimuons]/D");
  saveTree->Branch("angle_pair", angle_pair, "angle_pair[nDimuons]/D");

  saveTree->Branch("mass_pair", mass_pair, "mass_pair[nDimuons]/D");
  saveTree->Branch("pT_pair", pT_pair, "pT_pair[nDimuons]/D");
  saveTree->Branch("xF_pair", xF_pair, "xF_pair[nDimuons]/D");
  saveTree->Branch("x1_pair", x1_pair, "x1_pair[nDimuons]/D");
  saveTree->Branch("x2_pair", x2_pair, "x2_pair[nDimuons]/D");

  saveTree->Branch("mass_single", mass_single, "mass_single[nDimuons]/D");
  saveTree->Branch("pT_single", pT_single, "pT_single[nDimuons]/D");
  saveTree->Branch("xF_single", xF_single, "xF_single[nDimuons]/D");
  saveTree->Branch("x1_single", x1_single, "x1_single[nDimuons]/D");
  saveTree->Branch("x2_single", x2_single, "x2_single[nDimuons]/D");

  //Initialize track finder
  Log("Initializing the track finder and kalman filter ... ");
  VertexFit *vtxfit = new VertexFit();

  int offset = argc > 3 ? atoi(argv[3]) : 0;
  int nEvtMax = argc > 4 ? atoi(argv[4]) + offset : dataTree->GetEntries();
  if(nEvtMax > dataTree->GetEntries()) nEvtMax = dataTree->GetEntries();
  Log("Running from event " << offset << " through to event " << nEvtMax);
  for(int i = offset; i < nEvtMax; i++)
    {
      dataTree->GetEntry(i);
      Log("Processing event " << i << " with eventID = " << rawEvent->getEventID());
      
      int nTracklets = tracklets->GetEntries();
      vector<SRecTrack> tracks;
      vector<int> index_plus, index_minus;
      for(int j = 0; j < nTracklets; j++)
	{
	  Tracklet* tracklet = (Tracklet*)tracklets->At(j);
	  if(tracklet->getCharge() > 0)
	    {
	      index_plus.push_back(j);
	    }
	  else
	    {
	      index_minus.push_back(j);
	    }

	  SRecTrack track = tracklet->getSRecTrack();
	  track.setZVertex(vtxfit->findSingleMuonVertex(track));

	  tracks.push_back(track);
	}

      int nPlus = index_plus.size();
      int nMinus = index_minus.size();
      nDimuons = 0;
      for(int j = 0; j < nPlus; j++)
	{
	  SRecTrack _trackp = tracks[index_plus[j]];
	  for(int k = 0; k < nMinus; k++)
	    {
	      SRecTrack _trackm = tracks[index_minus[k]];

	      pID[nDimuons] = index_plus[j];
	      mID[nDimuons] = index_minus[k];

	      z_single[nDimuons] = (_trackp.getZVertex() + _trackm.getZVertex())/2.;
	      dz_single[nDimuons] = _trackp.getZVertex() - _trackm.getZVertex();
	      r1_single[nDimuons] = _trackp.getRVertex();
	      r2_single[nDimuons] = _trackm.getRVertex();
              chisq_single[nDimuons] = _trackp.getChisqVertex() + _trackm.getChisqVertex();

	      p1_single[nDimuons] = _trackp.getMomentumVertex(p1x_single[nDimuons], p1y_single[nDimuons], p1z_single[nDimuons]);
	      p2_single[nDimuons] = _trackm.getMomentumVertex(p2x_single[nDimuons], p2y_single[nDimuons], p2z_single[nDimuons]);
	      calc_variables(getMom(p1x_single[nDimuons], p1y_single[nDimuons], p1z_single[nDimuons]), getMom(p2x_single[nDimuons], p2y_single[nDimuons], p2z_single[nDimuons]), mass_single[nDimuons], pT_single[nDimuons], xF_single[nDimuons], x1_single[nDimuons], x2_single[nDimuons]);

	      //Vertex fitted dimuons
	      vtxfit->init();
	      vtxfit->addTrack(0, _trackp);
	      vtxfit->addTrack(1, _trackm);

	      vtxfit->addHypothesis(z_single[nDimuons], 50.);
	      hypoID[nDimuons] = vtxfit->processOneEvent();
	      z_pair[nDimuons] = vtxfit->getVertexZ0();
	      chisq_pair[nDimuons] = vtxfit->getKFChisq();
	      chisq_vx[nDimuons] = vtxfit->getVXChisq();

	      _trackp.setZVertex(z_pair[nDimuons]);
	      _trackm.setZVertex(z_pair[nDimuons]);
	      p1_pair[nDimuons] = _trackp.getMomentumVertex(p1x_pair[nDimuons], p1y_pair[nDimuons], p1z_pair[nDimuons]);
	      p2_pair[nDimuons] = _trackm.getMomentumVertex(p2x_pair[nDimuons], p2y_pair[nDimuons], p2z_pair[nDimuons]);
	      TLorentzVector p1_4mom = getMom(p1x_pair[nDimuons], p1y_pair[nDimuons], p1z_pair[nDimuons]);
	      TLorentzVector p2_4mom = getMom(p2x_pair[nDimuons], p2y_pair[nDimuons], p2z_pair[nDimuons]);
	      
	      r1_pair[nDimuons] = _trackp.getRVertex();
	      r2_pair[nDimuons] = _trackm.getRVertex();

	      vtxfit->print(); 
	      
	      mass_pair[nDimuons] = (p1_4mom + p2_4mom).M();
              angle_pair[nDimuons] = p1_4mom.Angle(p2_4mom.Vect());
	      calc_variables(p1_4mom, p2_4mom, mass_pair[nDimuons], pT_pair[nDimuons], xF_pair[nDimuons], x1_pair[nDimuons], x2_pair[nDimuons]);

	      nDimuons++;
	    }
	}

      if(nDimuons > 0) saveTree->Fill();
      tracklets->Clear();
    }

  saveFile->cd();
  saveTree->Write();
  saveFile->Close();

  delete vtxfit;

  return 1;
}
