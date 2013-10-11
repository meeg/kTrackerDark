#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <stdio.h>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TRandom.h>
#include <TMatrixD.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include <TH1I.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TClonesArray.h>

#include "GeomSvc.h"
#include "SRawEvent.h"
#include "SRecEvent.h"
#include "FastTracklet.h"

using namespace std;

double findCenter(TH1D *hist, double spacing)
{
  if(hist->GetEntries() < 150) return 9999;
  if(hist->GetRMS() < spacing/5. && hist->GetEntries() < 1000) return 9999;

  int nBin = hist->GetNbinsX();
  int nBinInSize = nBin/2;
  double binWidth = hist->GetBinWidth(1);

  //cout << hist->GetName() << endl;
  //cout << nBin << "  " << nBinInSize << "  " << binWidth << endl;

  int nEvt_max = 0;
  int index_max_left = 0;
  for(int i = 1; i <= nBinInSize; ++i)
    {
      int nEvt_curr = hist->Integral(i, i + nBinInSize);
      //cout << i << " : " << hist->GetBinCenter(i) << " <===> " << hist->GetBinCenter(i + nBinInSize) << " : " << nEvt_curr << " === " << nEvt_max << endl;
    
      if(nEvt_curr > nEvt_max)
	{
	  nEvt_max = nEvt_curr;
	  index_max_left = i;
	}
    }

  nEvt_max = 0;
  int index_max_right = nBin;
  for(int i = nBin; i >= nBinInSize; --i)
    {
      int nEvt_curr = hist->Integral(i - nBinInSize, i);
      
      if(nEvt_curr > nEvt_max)
	{
	  nEvt_max = nEvt_curr;
	  index_max_right = i;
	}
    }

  return (hist->GetBinCenter(index_max_left) + hist->GetBinCenter(index_max_right))/2.;
}

int main(int argc, char *argv[])
{
  //remove("alignment_hodo.txt");

  GeomSvc *p_geomSvc = GeomSvc::instance();
  p_geomSvc->init(GEOMETRY_VERSION);

  SRawEvent *rawEvent = new SRawEvent();
  TClonesArray* tracklets = new TClonesArray("Tracklet");

  TFile *dataFile = new TFile(argv[1], "READ");
  TTree *dataTree = (TTree *)dataFile->Get("save");

  dataTree->SetBranchAddress("rawEvent", &rawEvent);
  dataTree->SetBranchAddress("tracklets", &tracklets);

  int m_hodoIDs[] = {25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40};
  vector<int> hodoIDs(m_hodoIDs, m_hodoIDs+sizeof(m_hodoIDs)/sizeof(int));
  const int nHodos = hodoIDs.size();

  double z_exp, x_exp, y_exp, pos_exp;
  int hodoID, elementID;
  int nFired;

  TFile *saveFile = new TFile(argv[2], "recreate");
  TTree *saveTree = new TTree("save", "save");

  saveTree->Branch("z_exp", &z_exp, "z_exp/D");
  saveTree->Branch("x_exp", &x_exp, "x_exp/D");
  saveTree->Branch("y_exp", &y_exp, "y_exp/D");
  saveTree->Branch("pos_exp", &pos_exp, "pos_exp/D");
  saveTree->Branch("hodoID", &hodoID, "hodoID/I");
  saveTree->Branch("elementID", &elementID, "elementID/I");
  saveTree->Branch("nFired", &nFired, "nFired/I");

  //Initialization of container and hists
  vector<TH1D*> hist[nHodos];
  vector<double> offset[nHodos];
  int nElement[nHodos];
  double spacing[nHodos];
  for(int i = 0; i < nHodos; i++)
    {
      nElement[i] = p_geomSvc->getPlaneNElements(hodoIDs[i]);
      spacing[i] = p_geomSvc->getPlaneSpacing(hodoIDs[i]);
      string detectorName = p_geomSvc->getDetectorName(hodoIDs[i]);
      for(int j = 0; j < nElement[i]; j++)
	{
	  stringstream suffix;
	  suffix << j;
	  string histName = detectorName + "_" + suffix.str();

	  TH1D *hist_temp = new TH1D(histName.c_str(), histName.c_str(), 200, -spacing[i], spacing[i]);
	  hist[i].push_back(hist_temp);
	  offset[i].push_back(9999.);
	}
    }

  //User tracks to fill residual distributions
  int nEvtMax = dataTree->GetEntries();
  for(int i = 0; i < nEvtMax; i++)
    {
      dataTree->GetEntry(i);
      if(tracklets->GetEntries() < 1) continue;
 
      rawEvent->reIndex("oa");
      vector<Hit> hitAll = rawEvent->getAllHits();

      //Only the first track is used, for simplicity
      Tracklet* _track = (Tracklet*)tracklets->At(0);
      for(int j = 0; j < nHodos; j++)
	{
	  if((rawEvent->getNHitsInDetector(hodoIDs[j]) != 1) && (rawEvent->getNHitsInDetector(hodoIDs[j]) != 2)) continue;

	  //Expected hit position on one hodo plane
	  z_exp = p_geomSvc->getPlanePosition(hodoIDs[j]);
	  x_exp = _track->getExpPositionX(z_exp);
	  y_exp = _track->getExpPositionY(z_exp);
	  pos_exp = p_geomSvc->getUinStereoPlane(hodoIDs[j], x_exp, y_exp);

	  //if the hit is outside of that hodo, skip
	  if(!p_geomSvc->isInPlane(hodoIDs[j], x_exp, y_exp)) continue;

	  //Get the hits on that plane, in the case of 2 hits if they are not next to each other, skip 
	  list<int> hitlist = rawEvent->getHitsIndexInDetector(hodoIDs[j]);
	  nFired = hitlist.size();
	  if(nFired == 2 && abs(hitAll[hitlist.front()].elementID - hitAll[hitlist.back()].elementID) > 1) continue;

	  //Fill the first or the only hit
	  hodoID = hitAll[hitlist.front()].detectorID;
	  elementID = hitAll[hitlist.front()].elementID;     
	  hist[j][elementID-1]->Fill(pos_exp - p_geomSvc->getMeasurement(hodoIDs[j], elementID)); 
	  saveTree->Fill();

	  //Fill the second hit if available
	  if(nFired != 2) continue;
	  hodoID = hitAll[hitlist.back()].detectorID;
	  elementID = hitAll[hitlist.back()].elementID;
	  hist[j][elementID-1]->Fill(pos_exp - p_geomSvc->getMeasurement(hodoIDs[j], elementID)); 
	
	  saveTree->Fill();
	}

      rawEvent->clear();
      tracklets->Clear();
    }

  //Process the residual distributions
  double offset_plane[nHodos];
  int nValidEntries[nHodos];
  for(int i = 0; i < nHodos; i++)
    {
      offset_plane[i] = 0.;
      nValidEntries[i] = 0;
      for(int j = 0; j < nElement[i]; j++)
	{
	  int nEntries = hist[i][j]->GetEntries();
	  if(nEntries < 500)
	    {
	      hist[i][j]->Rebin(2);
	    }
	 
	  offset[i][j] = findCenter(hist[i][j], spacing[i]);
	  if(offset[i][j] < 100.)
	    {
	      offset_plane[i] += offset[i][j]*hist[i][j]->GetEntries();
	      nValidEntries[i] += hist[i][j]->GetEntries();
	    }
	}

      offset_plane[i] = offset_plane[i]/nValidEntries[i];
    }

  //Output alignment parameters into the txt file
  ofstream fout(argv[3], ios::out);
  for(int i = 0; i < nHodos; i++)
    {
      cout << " === " << p_geomSvc->getDetectorName(hodoIDs[i]) << "  " << offset_plane[i] << endl;
      fout << offset_plane[i] + p_geomSvc->getPlaneWOffset(hodoIDs[i]) << endl;
      for(int j = 0; j < nElement[i]; j++)
	{
	  if(offset[i][j] > 100.) offset[i][j] = offset_plane[i];
	  cout << p_geomSvc->getDetectorName(hodoIDs[i]) << "  " << j << "  " << hist[i][j]->GetEntries() << "  " << offset[i][j] << endl;
	}
    }

  //Save the temporary results into the ROOT file
  saveFile->cd();
  saveTree->Write();
  for(int i = 0; i < nHodos; i++)
    {
      for(int j = 0; j < nElement[i]; j++)
	{
	  hist[i][j]->Write();
          //cout << p_geomSvc->getDetectorName(hodoIDs[i]) << "  " << j << "  " << hist[i][j]->GetEntries() << "  " << findCenter(hist[i][j]) << endl;
	}
    }
  saveFile->Close();
  
  return 1;
}