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
#include <TF1.h>
#include <TClonesArray.h>

#include "GeomSvc.h"
#include "SRawEvent.h"
#include "SRecEvent.h"
#include "FastTracklet.h"

using namespace std;

double findCenter(TH1D *hist)
{
  if(hist->GetEntries() < 1000) return 9999;

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

void linearFit(double x[], double y[], double w[], int n, double& a, double& b)
{
  if(n < 2)
    {
      std::cout << "Should have at least two points!!!" << std::endl;
      return; 
    }

  double sum, sx, sy, sxx, sxy, syy, det;
  sum = 0.;
  sx = 0.;
  sy = 0.;
  sxx = 0.;
  syy = 0.;
  sxy = 0.;

  for(int i = 0; i < n; i++)
    {
      sum += w[i];
      sx += w[i]*x[i];
      sy += w[i]*y[i];
      sxx += w[i]*x[i]*x[i];
      syy += w[i]*y[i]*y[i];
      sxy += w[i]*x[i]*y[i];
    }

  det = sum*sxx - sx*sx;
  if(fabs(det) < 1.0e-20)
    {
      a = 1.0e20;
      b = x[0];
      
      return;
    }

  a = (sum*sxy - sx*sy)/det;
  b = (sy*sxx - sxy*sx)/det;
}

int main(int argc, char *argv[])
{
  GeomSvc *p_geomSvc = GeomSvc::instance();
  p_geomSvc->init(GEOMETRY_VERSION);

  SRawEvent *rawEvent = new SRawEvent();
  TClonesArray* tracklets = new TClonesArray("Tracklet");

  TFile *dataFile = new TFile(argv[1], "READ");
  TTree *dataTree = (TTree *)dataFile->Get("save");

  dataTree->SetBranchAddress("rawEvent", &rawEvent);
  dataTree->SetBranchAddress("tracklets", &tracklets);

  int m_propIDs[] = {41, 43, 45, 47};
  vector<int> propIDs(m_propIDs, m_propIDs+sizeof(m_propIDs)/sizeof(int));
  const int nProps = propIDs.size();

  double res, z_exp, x_exp, y_exp, pos_exp;
  int propID, moduleID, elementID;

  TFile *saveFile = new TFile(argv[2], "recreate");
  TTree *saveTree = new TTree("save", "save");

  saveTree->Branch("z_exp", &z_exp, "z_exp/D");
  saveTree->Branch("x_exp", &x_exp, "x_exp/D");
  saveTree->Branch("y_exp", &y_exp, "y_exp/D");
  saveTree->Branch("pos_exp", &pos_exp, "pos_exp/D");
  saveTree->Branch("res", &res, "res/D");
  saveTree->Branch("propID", &propID, "propID/I");
  saveTree->Branch("moduleID", &moduleID, "moduleID/I");
  saveTree->Branch("elementID", &elementID, "elementID/I");

  //Intialization of arrays and hists  
  TH1D* hist[4][9];
  double offset[4][9];
  int nEffective[4];
  double a[4], b[4];
  for(int i = 0; i < nProps; i++)
    {
      string detectorName = p_geomSvc->getDetectorName(propIDs[i]);
      for(int j = 0; j < 9; j++)
	{
	  stringstream suffix;
	  suffix << j;
	  string histName = detectorName + "_" + suffix.str();

	  hist[i][j] = new TH1D(histName.c_str(), histName.c_str(), 200, -5.08, 5.08);
	  offset[i][j] = 9999.;
	}

      nEffective[i] = 0;
    }

  //Use tracks to fill residual distributions
  int nEvtMax = dataTree->GetEntries();
  for(int i = 0; i < nEvtMax; i++)
    {
      dataTree->GetEntry(i);
      if(tracklets->GetEntries() < 1) continue;

      rawEvent->reIndex("oa");
      vector<Hit> hitAll = rawEvent->getAllHits();

      //Only the first track is used, for simplicity
      Tracklet* _track = (Tracklet*)tracklets->At(0);
      for(int j = 0; j < nProps; j++)
	{
	  //If and only if one paired prop. tube hits are found, we proceed
	  list<SRawEvent::hit_pair> hit_pairs = rawEvent->getHitPairsInSuperDetector((propIDs[j]+1)/2);
	  if(hit_pairs.size() != 1) continue;

	  //Expected hit position at first plane
	  z_exp = p_geomSvc->getPlanePosition(propIDs[j]);
	  x_exp = _track->getExpPositionX(z_exp);
	  y_exp = _track->getExpPositionY(z_exp);
	  pos_exp = p_geomSvc->getUinStereoPlane(propIDs[j], x_exp, y_exp);
	  if(!p_geomSvc->isInPlane(propIDs[j], x_exp, y_exp)) continue;

	  propID = hitAll[hit_pairs.front().first].detectorID;
	  elementID = hitAll[hit_pairs.front().first].elementID;     
	  moduleID = Int_t((elementID-1)/8);                //Need to re-define moduleID for run2
	  res = pos_exp - hitAll[hit_pairs.front().first].pos;//p_geomSvc->getMeasurement(propID, elementID);
	  hist[j][moduleID]->Fill(res);
	  saveTree->Fill();

	  //Expected hit position at second plane
	  z_exp = p_geomSvc->getPlanePosition(propIDs[j]+1);
	  x_exp = _track->getExpPositionX(z_exp);
	  y_exp = _track->getExpPositionY(z_exp);
	  pos_exp = p_geomSvc->getUinStereoPlane(propIDs[j]+1, x_exp, y_exp);
	  if(!p_geomSvc->isInPlane(propIDs[j]+1, x_exp, y_exp)) continue;

	  propID = hitAll[hit_pairs.front().second].detectorID;
	  elementID = hitAll[hit_pairs.front().second].elementID;     
	  moduleID = Int_t((elementID-1)/8);                //Need to re-define moduleID for run2
	  res = pos_exp - hitAll[hit_pairs.front().second].pos;//p_geomSvc->getMeasurement(propID, elementID);
	  hist[j][moduleID]->Fill(res);
	  saveTree->Fill();
	}

      rawEvent->clear();
      tracklets->Clear();
    }

  //Process the residual hists, use linear fit to extrapolate
  for(int i = 0; i < nProps; i++)
    {
      double x[9], y[9], w[9];
      for(int j = 0; j < 9; j++)
	{
	  offset[i][j] = findCenter(hist[i][j]) + p_geomSvc->getPlaneWOffset(propIDs[i], j);
	  if(offset[i][j] < 100.)
	    {
	      x[nEffective[i]] = j;
	      y[nEffective[i]] = offset[i][j];
	      w[nEffective[i]] = hist[i][j]->GetEntries();

	      ++nEffective[i];
	    }
	}

      linearFit(x, y, w, nEffective[i], a[i], b[i]);
    }

  //Output alignment numbers to txt file
  ofstream fout(argv[3], ios::out);
  for(int i = 0; i < nProps; i++)
    {
      cout << " === " << p_geomSvc->getDetectorName(propIDs[i]) << " === " << endl;
      cout << " = " << nEffective[i] << "  " << a[i] << "  " << b[i] << endl;
      
      double x_real[9], x_ext[9];
      double y_real[9], y_ext[9];
      int n_real = 0;
      int n_ext = 0;
      for(int j = 0; j < 9; j++)
	{
	  if(offset[i][j] > 100.) 
	    {
	      offset[i][j] = a[i]*j + b[i];
              x_ext[n_ext] = j;
	      y_ext[n_ext] = offset[i][j];
	      n_ext++;
	      cout << "Extrapolated!  "; 
	    }
	  else
	    {
	      x_real[n_real] = j;
	      y_real[n_real] = offset[i][j];
	      n_real++;
	    }
	  cout << p_geomSvc->getDetectorName(propIDs[i]) << "  " << j << "  " << hist[i][j]->GetEntries() << "  " << offset[i][j] << endl;
	  fout << p_geomSvc->getPlaneWOffset(propIDs[i], j) + offset[i][j] << endl;
	}

      //Make plots
      TF1 line("line", "[0] + [1]*x", -1, 9);
      line.SetParameter(0, b[i]);
      line.SetParameter(1, a[i]);

      line.SetLineWidth(2);
      line.SetLineStyle(kDashed);

      TGraph gr_real(n_real, x_real, y_real);
      TGraph gr_ext(n_ext, x_ext, y_ext);

      gr_real.SetMarkerColor(kRed);
      gr_real.SetMarkerStyle(8);
      gr_ext.SetMarkerColor(kBlue);
      gr_ext.SetMarkerStyle(4);

      TCanvas c;
      c.cd();
      line.Draw();
      line.GetYaxis()->SetRangeUser(-3., 3.);
      c.Update();
      gr_ext.Draw("P same");
      gr_real.Draw("P same");
      c.SaveAs(p_geomSvc->getDetectorName(propIDs[i]).c_str());
    }

  //Save temporary results to ROOT file
  saveFile->cd();
  saveTree->Write();
  for(int i = 0; i < nProps; i++)
    {
      for(int j = 0; j < 9; j++)
	{
	  hist[i][j]->Write();
	}
    }
  saveFile->Close();

  return 1;
}
