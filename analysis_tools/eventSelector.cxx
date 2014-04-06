#include <iostream>
#include <stdlib.h>
#include <string>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>

#include "SRawEvent.h"
#include "SRecEvent.h"

using namespace std;

int main(int argc, char* argv[])
{
  double lo = atof(argv[3]);
  double hi = atof(argv[4]);

  TFile* dataFile = new TFile(argv[1], "READ");
  TTree* dataTree = (TTree*)dataFile->Get("save");

  SRawEvent* rawEvent = new SRawEvent();
  SRecEvent* recEvent = new SRecEvent();
  dataTree->SetBranchAddress("rawEvent", &rawEvent);
  dataTree->SetBranchAddress("recEvent", &recEvent);

  TFile* saveFile = new TFile(argv[2], "recreate");
  TTree* saveTree = new TTree("save", "save");
  
  double mass;
  saveTree->Branch("x", &mass, "mass/D");

  int nEntry = 0;
  for(int i = 0; i < dataTree->GetEntries(); ++i)
    {
      dataTree->GetEntry(i);
      
      double intensity = rawEvent->getIntensity();
      if(intensity < lo || intensity > hi)
	{
	  rawEvent->clear();
	  recEvent->clear();
  	  continue;
	}

      for(int j = 0; j < recEvent->getNDimuons(); ++j)
	{
	  SRecDimuon dimuon = recEvent->getDimuon(j);
	  double chisq_kf = dimuon.chisq_kf;
	  double p1x = dimuon.p_pos.Px();
	  double p2x = dimuon.p_neg.Px();
	  double xF = dimuon.xF;
	  double dz = dimuon.vtx_pos.Z() - dimuon.vtx_neg.Z();
	  mass = dimuon.mass;
	  if(chisq_kf < 5. && xF < 1 && p1x > p2x && abs(dz) < 150.)
	    {
    	      saveTree->Fill();
	    }
	}
      
      rawEvent->clear();
      recEvent->clear();
    }

  saveFile->cd();
  saveTree->Write();
  saveFile->Close();

  return 1;
}
