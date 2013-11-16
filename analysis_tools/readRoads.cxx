#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <string>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>

#include "GeomSvc.h"
#include "TriggerRoad.h"
#include "TriggerAnalyzer.h"

using namespace std;

int main(int argc, char* argv[])
{
  ////control variables
  string dy_input = argv[1];
  string gun_input = argv[2];
  string road_output = argv[3];
  double cut_td = atof(argv[4]);
  double cut_gun = atof(argv[5]);
  int nPtBin = atoi(argv[6]);
  double pt_min = atof(argv[7]);
  double pt_max = atof(argv[8]);

  GeomSvc* p_geomSvc = GeomSvc::instance();
  p_geomSvc->init(GEOMETRY_VERSION);

  ////Raw MC hodo hits info.
  Double_t mass;
  Double_t weight;
  Double_t z0;

  Double_t px1, py1;
  Double_t px2, py2;

  Int_t nHodoHits1;
  Int_t detectorID1[100];
  Int_t elementID1[100];

  Int_t nHodoHits2;
  Int_t detectorID2[100];
  Int_t elementID2[100];

  TFile* dataFile = new TFile(dy_input.c_str(), "READ");
  TTree* dataTree = (TTree*)dataFile->Get("save");

  dataTree->SetBranchAddress("weight", &weight);
  dataTree->SetBranchAddress("z0", &z0);
  dataTree->SetBranchAddress("mass", &mass);

  dataTree->SetBranchAddress("nHodoHits1", &nHodoHits1);
  dataTree->SetBranchAddress("detectorID1", detectorID1);
  dataTree->SetBranchAddress("elementID1", elementID1);

  dataTree->SetBranchAddress("px1", &px1);
  dataTree->SetBranchAddress("py1", &py1);

  dataTree->SetBranchAddress("nHodoHits2", &nHodoHits2);
  dataTree->SetBranchAddress("detectorID2", detectorID2);
  dataTree->SetBranchAddress("elementID2", elementID2);

  dataTree->SetBranchAddress("px2", &px2);
  dataTree->SetBranchAddress("py2", &py2);

  TFile* saveFile = new TFile(road_output.c_str(), "recreate");
  TTree* saveTree1 = dataTree->CloneTree(0);

  ////Read the hit combination and make TriggerRoad list for every possible entry
  list<TriggerRoad> p_roads, m_roads;
  for(int i = 0; i < dataTree->GetEntries(); i++)
    {
      dataTree->GetEntry(i);
      
      list<TriggerRoad> p_road_add = TriggerRoad::makeRoadList(nHodoHits1, detectorID1, elementID1, z0, mass, sqrt(px1*px1 + py1*py1), weight);
      list<TriggerRoad> m_road_add = TriggerRoad::makeRoadList(nHodoHits2, detectorID2, elementID2, z0, mass, sqrt(px2*px2 + py2*py2), weight);
     
      p_roads.insert(p_roads.end(), p_road_add.begin(), p_road_add.end());
      m_roads.insert(m_roads.end(), m_road_add.begin(), m_road_add.end());
      if(p_road_add.empty() || m_road_add.empty()) continue;  

      saveTree1->Fill();
    }

  ////Make unique road lists 
  //Make unique positive roads 
  list<TriggerRoad> p_roads_unique; p_roads_unique.clear();
  for(list<TriggerRoad>::iterator iter = p_roads.begin(); iter != p_roads.end(); ++iter)
    {
      list<TriggerRoad>::iterator road = find(p_roads_unique.begin(), p_roads_unique.end(), *iter);
      if(road == p_roads_unique.end())
	{
	  p_roads_unique.push_back(*iter);
	}
      else
	{
	  (*road) += (*iter);
	}
    }

  //Make unique negative roads
  list<TriggerRoad> m_roads_unique; m_roads_unique.clear();
  for(list<TriggerRoad>::iterator iter = m_roads.begin(); iter != m_roads.end(); ++iter)
    {
      list<TriggerRoad>::iterator road = find(m_roads_unique.begin(), m_roads_unique.end(), *iter);
      if(road == m_roads_unique.end())
	{
	  m_roads_unique.push_back(*iter);
	}
      else
	{
	  (*road) += (*iter);
	}
    }

  //Remove duplicate roads that are present in both plus and minus roads
  for(list<TriggerRoad>::iterator p_road = p_roads_unique.begin(); p_road != p_roads_unique.end(); )
    {
      list<TriggerRoad>::iterator m_road = find(m_roads_unique.begin(), m_roads_unique.end(), *p_road);
      if(m_road != m_roads_unique.end())
	{
	  if(p_road->weight() < m_road->weight())
	    {
	      p_road = p_roads_unique.erase(p_road);
	    }
	  else
	    {
	      m_roads_unique.erase(m_road);
	      ++p_road;
	    }
	}
      else
	{
	  ++p_road;
	}
    }

  ////Assign groupID and unique road ID, drop low-pT roads
  int uniqueID = 0;
  double binWidth = (pt_max - pt_min)/nPtBin;
  for(list<TriggerRoad>::iterator iter = p_roads_unique.begin(); iter != p_roads_unique.end(); ) 
    {
      if(iter->pT_mean < pt_min || iter->pT_mean > pt_max)
	{
	  iter = p_roads_unique.erase(iter);
	  continue;
	}
	  
      iter->roadID = uniqueID;
      iter->groupID = (int((iter->pT_mean - pt_min)/binWidth) + 1)*iter->getTB();
      
      ++iter;
      ++uniqueID;
    }

  uniqueID = 0;
  for(list<TriggerRoad>::iterator iter = m_roads_unique.begin(); iter != m_roads_unique.end(); ) 
    {
      if(iter->pT_mean < pt_min || iter->pT_mean > pt_max)
	{
	  iter = m_roads_unique.erase(iter);
	  continue;
	}
	  
      iter->roadID = uniqueID;
      iter->groupID = (int((iter->pT_mean - pt_min)/binWidth) + 1)*iter->getTB();
      
      ++iter;
      ++uniqueID;
    }

  ////Calculate gun frequency
  int nHits;
  int detectorID[500];
  int elementID[500];

  TFile* gunFile = new TFile(gun_input.c_str(), "READ");
  TTree* gunTree = (TTree*)gunFile->Get("save");

  gunTree->SetBranchAddress("nHits", &nHits);
  gunTree->SetBranchAddress("detectorID", detectorID);
  gunTree->SetBranchAddress("elementID", elementID); 

  TriggerAnalyzer* triggerAna = new TriggerAnalyzer();
  triggerAna->init(p_roads_unique, m_roads_unique);
  triggerAna->buildTriggerTree();

  double factor = 50E6/gunTree->GetEntries();
  for(int i = 0; i < gunTree->GetEntries(); ++i)
    {
      gunTree->GetEntry(i);
      triggerAna->acceptEvent(nHits, detectorID, elementID);

      list<TriggerRoad> p_roads_found = triggerAna->getRoadsFound(+1);
      for(list<TriggerRoad>::iterator iter = p_roads_found.begin(); iter != p_roads_found.end(); ++iter)
	{
	  list<TriggerRoad>::iterator p_road_found = find(p_roads_unique.begin(), p_roads_unique.end(), *iter);
	  if(p_road_found != p_roads_unique.end()) p_road_found->rndf += 1.;
	}

      list<TriggerRoad> m_roads_found = triggerAna->getRoadsFound(-1);
      for(list<TriggerRoad>::iterator iter = m_roads_found.begin(); iter != m_roads_found.end(); ++iter)
	{
	  list<TriggerRoad>::iterator m_road_found = find(m_roads_unique.begin(), m_roads_unique.end(), *iter);
	  if(m_road_found != m_roads_unique.end()) m_road_found->rndf += 1.;
	}
    }

  ////Output to ROOT file that stores the roads
  TriggerRoad* p_road = new TriggerRoad(); p_road->clear();
  TriggerRoad* m_road = new TriggerRoad(); m_road->clear();

  saveFile->cd();
  TTree* saveTree2 = new TTree("single_p", "single_p");
  saveTree2->Branch("road", &p_road, 256000, 99);
  for(list<TriggerRoad>::iterator iter = p_roads_unique.begin(); iter != p_roads_unique.end(); ++iter)
    {
      *p_road = *iter;
      p_road->rndf = p_road->rndf*factor;
      saveTree2->Fill();
    }

  TTree* saveTree3 = new TTree("single_m", "single_m");
  saveTree3->Branch("road", &m_road, 256000, 99);
  for(list<TriggerRoad>::iterator iter = m_roads_unique.begin(); iter != m_roads_unique.end(); ++iter)
    {
      *m_road = *iter;
      m_road->rndf = m_road->rndf*factor;
      saveTree3->Fill();
    }

  saveFile->cd();
  saveTree1->Write();
  saveTree2->Write();
  saveTree3->Write();
  saveFile->Close();

  //Output to txt file
  triggerAna->init(road_output.c_str(), cut_td, cut_gun);
  triggerAna->outputEnabled();

  delete triggerAna;
  return 1;
}
