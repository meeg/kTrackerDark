#include <fstream>
#include <sstream>

#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>

#include "TriggerAnalyzer.h"
#include "JobOptsSvc.h"

#define REQUIRE_TB

TNode::TNode(int uID)
{
  uniqueID = uID;
  children.clear();
}

void TNode::add(TNode* child)
{
  children.push_back(child);
}

TriggerAnalyzer::TriggerAnalyzer()
{
  GeomSvc* p_geomSvc = GeomSvc::instance();

  detectorIDs_trigger = p_geomSvc->getDetectorIDs("H1[TB]");
  std::vector<int> H2X_trigger = p_geomSvc->getDetectorIDs("H2[TB]");
  std::vector<int> H3X_trigger = p_geomSvc->getDetectorIDs("H3[TB]");
  std::vector<int> H4X_trigger = p_geomSvc->getDetectorIDs("H4[TB]");

  detectorIDs_trigger.insert(detectorIDs_trigger.end(), H2X_trigger.begin(), H2X_trigger.end());
  detectorIDs_trigger.insert(detectorIDs_trigger.end(), H3X_trigger.begin(), H3X_trigger.end());
  detectorIDs_trigger.insert(detectorIDs_trigger.end(), H4X_trigger.begin(), H4X_trigger.end());

  root[0] = NULL;
  root[1] = NULL;
}

TriggerAnalyzer::~TriggerAnalyzer()
{
  clear(1); clear(-1);
}

bool TriggerAnalyzer::init(std::string schemaName)
{
  JobOptsSvc* p_jobOptsSvc = JobOptsSvc::instance();
  GeomSvc* p_geomSvc = GeomSvc::instance();

  char query[300];
  sprintf(query, "SELECT charge,St1DetectorName,St1ElementID,St2DetectorName,St2ElementID,"
	  "St3DetectorName,St3ElementID,St4DetectorName,St4ElementID FROM %s.TriggerRoads", schemaName.c_str());
  
  TSQLServer* server = TSQLServer::Connect(p_jobOptsSvc->m_mySQLurl.c_str(), "seaguest","qqbar2mu+mu-");
  if(server == NULL) return false;

  TSQLResult* res = server->Query(query);
  if(res == NULL) return false;

  unsigned int nRoads = res->GetRowCount();
  for(unsigned int i = 0; i < nRoads; ++i)
    {
      TSQLRow* row = res->Next();

      TriggerRoad road_new;
      for(int j = 1; j <= 4; ++j)
	{
	  road_new.addElement(p_geomSvc->getDetectorID(row->GetField(2*j-1)), atoi(row->GetField(2*j)));
	}
      road_new.enable();

      int charge = atoi(row->GetField(0));
      if(charge > 0)
	{
	  roads_enabled[0].push_back(road_new);
	}
      else
	{
	  roads_enabled[1].push_back(road_new);
	}

      delete row;
    }

  delete res;
  delete server;

  return true;
}

bool TriggerAnalyzer::init()
{
  using namespace std;

  JobOptsSvc* p_jobOptsSvc = JobOptsSvc::instance();

  std::string fileNames[4] = {p_jobOptsSvc->m_roadFile_pt, p_jobOptsSvc->m_roadFile_pb, p_jobOptsSvc->m_roadFile_mt, p_jobOptsSvc->m_roadFile_mb};
  char buffer[300];
  int pRoads = 0;
  int mRoads = 0;
  for(int i = 0; i < 4; ++i)
    {
      fstream fin(fileNames[i].c_str(), ios::in);

      while(fin.getline(buffer, 300))
	{
	  istringstream stringBuf(buffer);

	  int elementIDs[4];
	  int roadID;
	  int groupID;
	  int charge;
	  stringBuf >> roadID >> elementIDs[0] >> elementIDs[1] >> elementIDs[2] >> elementIDs[3] >> charge >> groupID;

	  TriggerRoad road_new;
	  road_new.groupID = groupID;
	  if(i == 0 || i == 2)
	    {
	      road_new.addElement(26, elementIDs[0]);
	      road_new.addElement(32, elementIDs[1]);
	      road_new.addElement(34, elementIDs[2]);
	      road_new.addElement(40, elementIDs[3]);
	    }
	  else
	    {
	      road_new.addElement(25, elementIDs[0]);
	      road_new.addElement(31, elementIDs[1]);
	      road_new.addElement(33, elementIDs[2]);
	      road_new.addElement(39, elementIDs[3]);
	    
	      road_new.groupID *= -1;
	    }

	  if(i < 2)
	    {
	      road_new.roadID = pRoads++;
	      roads_enabled[0].push_back(road_new);
	    }
	  else
	    {
	      road_new.roadID = mRoads++;
	      roads_enabled[1].push_back(road_new);
	    }
	}
    }
  makeRoadPairs();

  std::cout << "TriggerAnalyzer: " << roads_enabled[0].size() << " positive roads and " << roads_enabled[1].size() << " negative roads are activated." << std::endl;
  return roads_enabled[0].size() > 0 && roads_enabled[1].size();
}

bool TriggerAnalyzer::init(std::string fileName, double cut_td, double cut_gun)
{
  TriggerRoad* road = new TriggerRoad(); road->clear();

  TFile dataFile(fileName.c_str(), "READ");
  TTree* p_dataTree = (TTree*)dataFile.Get("single_p");
  TTree* m_dataTree = (TTree*)dataFile.Get("single_m");

  p_dataTree->SetBranchAddress("road", &road);
  m_dataTree->SetBranchAddress("road", &road);

  roads[0].clear();
  for(int i = 0; i < p_dataTree->GetEntries(); i++)
    {
      p_dataTree->GetEntry(i);
      if(road->isValid()) roads[0].push_back(*road);
      road->clear();
    }

  roads[1].clear();
  for(int i = 0; i < m_dataTree->GetEntries(); i++)
    {
      m_dataTree->GetEntry(i);
      if(road->isValid()) roads[1].push_back(*road);
      road->clear();
    }

  filterRoads(cut_td, cut_gun);
  makeRoadPairs();
  return true;
}

bool TriggerAnalyzer::init(std::list<TriggerRoad> p_roads, std::list<TriggerRoad> m_roads, double cut_td, double cut_gun)
{
  roads[0].clear();
  roads[0].assign(p_roads.begin(), p_roads.end());

  roads[1].clear();
  roads[1].assign(m_roads.begin(), m_roads.end());

  filterRoads(cut_td, cut_gun);
  makeRoadPairs();
  return true;
}

void TriggerAnalyzer::filterRoads(double cut_td, double cut_gun)
{
  //Filter road list
  roads_enabled[0].clear();
  roads_disabled[0].clear();
  for(std::list<TriggerRoad>::iterator road = roads[0].begin(); road != roads[0].end(); ++road)
    {
      if(road->ratio() < cut_td || road->rndf > cut_gun)
	{
	  road->disable();
	  roads_disabled[0].push_back(*road);
	}
      else
	{
	  road->enable();
	  roads_enabled[0].push_back(*road);
	}
    }

  roads_enabled[1].clear();
  roads_disabled[1].clear();
  for(std::list<TriggerRoad>::iterator road = roads[1].begin(); road != roads[1].end(); ++road)
    {
      if(road->ratio() < cut_td || road->rndf > cut_gun)
	{
	  road->disable();
	  roads_disabled[1].push_back(*road);
	}
      else
	{
	  road->enable();
	  roads_enabled[1].push_back(*road);
	}
    }

  std::cout << "Loaded " << roads[0].size() << " positive roads and " << roads[1].size() << " negative roads" << std::endl;
  std::cout << roads_enabled[0].size() << " positive roads and " << roads_enabled[1].size() << " negative roads are activated." << std::endl;
}

void TriggerAnalyzer::makeRoadPairs()
{
  //Form accepted trigger pair and rejected trigger pair
  triggers.clear();
  for(int i = 1; i <= 7; i++)
    {
      for(int j = 1; j <= 7; j++)
	{
	  //if(abs(i-j) < 4 && i+j > 10) 
	    {
#ifndef REQUIRE_TB
	      //Don't separate top/bottom
	      triggers.push_back(std::make_pair(i, j));
#else
	      //Separate top/bottom, i.e. one must come from top and the other should be from bottom
	      triggers.push_back(std::make_pair(-i, j));
	      triggers.push_back(std::make_pair(i, -j));
#endif
	    } 
	}
    }
}

bool TriggerAnalyzer::acceptEvent(TriggerRoad& p_road, TriggerRoad& m_road)
{
  std::list<TriggerRoad>::iterator p_iter = std::find(roads[0].begin(), roads[0].end(), p_road);
  std::list<TriggerRoad>::iterator m_iter = std::find(roads[1].begin(), roads[1].end(), m_road);

  if(p_iter != roads_enabled[0].end() && m_iter != roads_enabled[1].end())
    {
      return true;
    }
  return false;
}

bool TriggerAnalyzer::buildData(int nHits, int detectorIDs[], int elementIDs[])
{
  GeomSvc* p_geomSvc = GeomSvc::instance();
  data.clear();

  //Form data matrix
  std::set<int> vertex;
  std::set<int> HX[4];

  vertex.insert(-1);
  for(int i = 0; i < nHits; i++)
    {
      std::string detectorName = p_geomSvc->getDetectorName(detectorIDs[i]);
      int uniqueID = detectorIDs[i]*100 + elementIDs[i];
      if(detectorName.find("H1T") != std::string::npos || detectorName.find("H1B") != std::string::npos)
	{
	  HX[0].insert(uniqueID);
	}
      else if(detectorName.find("H2T") != std::string::npos || detectorName.find("H2B") != std::string::npos)
	{
	  HX[1].insert(uniqueID);
	}
      else if(detectorName.find("H3T") != std::string::npos || detectorName.find("H3B") != std::string::npos)
	{
	  HX[2].insert(uniqueID);
	}
      else if(detectorName.find("H4T") != std::string::npos || detectorName.find("H4B") != std::string::npos)
	{
	  HX[3].insert(uniqueID);
	}
    }

  for(int i = 0; i < 4; i++)
    {
      if(HX[i].empty()) return false;
    }

  data.push_back(vertex);
  data.push_back(HX[0]);
  data.push_back(HX[1]);
  data.push_back(HX[2]);
  data.push_back(HX[3]); 

  return true;
}

bool TriggerAnalyzer::acceptEvent(SRawEvent* rawEvent, int mode)
{
  int nHits = 0;
  int detectorIDs[10000];
  int elementIDs[10000];

  std::vector<Hit> triggerHits;
  if(mode == USE_TRIGGER_HIT)
    {
      triggerHits = rawEvent->getTriggerHits();
    }
  else if(mode == USE_HIT)
    {
      triggerHits = rawEvent->getAllHits();
    }

  for(std::vector<Hit>::iterator iter = triggerHits.begin(); iter != triggerHits.end(); ++iter)
    {
      if(iter->detectorID < 25 || iter->detectorID > 40) continue;
      if(iter->inTime != 1) continue;

      detectorIDs[nHits] = iter->detectorID;
      elementIDs[nHits] = iter->elementID;
    
      ++nHits;
    }

  return acceptEvent(nHits, detectorIDs, elementIDs);
}

bool TriggerAnalyzer::acceptEvent(int nHits, int detectorIDs[], int elementIDs[])
{
  //initialize the container and counter
  roads_found[0].clear();
  roads_found[1].clear();
  nRoads[0][0] = 0; nRoads[0][1] = 0;
  nRoads[1][0] = 0; nRoads[1][1] = 0;
  
  //Build data structure
  if(!buildData(nHits, detectorIDs, elementIDs)) return false;

  //search for patterns
  for(int i = 0; i < 2; i++)
    {
      roads_temp.clear();
      search(root[i], data, 0, i);
    }

  //road selection criteria: remove invalid roads
  std::set<int> groupIDs[2];
  for(int i = 0; i < 2; i++)
    {
      for(std::list<TriggerRoad>::iterator iter = roads_found[i].begin(); iter != roads_found[i].end(); )
	{
    	  std::list<TriggerRoad>::iterator road = std::find(roads_enabled[i].begin(), roads_enabled[i].end(), *iter);
          if(road == roads_enabled[i].end())
	    {
	      iter = roads_found[i].erase(iter);
	      continue;
	    }
	  else
	    {
	      *iter = *road;
              iter->groupID > 0 ? ++nRoads[i][0] : ++nRoads[i][1];	

#ifndef REQUIRE_TB 
	      //Don't separate top/bottom
	      groupIDs[i].insert(abs(road->groupID));
#else
	      //Separate top/bottom
	      groupIDs[i].insert(road->groupID);
#endif
	      ++iter;
	    }
	}
    }
  if(groupIDs[0].empty() || groupIDs[1].empty()) return false;
  
  //Lv-2 selection based on group id selection
  for(std::set<int>::iterator p_groupID = groupIDs[0].begin(); p_groupID != groupIDs[0].end(); ++p_groupID)
    {
      for(std::set<int>::iterator m_groupID = groupIDs[1].begin(); m_groupID != groupIDs[1].end(); ++m_groupID)
	{
	  std::list<Trigger>::iterator trigger = std::find(triggers.begin(), triggers.end(), std::make_pair(*p_groupID, *m_groupID));
	  if(trigger != triggers.end()) return true;
	}
    }

  return false;
}

void TriggerAnalyzer::search(TNode* root, DataMatrix& data, int level, int charge)
{
  roads_temp.push_back(root->uniqueID);
  if(root->children.empty())
    {
      //printRoadFound();
      roads_found[charge].push_back(TriggerRoad(roads_temp));
      roads_temp.pop_back();

      return;
    }

  for(std::list<TNode*>::iterator iter = root->children.begin(); iter != root->children.end(); ++iter)
    {
      if(data[level+1].find((*iter)->uniqueID) == data[level+1].end()) continue;
      search(*iter, data, level+1, charge);
    }
  roads_temp.pop_back();
}

void TriggerAnalyzer::printTree(TNode* root)
{
  if(root == NULL) return;

  roads_temp.push_back(root->uniqueID);
  if(root->children.empty())
    {
      printRoadFound();
      roads_temp.pop_back();

      return;
    }

  for(std::list<TNode*>::iterator iter = root->children.begin(); iter != root->children.end(); ++iter)
    {
      clearTree(*iter);
    }
  roads_temp.pop_back();
}

void TriggerAnalyzer::clearTree(TNode* root)
{
  if(root == NULL) return;

  if(root->children.empty())
    {
      delete root;
      return;
    }

  for(std::list<TNode*>::iterator iter = root->children.begin(); iter != root->children.end(); ++iter)
    {
      clearTree(*iter);
    }
  delete root;
}

void TriggerAnalyzer::buildTriggerTree()
{
  for(int i = 0; i < 2; i++)
    {
      root[i] = new TNode(-1);
      for(std::list<TriggerRoad>::iterator iter = roads_enabled[i].begin(); iter != roads_enabled[i].end(); ++iter)
	{
	  if(!iter->isEnabled()) continue;

	  TNode* parentNode[5]; //note: index 4 is useless, just to keep the following code simpler
	  parentNode[0] = root[i];
	  for(int j = 0; j < 4; j++)
	    {
	      int uniqueID = iter->getDetectorID(j)*100 + iter->getElementID(j);
	      bool isNewNode = true;
	      for(std::list<TNode*>::iterator iter = parentNode[j]->children.begin(); iter != parentNode[j]->children.end(); ++iter)
		{
		  if(uniqueID == (*iter)->uniqueID) 
		    {
		      parentNode[j+1] = *iter;
		      isNewNode = false;

		      break;
		    }
		}

	      if(isNewNode)
		{
		  TNode* node_new = new TNode(uniqueID);
		  parentNode[j]->add(node_new);

		  parentNode[j+1] = node_new;
		}
	    }
	}
    }
}

void TriggerAnalyzer::printRoadFound()
{
  std::cout << "Found one road: " << std::endl;
  for(std::list<int>::iterator iter = roads_temp.begin(); iter != roads_temp.end(); ++iter)
    {
      std::cout << *iter << " --> ";
    }
  std::cout << std::endl; 
}

void TriggerAnalyzer::printData(DataMatrix& data)
{
  for(unsigned int i = 0; i < data.size(); i++)
    {
      std::cout << i << ": ";
      for(std::set<int>::iterator iter = data[i].begin(); iter != data[i].end(); ++iter)
	{
	  std::cout << *iter << "  ";
	}
      std::cout << std::endl;
    }
}

void TriggerAnalyzer::outputEnabled()
{
  using namespace std;

  //Output single roads enabled
  fstream fout_road[2][2];
  fout_road[0][0].open("roads_plus_top.txt", ios::out);
  fout_road[0][1].open("roads_plus_bottom.txt", ios::out);
  fout_road[1][0].open("roads_minus_top.txt", ios::out);
  fout_road[1][1].open("roads_minus_bottom.txt", ios::out);

  roads_enabled[0].sort(TriggerRoad::byPt);
  roads_enabled[1].sort(TriggerRoad::byPt);
  for(int i = 0; i < 2; i++)
    {
      for(std::list<TriggerRoad>::iterator iter = roads_enabled[i].begin(); iter != roads_enabled[i].end(); ++iter)
	{
	  if(!iter->isEnabled()) continue;

	  int tb = iter->getTB() > 0 ? 0 : 1;
	  fout_road[i][tb] << *iter << endl;
	}
    }

  fout_road[0][0].close();
  fout_road[0][1].close();
  fout_road[1][0].close();
  fout_road[1][1].close();

  //Output road combination
  fstream fout_pair;
  fout_pair.open("road_pairs.txt", ios::out);

  for(std::list<Trigger>::iterator iter = triggers.begin(); iter != triggers.end(); ++iter)
    {
      fout_pair << iter->first << "  " << iter->second << endl;
    }
  fout_pair.close();
}

void TriggerAnalyzer::trimEvent(SRawEvent* rawEvent)
{
  acceptEvent(rawEvent, USE_HIT);

  for(int i = 0; i < 2; ++i)
    {
      for(std::list<TriggerRoad>::iterator iter = roads_found[i].begin(); iter != roads_found[i].end(); ++iter)
	{
	  for(int j = 0; j < 4; ++j)
	    {
	      rawEvent->setHitFlag(iter->detectorIDs[j], iter->elementIDs[j], 2);
	    }
	}
    }
}
