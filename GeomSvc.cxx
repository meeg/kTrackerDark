/*
GeomSvc.cxx

Implementation of class GeomSvc.

Author: Kun Liu, liuk@fnal.gov
Created: 10-19-2011
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

#include <TROOT.h>
#include <TTree.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TMath.h>
#include <TString.h>

#include "GeomSvc.h"

GeomSvc* GeomSvc::p_geometrySvc = NULL;

GeomSvc* GeomSvc::instance()
{
  if(p_geometrySvc == NULL)
    {
      p_geometrySvc = new GeomSvc;
    }

  return p_geometrySvc;
}

void GeomSvc::close()
{
  if(!p_geometrySvc)
    {
      for(int i = 0; i < 24; i++)
	{
	  delete rtprofile[i];
	}

      delete p_geometrySvc;
    }
  else
    {
      std::cout << "Error: no instance of geometry service found! " << std::endl;
    }
}

void GeomSvc::init(std::string geometrySchema)
{
  using namespace std;

  //Initialize the detectorID --- detectorName convention
  typedef std::map<std::string, int>::value_type nameToID;

  map_detectorID.insert(nameToID("D1U", 1));
  map_detectorID.insert(nameToID("D1Up", 2));
  map_detectorID.insert(nameToID("D1X", 3));
  map_detectorID.insert(nameToID("D1Xp", 4));
  map_detectorID.insert(nameToID("D1V", 5));
  map_detectorID.insert(nameToID("D1Vp", 6));
  map_detectorID.insert(nameToID("D2V", 7));
  map_detectorID.insert(nameToID("D2Vp", 8));
  map_detectorID.insert(nameToID("D2Xp", 9));
  map_detectorID.insert(nameToID("D2X", 10));
  map_detectorID.insert(nameToID("D2U", 11));
  map_detectorID.insert(nameToID("D2Up", 12));
  map_detectorID.insert(nameToID("D3pVp", 13));
  map_detectorID.insert(nameToID("D3pV", 14));
  map_detectorID.insert(nameToID("D3pXp", 15));
  map_detectorID.insert(nameToID("D3pX", 16));
  map_detectorID.insert(nameToID("D3pUp", 17));
  map_detectorID.insert(nameToID("D3pU", 18));
  map_detectorID.insert(nameToID("D3mU", 19));
  map_detectorID.insert(nameToID("D3mUp", 20));
  map_detectorID.insert(nameToID("D3mX", 21));
  map_detectorID.insert(nameToID("D3mXp", 22));
  map_detectorID.insert(nameToID("D3mV", 23));
  map_detectorID.insert(nameToID("D3mVp", 24));
  
  map_detectorID.insert(nameToID("H1X1", 25));
  map_detectorID.insert(nameToID("H1X2", 26));
  map_detectorID.insert(nameToID("H1Y1", 27));
  map_detectorID.insert(nameToID("H1Y2", 28));
  map_detectorID.insert(nameToID("H2Y1", 29));
  map_detectorID.insert(nameToID("H2Y2", 30));
  map_detectorID.insert(nameToID("H2X1", 31));
  map_detectorID.insert(nameToID("H2X2", 32));
  map_detectorID.insert(nameToID("H3X1", 33));
  map_detectorID.insert(nameToID("H3X2", 34));
  map_detectorID.insert(nameToID("H4Y11", 35));
  map_detectorID.insert(nameToID("H4Y12", 36));
  map_detectorID.insert(nameToID("H4Y21", 37));
  map_detectorID.insert(nameToID("H4Y22", 38));
  map_detectorID.insert(nameToID("H4X1", 39));
  map_detectorID.insert(nameToID("H4X2", 40));

  map_detectorID.insert(nameToID("P1Y1", 41));
  map_detectorID.insert(nameToID("P1Y2", 42));
  map_detectorID.insert(nameToID("P1X1", 43));
  map_detectorID.insert(nameToID("P1X2", 44));
  map_detectorID.insert(nameToID("P2X1", 45));
  map_detectorID.insert(nameToID("P2X2", 46));
  map_detectorID.insert(nameToID("P2Y1", 47));
  map_detectorID.insert(nameToID("P2Y2", 48));

  typedef std::map<int, std::string>::value_type idToName;
  for(std::map<std::string, int>::iterator iter = map_detectorID.begin(); iter != map_detectorID.end(); ++iter)
    {
      map_detectorName.insert(idToName(iter->second, iter->first));
    }

  //Initialize the geometrical variables which should be from MySQL database
  for(int i = 1; i <= nChamberPlanes+nHodoPlanes+nPropPlanes; i++)
    {
      nElements[i] = 0;
      x0[i] = 0.;
      y0[i] = 0.;
      z0[i] = 0.;

      x1[i] = 1E6;
      y1[i] = 1E6;
      x2[i] = -1E6;
      y2[i] = -1E6;
    }

  //Connect server
  char serverName[200];
  sprintf(serverName, "mysql://%s", MYSQL_SERVER);
  TSQLServer* con = TSQLServer::Connect(serverName, "seaguest","qqbar2mu+mu-");
  
  //Make query to Planes table
  char query[300];
  const char* buf_planes = "SELECT detectorName,spacing,cellWidth,overlap,numElements,angleFromVert,"
    "xPrimeOffset,x0,y0,z0,planeWidth,planeHeight,theta_x,theta_y,theta_z from %s.Planes WHERE"
    " detectorName LIKE 'D%%' OR detectorName LIKE 'H__' OR detectorName LIKE 'H____' OR "
    "detectorName LIKE 'P____'";
  sprintf(query, buf_planes, geometrySchema.c_str());
  TSQLResult *res = con->Query(query);

  unsigned int nRows = res->GetRowCount();
  int dummy = 0;
  for(unsigned int i = 0; i < nRows; i++)
    {
      TSQLRow *row = res->Next();
      string detectorName(row->GetField(0));
      toLocalDetectorName(detectorName, dummy);   
       
      int detectorID = map_detectorID[detectorName];
      spacing[detectorID] = atof(row->GetField(1));
      cellWidth[detectorID] = atof(row->GetField(2));
      overlap[detectorID] = atof(row->GetField(3));
      angleFromVert[detectorID] = atof(row->GetField(5));
      xoffset[detectorID] = atof(row->GetField(6));
      theta_x[detectorID] = atof(row->GetField(12));
      theta_y[detectorID] = atof(row->GetField(13));
      theta_z[detectorID] = atof(row->GetField(14));

      //Following items need to be sumed or averaged over all modules
      nElements[detectorID] += atoi(row->GetField(4));
      double x0_i = atof(row->GetField(7));
      double y0_i = atof(row->GetField(8));
      double z0_i = atof(row->GetField(9));
      double width_i = atof(row->GetField(10));
      double height_i = atof(row->GetField(11));

      double x1_i = x0_i - 0.5*width_i;
      double x2_i = x0_i + 0.5*width_i;
      double y1_i = y0_i - 0.5*height_i;
      double y2_i = y0_i + 0.5*height_i;

      x0[detectorID] += x0_i;
      y0[detectorID] += y0_i;
      z0[detectorID] += z0_i;
      if(x1[detectorID] > x1_i) x1[detectorID] = x1_i;
      if(x2[detectorID] < x2_i) x2[detectorID] = x2_i;
      if(y1[detectorID] > y1_i) y1[detectorID] = y1_i;
      if(y2[detectorID] < y2_i) y2[detectorID] = y2_i;

      //Calculated value
      resolution[detectorID] = cellWidth[detectorID]/sqrt(12.);
      sintheta[detectorID] = sin(angleFromVert[detectorID] + theta_z[detectorID]);
      costheta[detectorID] = cos(angleFromVert[detectorID] + theta_z[detectorID]);
      tantheta[detectorID] = tan(angleFromVert[detectorID] + theta_z[detectorID]);

      delete row;      
    }
  delete res;
 
  for(int i = 41; i <= nChamberPlanes+nHodoPlanes+nPropPlanes; i++)
    {
      x0[i] = x0[i]/9.;
      y0[i] = y0[i]/9.;
      z0[i] = z0[i]/9.;
    }

  ///Initialize the alignment parameters
  for(int i = 1; i <= nChamberPlanes+nHodoPlanes+nPropPlanes; i++)
    {
      offset_pos[i] = 0.;
      offset_z0[i] = 0.;
      offset_phi[i] = 0.;
    }

  for(int i = 0; i < 4; i++)
    {
      for(int j = 0; j < 9; j++)
	{
	  offset_pos_prop[i][j] = 0.;
	}
    }

  //load the initial value in the planeOffsets table
  const char* buf_offsets = "SELECT detectorName,deltaX,deltaY,deltaZ,rotateAboutZ FROM %s.PlaneOffsets WHERE"
    " detectorName LIKE 'D%%' OR detectorName LIKE 'H__' OR detectorName LIKE 'H____' OR detectorName LIKE 'P____'";
  sprintf(query, buf_offsets, geometrySchema.c_str());
  res = con->Query(query);

  nRows = res->GetRowCount();
  for(unsigned int i = 0; i < nRows; ++i)
    {
      TSQLRow* row = res->Next();
      string detectorName(row->GetField(0));
      toLocalDetectorName(detectorName, dummy);

      double deltaX = atof(row->GetField(1));
      double deltaY = atof(row->GetField(2));
      double deltaZ = atof(row->GetField(3));
      double rotateAboutZ = atof(row->GetField(4));

      int detectorID = map_detectorID[detectorName];
      if(detectorID > 40) continue; //Temporarily, prop.tubes not implemented

      offset_z0[detectorID] = deltaZ;
      offset_phi[detectorID] = rotateAboutZ;
      
      sintheta[detectorID] = sin(angleFromVert[detectorID] + theta_z[detectorID] + offset_phi[detectorID]);
      costheta[detectorID] = cos(angleFromVert[detectorID] + theta_z[detectorID] + offset_phi[detectorID]);
      tantheta[detectorID] = tan(angleFromVert[detectorID] + theta_z[detectorID] + offset_phi[detectorID]);
      
      offset_pos[detectorID] = deltaX*costheta[detectorID] + deltaY*sintheta[detectorID];
      if(detectorID <= 24) resolution[detectorID] = RESOLUTION_DC; 
      
      delete row;
    }

  delete res;
  delete con;

  /////Here starts the user-defined part
  //load alignment parameters
  loadAlignment("alignment.txt", "alignment_hodo.txt", "alignment_prop.txt");
  loadMilleAlignment("align_mille.txt");
  calibration_loaded = false;

  ///Initialize the position look up table for all wires, hodos, and tubes
  typedef std::map<std::pair<int, int>, double>::value_type posType;
  for(int i = 1; i <= nChamberPlanes; i++)
    {
      for(int j = 1; j <= nElements[i]; j++)
	{
	  double pos = (j - (nElements[i]+1.)/2.)*spacing[i] + xoffset[i] + x0[i]*costheta[i] + y0[i]*sintheta[i] + offset_pos[i];
	  map_wirePosition.insert(posType(make_pair(i, j), pos));
	}
    }
  

  // 2. for hodoscopes and prop. tubes
  for(int i = nChamberPlanes + 1; i <= nChamberPlanes+nHodoPlanes+nPropPlanes; i++)
    {
      for(int j = 1; j <= nElements[i]; j++)
	{
	  double pos;
	  if(i <= nChamberPlanes+nHodoPlanes)
	    {
	      pos = x0[i]*costheta[i] + y0[i]*sintheta[i] + xoffset[i] + (j - (nElements[i]+1)/2.)*spacing[i] + offset_pos[i];
	    }
	  else
	    {
	      int splaneID = int((i - nChamberPlanes - nHodoPlanes - 1)/2);
	      int moduleID = int((j - 1)/8);        //Need to re-define moduleID for run2, note it's reversed compared to elementID
	      pos = x0[i]*costheta[i] + y0[i]*sintheta[i] + xoffset[i] + (j - (nElements[i]+1)/2.)*spacing[i] + offset_pos_prop[splaneID][moduleID];
	    }
	  map_wirePosition.insert(posType(make_pair(i, j), pos));
	}
    } 

  ///Initialize channel mapping  --- not needed at the moment
  //xmin_kmag = -304.8;
  //xmax_kmag = 304.8;
  //ymin_kmag = -228.6;
  //ymax_kmag = 228.6;
  xmin_kmag = -57.*2.54;
  xmax_kmag = 57.*2.54;
  ymin_kmag = -40.*2.54;
  ymax_kmag = 40.*2.54;

  zmin_kmag = 1064.26 - 120.*2.54;
  zmax_kmag = 1064.26 + 120.*2.54;
}

std::vector<int> GeomSvc::getDetectorIDs(std::string detectorName)
{
  std::vector<int> detectorIDs;
  detectorIDs.clear();

  for(std::map<std::string, int>::iterator iter = map_detectorID.begin(); iter != map_detectorID.end(); ++iter)
    {
      if((*iter).first.find(detectorName) != std::string::npos)
	{
	  detectorIDs.push_back((*iter).second);
	}
    }

  return detectorIDs;
}
     
bool GeomSvc::isInPlane(int detectorID, double x, double y)
{
  if(x < x1[detectorID] || x > x2[detectorID]) return false;
  if(y < y1[detectorID] || y > y2[detectorID]) return false;

  return true;
}

bool GeomSvc::isInKMAG(double x, double y)
{
  if(x < xmin_kmag || x > xmax_kmag) return false;
  if(y < ymin_kmag || y > ymax_kmag) return false;

  return true;
}

void GeomSvc::getMeasurement(int detectorID, int elementID, double& measurement, double& dmeasurement)
{
  measurement = map_wirePosition[std::make_pair(detectorID, elementID)];
  dmeasurement = resolution[detectorID];
}

double GeomSvc::getMeasurement(int detectorID, int elementID)
{
  return map_wirePosition[std::make_pair(detectorID, elementID)];
}

int GeomSvc::getExpElementID(int detectorID, double pos_exp)
{
  int elementID = -1;
  for(int i = 1; i < nElements[detectorID]; i++)
    {
      double pos = map_wirePosition[std::make_pair(detectorID, i)];
      if(fabs(pos - pos_exp) < 0.5*spacing[detectorID])
	{
	  elementID = i;
	  break;
	}
    }

  return elementID;
}

void GeomSvc::get2DBoxSize(int detectorID, int elementID, double& x_min, double& x_max, double& y_min, double& y_max)
{
  std::string detectorName = getDetectorName(detectorID);
  if(detectorName.find("X") != std::string::npos)
    {
      double x_center = map_wirePosition[std::make_pair(detectorID, elementID)];
      double x_width = 0.5*(spacing[detectorID] + overlap[detectorID]);
      x_min = x_center - x_width;
      x_max = x_center + x_width;

      y_min = y1[detectorID];
      y_max = y2[detectorID];
    }
  else
    {
      double y_center = map_wirePosition[std::make_pair(detectorID, elementID)];
      double y_width = 0.5*(spacing[detectorID] + overlap[detectorID]);
      y_min = y_center - y_width;
      y_max = y_center + y_width;

      x_min = x1[detectorID];
      x_max = x2[detectorID];
    }
}

int GeomSvc::getPlaneType(int planeID)
{
  std::string detectorName = getDetectorName(planeID);
  if(detectorName.find("X") != std::string::npos) return 1;
  if(detectorName.find("U") != std::string::npos) return 2;
  if(detectorName.find("V") != std::string::npos) return 3;

  return 0;
}

void GeomSvc::toLocalDetectorName(std::string& detectorName, int& eID)
{
  using namespace std;

  if(detectorName.find("P") != string::npos)
    {
      string XY = detectorName[2] == 'H' ? "Y" : "X";
      string FB = (detectorName[3] == 'f' || detectorName[4] == 'f') ? "1" : "2"; //temporary solution
      
      detectorName.replace(2, detectorName.length(), "");
      detectorName += XY;
      detectorName += FB;
     
      int moduleID = std::isdigit(detectorName[3]) == 1 ? atoi(&detectorName[3]) : atoi(&detectorName[4]); 
      if(eID <= 8) //Means either it's at lowest module or elementID is only from 1 to 8, repeatedly
	{
	  eID = (9 - moduleID)*8 + eID;
	}	  
    }
  else if(detectorName.find("H") != string::npos)
    {
      string suffix;
      if(detectorName.find("B") != string::npos || detectorName.find("L") != string::npos)
	{
	  suffix = "1";
	}
      else if(detectorName.find("T") != string::npos || detectorName.find("R") != string::npos)
	{
	  suffix = "2";
	}

      if(detectorName.find("H4Y") != string::npos)
	{
	  detectorName.replace(4, detectorName.length(), "");
	}
      else if(detectorName.find("T") != string::npos || detectorName.find("B") != string::npos)
	{
	  detectorName.replace(2, detectorName.length(), "X");
	}
      else if(detectorName.find("L") != string::npos || detectorName.find("R") != string::npos)
	{
	  detectorName.replace(2, detectorName.length(), "Y");
	}

      detectorName += suffix;
    }
}

double GeomSvc::getDriftDistance(int planeID, double tdcTime)
{
  if(!calibration_loaded) return -1.;
  if(planeID <= 24)
    {
      return rtprofile[planeID-1]->Eval(tdcTime);
    }
  
  return 0.;
}

void GeomSvc::loadAlignment(std::string alignmentFile_chamber, std::string alignmentFile_hodo, std::string alignmentFile_prop)
{
  using namespace std;

  //load alignment numbers for chambers
  char filename[300];
  sprintf(filename, "%s/%s", KTRACKER_ROOT, alignmentFile_chamber.c_str());
  fstream _align_chamber;
  _align_chamber.open(filename, ios::in);
  
  char buf[300];
  if(_align_chamber)
    {
      for(int i = 1; i <= nChamberPlanes; i++)
	{
	  _align_chamber.getline(buf, 100);
	  istringstream stringBuf(buf);

	  stringBuf >> offset_pos[i] >> resolution[i];
          if(resolution[i] < RESOLUTION_DC) resolution[i] = RESOLUTION_DC;
	}	  
    }
  else
    {
      cout << "Alignment file for chamber is not found! " << endl;
      cout << "Will use survey numbers instead! " << endl;
    }
  _align_chamber.close();

  for(int i = 1; i <= nChamberPlanes; i += 2)
    {
      double resol = resolution[i] > resolution[i+1] ? resolution[i] : resolution[i+1];
      resolution[i] = resol;
      resolution[i+1] = resol;
    }

  //load alignment numbers for hodos
  fstream _align_hodo;
  sprintf(filename, "%s/%s", KTRACKER_ROOT, alignmentFile_hodo.c_str());
  _align_hodo.open(filename, ios::in);
  
  if(_align_hodo)
    {
      for(int i = nChamberPlanes+1; i <= 40; i++)
	{
	  _align_hodo.getline(buf, 100);
	  istringstream stringBuf(buf);

	  stringBuf >> offset_pos[i];
	}	  
    }
  else
    {
      cout << "Alignment file for hodo is not found! " << endl;
      cout << "Will use ideal numbers instead! " << endl;
    }
  _align_hodo.close();

  //load alignment numbers for prop. tubes
  fstream _align_prop;
  sprintf(filename, "%s/%s", KTRACKER_ROOT, alignmentFile_prop.c_str());
  _align_prop.open(filename, ios::in);

  if(_align_prop)
    {
      for(int i = 0; i < 4; i++)
	{
	  for(int j = 0; j < 9; j++)
	    {
	      _align_prop.getline(buf, 100);
	      istringstream stringBuf(buf);

	      stringBuf >> offset_pos_prop[i][j];
	    }
	}

      for(int i = 41; i <= 48; i++)
	{
	  spacing[i] = 5.08;
	}
    }
  else
    {
      cout << "Alignment file for prop. tubes is not found! " << endl;
      cout << "Will use survey numbers instead! " << endl;
    }
}

void GeomSvc::loadMilleAlignment(std::string alignmentFile_mille)
{
  using namespace std;

  //load alignment numbers for chambers
  char filename[300];
  sprintf(filename, "%s/%s", KTRACKER_ROOT, alignmentFile_mille.c_str());
  fstream _align_mille;
  _align_mille.open(filename, ios::in);
  
  char buf[300];
  if(_align_mille)
    {
      for(int i = 1; i <= nChamberPlanes; i++)
	{
	  _align_mille.getline(buf, 100);
	  istringstream stringBuf(buf);

	  stringBuf >> offset_z0[i] >> offset_phi[i] >> offset_pos[i] >> resolution[i];

	  z0[i] += offset_z0[i];
	  sintheta[i] = sin(angleFromVert[i] + theta_z[i] + offset_phi[i]);
	  costheta[i] = cos(angleFromVert[i] + theta_z[i] + offset_phi[i]);
	  tantheta[i] = tan(angleFromVert[i] + theta_z[i] + offset_phi[i]);

	  //if(resolution[i] < spacing[i]/sqrt(12.)/10.) resolution[i] = spacing[i]/sqrt(12.)/10.;
	  if(resolution[i] < RESOLUTION_DC) resolution[i] = RESOLUTION_DC;
	}	  
    }
  else
    {
      cout << "Mille alignment file for chamber is not found! " << endl;
      cout << "Will use old numbers instead! " << endl;
    }
  _align_mille.close();

  for(int i = 1; i <= nChamberPlanes; i+=2)
    {
      resolution[i] = 0.5*(resolution[i] + resolution[i+1]);
      resolution[i+1] = resolution[i];
    }
}

void GeomSvc::loadCalibration(std::string calibrationFile)
{
  using namespace std;
 
  char filename[300];
  sprintf(filename, "%s/%s", KTRACKER_ROOT, calibrationFile.c_str()); 
  fstream _cali_file;
  _cali_file.open(filename, ios::in);
 
  char buf[300];
  int iBin, nBin, detectorID;
  double R[200], T[200];
  if(_cali_file)
    {
      calibration_loaded = true;

      while(_cali_file.getline(buf, 100))
	{
	  istringstream detector_info(buf);
	  detector_info >> detectorID >> nBin;

	  for(int i = 0; i < nBin; i++)
	    {
	      _cali_file.getline(buf, 100);
	      istringstream cali_line(buf);

	      cali_line >> iBin >> T[i] >> R[i];
	    }

	  rtprofile[detectorID-1] = new TSpline3(getDetectorName(detectorID).c_str(), T, R, nBin, "b1e1");
	}
    }
  _cali_file.close();
}

double GeomSvc::getPlaneWOffset(int planeID, int moduleID)
{
  int pID = (planeID - nChamberPlanes - nHodoPlanes - 1)/2;
  return offset_pos_prop[pID][moduleID];
}

void GeomSvc::print()
{
  std::cout << "Detector Name ----- Detector ID mapping" << std::endl;
  for(std::map<std::string, int>::iterator iter = map_detectorID.begin(); iter != map_detectorID.end(); ++iter)
    {
      std::cout << (*iter).first << " <-------> " << (*iter).second << std::endl;
      std::cout << "Spacing:       " << spacing[(*iter).second] << std::endl;
      std::cout << "CellWidth:     " << cellWidth[(*iter).second] << std::endl;
      std::cout << "overlap:       " << overlap[(*iter).second] << std::endl;
      std::cout << "Xoffset:       " << xoffset[(*iter).second] << std::endl;
      std::cout << "width:         " << getPlaneScaleX((*iter).second) << std::endl;
      std::cout << "height:        " << getPlaneScaleY((*iter).second) << std::endl;
      std::cout << "ScaleX:        " << x1[(*iter).second] << " === " << x2[(*iter).second] << std::endl;
      std::cout << "ScaleY:        " << y1[(*iter).second] << " === " << y2[(*iter).second] << std::endl;
      std::cout << "nElements:     " << nElements[(*iter).second] << std::endl;
      std::cout << "angleFromVert: " << angleFromVert[(*iter).second] << std::endl;
      std::cout << "Z position:    " << z0[(*iter).second] << std::endl;
      std::cout << "X position:    " << x0[(*iter).second] << std::endl;
      std::cout << "Y position:    " << y0[(*iter).second] << std::endl;
    }
}

void GeomSvc::printAlignPar()
{
  std::cout << "detectorID         DetectorName            offset_pos             offset_z             offset_phi" << std::endl;
  for(std::map<std::string, int>::iterator iter = map_detectorID.begin(); iter != map_detectorID.end(); ++iter)
    {
      std::cout << iter->second << "     " << iter->first << "    " << offset_pos[iter->second] << "     " << offset_z0[iter->second] << "      " << offset_phi[iter->second] << std::endl;
    }
}

void GeomSvc::printTable()
{
  std::cout << "detectorName    detectorID      Spacing     Xoffset     overlap     width       height       nElement       angleFromVertical       Z" << std::endl;
  for(std::map<std::string, int>::iterator iter = map_detectorID.begin(); iter != map_detectorID.end(); ++iter)
    {
      std::cout << (*iter).first << "     ";
      std::cout << (*iter).second << "     ";
      std::cout << spacing[(*iter).second] << "     ";
      std::cout << xoffset[(*iter).second] << "     ";
      std::cout << overlap[(*iter).second] << "     ";
      std::cout << getPlaneScaleX((*iter).second) << "     ";
      std::cout << getPlaneScaleY((*iter).second) << "     ";
      std::cout << nElements[(*iter).second] << "     ";
      std::cout << angleFromVert[(*iter).second] << "     ";
      std::cout << z0[(*iter).second] << "     ";
      std::cout << x0[(*iter).second] << "     ";
      std::cout << y0[(*iter).second] << std::endl;
    }
}
