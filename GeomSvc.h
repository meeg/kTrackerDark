/*
GeomSvc.h

Definition of the class GeomSvc, which provide geometry service to the 
whole tracking system. 

For now the drift time is not used so as to avoid left-right ambiguity

Author: Kun Liu, liuk@fnal.gov
Created: 10-18-2011

Updated by Kun Liu on 07-03-2012
---  Support all detector planes including chambers, hodoscopes, and prop.
     tubes
---  The convention of detector plane ID is: firstly chambers, then hodos,
     and then prop. tubes
*/

#ifndef _GEOMSVC_H
#define _GEOMSVC_H

#include "MODE_SWITCH.h"

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>

#include <TSpline.h>

class GeomSvc
{
public:
  ///singlton instance
  static GeomSvc* instance();

  ///Initialization, either from MySQL or from ascii file
  void init(std::string geometrySchema);
  void loadCalibration(std::string calibrateFile);
  void loadAlignment(std::string alignmentFile_chamber, std::string alignmentFile_hodo, std::string alignmentFile_prop);
  void loadMilleAlignment(std::string alignmentFile_mille);

  ///Close the geometry service before exit or starting a new one
  void close();

  ///Convert the official detectorName to local detectorName
  void toLocalDetectorName(std::string& detectorName, int& eID);

  ///Get the plane position
  int getDetectorID(std::string detectorName) { return map_detectorID[detectorName]; }
  std::string getDetectorName(int detectorID) { return map_detectorName[detectorID]; }
  std::vector<int> getDetectorIDs(std::string pattern);
  bool findPatternInDetector(int detectorID, std::string pattern);

  double getPlanePosition(int planeID) { return z0[planeID]; }
  double getPlaneSpacing(int planeID) { return spacing[planeID]; }
  double getCellWidth(int planeID) { return cellWidth[planeID]; }
  double getCostheta(int planeID) { return costheta[planeID]; }
  double getSintheta(int planeID) { return sintheta[planeID]; }
  double getPlaneCenterX(int planeID) { return x0[planeID]; }
  double getPlaneCenterY(int planeID) { return y0[planeID]; }
  double getPlaneScaleX(int planeID) { return x2[planeID] - x1[planeID]; }
  double getPlaneScaleY(int planeID) { return y2[planeID] - y1[planeID]; }
  int getPlaneNElements(int planeID) { return nElements[planeID]; }
  double getPlaneResolution(int planeID) { return resolution[planeID]; }
  
  double getPlaneZOffset(int planeID) {return offset_z0[planeID]; }
  double getPlanePhiOffset(int planeID) { return offset_phi[planeID]; }
  double getPlaneWOffset(int planeID) { return offset_pos[planeID]; }
  double getPlaneWOffset(int planeID, int moduleID);

  int getPlaneType(int planeID);

  double getRotationInZ(int planeID) { return theta_z[planeID]; }
  double getRotationInX(int planeID) { return theta_x[planeID]; }
  double getRotationInY(int planeID) { return theta_y[planeID]; }

  double getKMAGCenter() { return (zmin_kmag + zmax_kmag)/2.; }
  double getKMAGUpstream() { return zmin_kmag; }
  double getKMAGDownstream() { return zmax_kmag; }

  ///Convert the planeID and elementID to the actual hit position
  void getMeasurement(int planeID, int elementID, double& measurement, double& dmeasurement);
  double getMeasurement(int planeID, int elementID);
  void get2DBoxSize(int planeID, int elementID, double& x_min, double& x_max, double& y_min, double& y_max);
  int getExpElementID(int planeID, double pos_exp);

  bool calibrationLoaded() { return calibration_loaded; }
  double getDriftDistance(int planeID, double tdcTime);
  TSpline3 *getRTCurve(int planeID) { return rtprofile[planeID-1]; }

  ///Convert the stereo hits to Y value
  double getYinStereoPlane(int planeID, double x, double u) { return u/sintheta[planeID] - x/tantheta[planeID]; }
  double getUinStereoPlane(int planeID, double x, double y) { return x*costheta[planeID] + y*sintheta[planeID]; } 
  double getXinStereoPlane(int planeID, double u, double y) { return u/costheta[planeID] - y*tantheta[planeID]; } 

  ///See if a point is in a plane
  bool isInPlane(int planeID, double x, double y);
  bool isInKMAG(double x, double y);

  ///Debugging print of the content
  void print();
  void printAlignPar();
  void printTable();
  void printWirePosition();

private:

  ///Parameters
  double spacing[nChamberPlanes+nHodoPlanes+nPropPlanes+1];
  double cellWidth[nChamberPlanes+nHodoPlanes+nPropPlanes+1];
  double xoffset[nChamberPlanes+nHodoPlanes+nPropPlanes+1];
  double overlap[nChamberPlanes+nHodoPlanes+nPropPlanes+1];
  int nElements[nChamberPlanes+nHodoPlanes+nPropPlanes+1];
  double angleFromVert[nChamberPlanes+nHodoPlanes+nPropPlanes+1];
  double sintheta[nChamberPlanes+nHodoPlanes+nPropPlanes+1];
  double costheta[nChamberPlanes+nHodoPlanes+nPropPlanes+1];
  double tantheta[nChamberPlanes+nHodoPlanes+nPropPlanes+1];
  double z0[nChamberPlanes+nHodoPlanes+nPropPlanes+1];
  double x1[nChamberPlanes+nHodoPlanes+nPropPlanes+1], x2[nChamberPlanes+nHodoPlanes+nPropPlanes+1], x0[nChamberPlanes+nHodoPlanes+nPropPlanes+1];
  double y1[nChamberPlanes+nHodoPlanes+nPropPlanes+1], y2[nChamberPlanes+nHodoPlanes+nPropPlanes+1], y0[nChamberPlanes+nHodoPlanes+nPropPlanes+1];

  ///Alignment parameters
  //For chambers and hodoscopes, it is defined by plane
  double offset_pos[nChamberPlanes+nHodoPlanes+nPropPlanes+1];
  //For prop. tubes, it is defined by module
  double offset_pos_prop[4][9];   //4 planes, each has 9 modules
  double resolution[nChamberPlanes+nHodoPlanes+nPropPlanes+1];

  double theta_x[nChamberPlanes+nHodoPlanes+nPropPlanes+1];
  double theta_y[nChamberPlanes+nHodoPlanes+nPropPlanes+1];
  double theta_z[nChamberPlanes+nHodoPlanes+nPropPlanes+1];

  //Alignment from millepede
  double offset_phi[nChamberPlanes+nHodoPlanes+nPropPlanes+1];
  double offset_z0[nChamberPlanes+nHodoPlanes+nPropPlanes+1];

  //Calibration parameters
  bool calibration_loaded;
  TSpline3 *rtprofile[24];

  double xmin_kmag, xmax_kmag;
  double ymin_kmag, ymax_kmag;
  double zmin_kmag, zmax_kmag;

  std::map<std::string, int> map_detectorID;
  std::map<int, std::string> map_detectorName;
  std::map<std::pair<int, int>, double> map_wirePosition;

  static GeomSvc* p_geometrySvc;
};

#endif
