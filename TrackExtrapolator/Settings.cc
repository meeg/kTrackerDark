#include "Settings.hh"
#include "../MODE_SWITCH.h"

Settings::Settings()
{
  //  These are the settings that will be used by the monte carlo if parameters are not specified.

  seed = 0;
  beamMomentum = 120*GeV;
  beamCurrent = 2e12;
  asciiFieldMap = true;
  generator = "gun";
  energyCut = 1.0*GeV;
  recordMethod = "hits";
  eventPos = "both";
  dimuonSource = "both";
  login = "seaguest";
  outputFileName = "test_default";
  password = "qqbar2mu+mu-";
  fMagName = "tab.Fmag";
  kMagName = "tab.Kmag";
  sqlServer = MYSQL_SERVER_ADDR;
  dimuonRepeat = 1;
  ironOn = true;
  trackingZCut = 400*cm;
  trackingEnergyCut = 1.0*GeV;
#if defined ALIGNMENT_MODE
  kMagMultiplier = 0.;
  fMagMultiplier = 0.;
#elif defined MC_MODE
  kMagMultiplier = 1.;
  fMagMultiplier = 1.;
#else
  kMagMultiplier = 1.;
  fMagMultiplier = 1.;
#endif
  geometrySchema = "geometry_R997";
  magnetSchema = "geometry_R996_magneticFields";
  target = 1;
  pythia_shower = true;
  bucket_size = 40000;
}

Settings::~Settings()
{
}
