#ifndef _MODE_SWITCH_H
#define _MODE_SWITCH_H

//-------------- kTracker ROOT -------------------
#define KTRACKER_ROOT "/Users/liuk/currentWork/kTracker_dev"

//--------------- Mode controls ------------------
#define GEOMETRY_VERSION "geometry_R997"

//==== Enable/disable dimuon mode
#define DIMUON_MODE 0

//==== Turn on/off KMag
#define KMAG_ON 1

//==== Enable Kalman fitting in fast tracking and alignme
//#define _ENABLE_KF

//==== Enable massive debugging output
#define _DEBUG_ON
//#define _DEBUG_ON_LEVEL_2

//--------------- Fast tracking configuration ----
#define TX_MAX 0.1
#define TY_MAX 0.12
#define X0_MAX 80.
#define Y0_MAX 150.
#define INVP_MIN 0.01
#define INVP_MAX 0.1
#define PROB_LOOSE 0.0
#define PROB_TIGHT 0.00001
#define HIT_REJECT 3.

//--------------- Geometry setup -----------------
#define nChamberPlanes 24
#define nHodoPlanes 16
#define nPropPlanes 8

#define Z_KMAG_BEND 1064.26
#define Z_FMAG_BEND 251.4
#define Z_KFMAG_BEND 375.
#define PT_KICK_FMAG 2.909
#define PT_KICK_KMAG 0.4067
#define ELOSS_KFMAG 8.12
#define ELOSS_ABSORBER 1.81
#define Z_ST2 1347.36
#define Z_ABSORBER 2028.19

#define Z_REF 0.

#define RESOLUTION_DC 0.1

//-------------- Useful marcros -----------------
#define Log(message) std::cout << "DEBUG: " << __FILE__ << "  " << __LINE__ << "  " << __FUNCTION__ << " :::  " << message << std::endl
#define varName(x) #x

#endif
