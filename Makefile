ROOTCONFIG   := root-config
ROOTCINT     := rootcint
ARCH         := $(shell $(ROOTCONFIG) --arch)
ROOTCFLAGS   := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS  := $(shell $(ROOTCONFIG) --ldflags)
ROOTGLIBS    := $(shell $(ROOTCONFIG) --libs)

G4CONFIG     := geant4-config
G4CFLAGS     := $(shell $(G4CONFIG) --cflags) 
G4LDFLAGS    := $(shell $(G4CONFIG) --libs) 

MYSQLCONIG   := mysql_config
MYSQLCFLAGS  := $(shell $(MYSQLCONIG) --include)
MYSQLLDFLAGS := $(shell $(MYSQLCONIG) --libs)

CXX           = g++
CXXFLAGS      = -O3 -Wall -fPIC
LD            = g++
LDFLAGS       = -O3
SOFLAGS       = -shared

ifeq ($(ARCH),macosx64)
FORTRAN       = gfortran
else
FORTRAN       = g++
endif
FFLAGS        = -fPIC -m64 -pthread
FLFLAGS       = -lgfortran

CXXFLAGS     += $(ROOTCFLAGS)
LDFLAGS      += $(ROOTLDFLAGS) $(ROOTGLIBS)

CXXFLAGS     += -I$(BOOST)

CXXFLAGS     += $(G4CFLAGS)
LDFLAGS      += $(G4LDFLAGS)

CXXFLAGS     += $(MYSQLCFLAGS)
LDFLAGS      += $(MYSQLLDFLAGS)

SRAWEVENTO    = SRawEvent.o SRawEventDict.o
SRECEVENTO    = SRecEvent.o SRecEventDict.o
GEOMSVCO      = GeomSvc.o
MYSQLSVCO     = MySQLSvc.o
JOBOPTSSVCO   = JobOptsSvc.o
KALMANUTILO   = KalmanUtil.o
KALMANFILTERO = KalmanFilter.o
KALMANTRACKO  = KalmanTrack.o
KALMANFITTERO = KalmanFitter.o
FASTTRACKLETO = FastTracklet.o FastTrackletDict.o
KALMANFASTO   = KalmanFastTracking.o
VERTEXFITO    = VertexFit.o

SMPUTILO      = SMillepedeUtil.o SMillepedeUtilDict.o
SMILLEPEDEO   = SMillepede.o
MILLEPEDEO    = millepede.o
MILLEPEDES    = millepede.f

TRIGGERROADO  = TriggerRoad.o TriggerRoadDict.o
TRIGGERANALYZERO = TriggerAnalyzer.o

KFASTTRACKO   = kFastTracking.o
KFASTTRACK    = kFastTracking

KONLINETRACKO   = kOnlineTracking.o
KONLINETRACK    = kOnlineTracking

KTRACKERO     = kTracker.o
KTRACKER      = kTracker

KVERTEXO      = kVertex.o 
KVERTEX       = kVertex

MILLEALIGNO   = milleAlign.o
MILLEALIGN    = milleAlign

KTRACKERSO    = $(KTRACKER_LIB)/libkTracker.so
SRAWEVENTSO   = $(KTRACKER_LIB)/libSRawEvent.so

TRKEXTOBJS    = TrackExtrapolator/TrackExtrapolator.o TrackExtrapolator/DetectorConstruction.o TrackExtrapolator/Field.o TrackExtrapolator/TabulatedField3D.o \
		TrackExtrapolator/Settings.o TrackExtrapolator/GenericSD.o TrackExtrapolator/MCHit.o TrackExtrapolator/TPhysicsList.o 
CLASSOBJS     = $(GEOMSVCO) $(JOBOPTSSVCO) $(SRAWEVENTO) $(SRECEVENTO) $(KALMANUTILO) $(KALMANFILTERO) $(KALMANTRACKO) $(KALMANFITTERO) $(VERTEXFITO) \
		$(KALMANFASTO) $(FASTTRACKLETO) $(MYSQLSVCO) $(TRIGGERROADO) $(TRIGGERANALYZERO)
ALIGNOBJS     = $(SMPUTILO) $(SMILLEPEDEO) $(MILLEPEDEO)
OBJS          = $(CLASSOBJS) $(ALIGNOBJS) $(KVERTEXO) $(KTRACKERMULO) $(KSEEDERO) $(KVERTEXMO) $(KFASTTRACKO) $(KONLINETRACKO) $(MILLEALIGNO) $(KTRACKERO)
SLIBS         = $(KTRACKERSO) $(SRAWEVENTSO)
PROGRAMS      = $(KVERTEX) $(MILLEALIGN) $(KFASTTRACK) $(KONLINETRACK) $(KTRACKER)

all:            $(PROGRAMS) $(SLIBS)

.SUFFIXES: .cxx .o

$(MILLEPEDEO): $(MILLEPEDES) 
	$(FORTRAN) $(FFLAGS) -c $^ -o $@

$(KTRACKERSO):  $(CLASSOBJS) $(TRKEXTOBJS)
	$(LD) $^ -o $@  $(SOFLAGS) $(LDFLAGS) 
	@echo "$@ done."

$(SRAWEVENTSO):  $(SRAWEVENTO)
	$(LD) $^ -o $@  $(SOFLAGS) $(LDFLAGS)
	@echo "$@ done."

$(KTRACKERMUL):   $(KTRACKERMULO) $(CLASSOBJS) $(TRKEXTOBJS)
	$(LD) $^ -o $@ $(LDFLAGS) 
	@echo "$@ done."

$(KFASTTRACK):   $(KFASTTRACKO) $(CLASSOBJS) $(TRKEXTOBJS)
	$(LD) $^ -o $@ $(LDFLAGS) 
	@echo "$@ done."

$(KONLINETRACK):   $(KONLINETRACKO) $(CLASSOBJS) $(TRKEXTOBJS)
	$(LD) $^ -o $@ $(LDFLAGS) 
	@echo "$@ done."

$(KTRACKER):   $(KTRACKERO) $(CLASSOBJS) $(TRKEXTOBJS)
	$(LD) $^ -o $@ $(LDFLAGS) 
	@echo "$@ done."

$(KVERTEX):   $(KVERTEXO) $(CLASSOBJS) $(TRKEXTOBJS)
	$(LD) $^ -o $@ $(LDFLAGS) 
	@echo "$@ done."

$(MILLEALIGN):   $(MILLEALIGNO) $(CLASSOBJS) $(ALIGNOBJS) $(TRKEXTOBJS)
	$(LD) $^ -o $@ $(LDFLAGS) $(FLFLAGS) 
	@echo "$@ done."

.SUFFIXES: .cxx

SRawEventDict.cxx: SRawEvent.h SRawEventLinkDef.h
	@echo "Generating dictionary for $@ ..."
	$(ROOTCINT) -f $@ -c $^

SRecEventDict.cxx: SRecEvent.h SRecEventLinkDef.h
	@echo "Generating dictionary for $@ ..."
	$(ROOTCINT) -f $@ -c $^

SMillepedeUtilDict.cxx: SMillepedeUtil.h SMillepedeUtilLinkDef.h
	@echo "Generating dictionary for $@ ..."
	$(ROOTCINT) -f $@ -c $^

FastTrackletDict.cxx: FastTracklet.h FastTrackletLinkDef.h
	@echo "Generating dictionary for $@ ..."
	$(ROOTCINT) -f $@ -c $^

TriggerRoadDict.cxx: TriggerRoad.h TriggerRoadLinkDef.h
	@echo "Generating dictionary for $@ ..."
	$(ROOTCINT) -f $@ -c $^

.cxx.o:
	$(CXX) $(CXXFLAGS) -c $<

.PHONY: clean

clean:
	@echo "Cleanning everything ... "
	@rm $(PROGRAMS) $(OBJS) $(SLIBS) *Dict.cxx *Dict.h
