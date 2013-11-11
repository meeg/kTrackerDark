ROOTCONFIG   := root-config
ROOTCINT     := rootcint
ARCH         := $(shell $(ROOTCONFIG) --arch)
ROOTCFLAGS   := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS  := $(shell $(ROOTCONFIG) --ldflags)
ROOTGLIBS    := $(shell $(ROOTCONFIG) --libs)

G4CONFIG     := geant4-config
G4CFLAGS     := $(shell $(G4CONFIG) --cflags) 
G4LDFLAGS    := $(shell $(G4CONFIG) --libs) 

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

LDFLAGS      += $(FLFLAGS)

CXXFLAGS     += $(ROOTCFLAGS)
LDFLAGS      += $(ROOTLDFLAGS) $(ROOTGLIBS)

CXXFLAGS     += $(G4CFLAGS)
LDFLAGS      += $(G4LDFLAGS)

CXXFLAGS     += -I$(MYSQL_INCLUDE)
LDFLAGS      += -lz -L$(MYSQL_LIB) -lmysqlclient

SRAWEVENTO    = SRawEvent.o SRawEventDict.o

SRECEVENTO    = SRecEvent.o SRecEventDict.o

GEOMSVCO      = GeomSvc.o

MYSQLSVCO     = MySQLSvc.o

SEEDFINDERO   = SeedFinder.o

KALMANUTILO   = KalmanUtil.o

KALMANFILTERO = KalmanFilter.o

KALMANTRACKO  = KalmanTrack.o

KALMANFINDERO = KalmanFinder.o

KALMANFITTERO = KalmanFitter.o

FASTTRACKLETO = FastTracklet.o FastTrackletDict.o

KALMANFASTO   = KalmanFastTracking.o

VERTEXFITO    = VertexFit.o

SMPUTILO      = SMillepedeUtil.o SMillepedeUtilDict.o
SMILLEPEDEO   = SMillepede.o
MILLEPEDEO    = millepede.o
MILLEPEDES    = millepede.f

KSEEDERO      = kSeeder.o
KSEEDER       = kSeeder

KTRACKERMULO  = kTracker.o
KTRACKERMUL   = kTracker

KFASTTRACKO   = kFastTracking.o
KFASTTRACK    = kFastTracking

KONLINETRACKO   = kOnlineTracking.o
KONLINETRACK    = kOnlineTracking

KVERTEXO      = kVertex.o 
KVERTEX       = kVertex

KVERTEXMO     = kVertex_mix.o 
KVERTEXM      = kVertex_mix

MILLEALIGNO   = milleAlign.o
MILLEALIGN    = milleAlign

KTRACKERSO    = libkTracker.so

TRKEXTOBJS    = TrackExtrapolator/TrackExtrapolator.o TrackExtrapolator/DetectorConstruction.o TrackExtrapolator/Field.o TrackExtrapolator/TabulatedField3D.o \
		TrackExtrapolator/Settings.o TrackExtrapolator/GenericSD.o TrackExtrapolator/MCHit.o TrackExtrapolator/TPhysicsList.o 
CLASSOBJS     = $(GEOMSVCO) $(SRAWEVENTO) $(SRECEVENTO) $(SEEDFINDERO) $(KALMANUTILO) $(KALMANFILTERO) $(KALMANTRACKO) $(KALMANFINDERO) $(KALMANFITTERO) $(VERTEXFITO) \
		$(SMPUTILO) $(SMILLEPEDEO) $(MILLEPEDEO) $(KALMANFASTO) $(FASTTRACKLETO) $(MYSQLSVCO)
OBJS          = $(CLASSOBJS) $(KVERTEXO) $(KTRACKERMULO) $(KSEEDERO) $(KVERTEXMO) $(KFASTTRACKO) $(KONLINETRACKO) $(MILLEALIGNO)
SLIBS         = $(KTRACKERSO)
PROGRAMS      = $(KSEEDER) $(KTRACKERMUL) $(KVERTEX) $(KVERTEXM) $(MILLEALIGN) $(KFASTTRACK) $(KONLINETRACK)

all:            $(PROGRAMS) $(SLIBS)

.SUFFIXES: .cxx .o

$(MILLEPEDEO): $(MILLEPEDES) 
	$(FORTRAN) $(FFLAGS) -c $^ -o $@

$(KTRACKERSO):  $(CLASSOBJS) $(TRKEXTOBJS)
	$(LD) $^ -o $@  $(SOFLAGS) $(LDFLAGS) 
	@echo "$@ done."

$(KSEEDER):   $(KSEEDERO) $(CLASSOBJS) $(TRKEXTOBJS)
	$(LD) $^ -o $@ $(LDFLAGS) 
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

$(KVERTEX):   $(KVERTEXO) $(CLASSOBJS) $(TRKEXTOBJS)
	$(LD) $^ -o $@ $(LDFLAGS) 
	@echo "$@ done."

$(KVERTEXM):   $(KVERTEXMO) $(CLASSOBJS) $(TRKEXTOBJS)
	$(LD) $^ -o $@ $(LDFLAGS) 
	@echo "$@ done."

$(MILLEALIGN):   $(MILLEALIGNO) $(CLASSOBJS) $(TRKEXTOBJS)
	$(LD) $^ -o $@ $(LDFLAGS) 
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

.cxx.o:
	$(CXX) $(CXXFLAGS) -c $<

.PHONY: clean

clean:
	@echo "Cleanning everything ... "
	@rm $(PROGRAMS) $(OBJS) $(SLIBS) *Dict.cxx *Dict.h
