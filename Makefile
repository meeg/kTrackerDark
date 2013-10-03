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
FORTRAN       = gfortran-4.2
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
SRAWEVENTS    = SRawEvent.cxx SRawEventDict.cxx

SRECEVENTO    = SRecEvent.o SRecEventDict.o
SRECEVENTS    = SRecEvent.cxx SRecEventDict.cxx

GEOMSVCO      = GeomSvc.o
GEOMSVCS      = GeomSvc.cxx

MYSQLSVCO     = MySQLSvc.o
MYSQLSVCS     = MySQLSvc.cxx

SEEDFINDERO   = SeedFinder.o
SEEDFINDERS   = SeedFinder.cxx

KALMANUTILO   = KalmanUtil.o
KALMANUTILS   = KalmanUtil.cxx

KALMANFILTERO = KalmanFilter.o
KALMANFILTERS = KalmanFilter.cxx

KALMANTRACKO  = KalmanTrack.o
KALMANTRACKS  = KalmanTrack.cxx

KALMANFINDERO = KalmanFinder.o
KALMANFINDERS = KalmanFinder.cxx

KALMANFITTERO = KalmanFitter.o
KALMANFITTERS = KalmanFitter.cxx

FASTTRACKLETO = FastTracklet.o FastTrackletDict.o
FASTTRACKLETS = FastTracklet.cxx FastTrackletDict.cxx

KALMANFASTO   = KalmanFastTracking.o
KALMANFASTS   = KalmanFastTracking.cxx

VERTEXFITO    = VertexFit.o
VERTEXFITS    = VertexFit.cxx

SMPUTILO      = SMillepedeUtil.o SMillepedeUtilDict.o
SMPUTILS      = SMillepedeUtil.cxx SMillepedeUtilDict.cxx
SMILLEPEDEO   = SMillepede.o
SMILLEPEDES   = SMillepede.cxx
MILLEPEDEO    = millepede.o
MILLEPEDES    = millepede.f

KSEEDERO      = kSeeder.o
KSEEDERS      = kSeeder.cxx
KSEEDER       = kSeeder

KTRACKERMULO  = kTracker.o
KTRACKERMULS  = kTracker.cxx
KTRACKERMUL   = kTracker

KFASTTRACKO   = kFastTracking.o
KFASTTRACKS   = kFastTracking.cxx
KFASTTRACK    = kFastTracking

KONLINETRACKO   = kOnlineTracking.o
KONLINETRACKS   = kOnlineTracking.cxx
KONLINETRACK    = kOnlineTracking

KVERTEXO      = kVertex.o 
KVERTEXS      = kVertex.cxx 
KVERTEX       = kVertex

KVERTEXFO     = kVertex_fast.o 
KVERTEXFS     = kVertex_fast.cxx 
KVERTEXF      = kVertex_fast

KVERTEXMO     = kVertex_mix.o 
KVERTEXMS     = kVertex_mix.cxx 
KVERTEXM      = kVertex_mix

MILLEALIGNO   = milleAlign.o
MILLEALIGNS   = milleAlign.cxx
MILLEALIGN    = milleAlign

KTRACKERSO    = libkTracker.so

TRKEXTOBJS    = TrackExtrapolator/TrackExtrapolator.o TrackExtrapolator/DetectorConstruction.o TrackExtrapolator/Field.o TrackExtrapolator/TabulatedField3D.o \
		TrackExtrapolator/Settings.o TrackExtrapolator/GenericSD.o TrackExtrapolator/MCHit.o TrackExtrapolator/TPhysicsList.o 
CLASSOBJS     = $(GEOMSVCO) $(SRAWEVENTO) $(SRECEVENTO) $(SEEDFINDERO) $(KALMANUTILO) $(KALMANFILTERO) $(KALMANTRACKO) $(KALMANFINDERO) $(KALMANFITTERO) $(VERTEXFITO) \
		$(SMPUTILO) $(SMILLEPEDEO) $(MILLEPEDEO) $(KALMANFASTO) $(FASTTRACKLETO) $(MYSQLSVCO)
OBJS          = $(CLASSOBJS) $(KVERTEXO) $(KTRACKERMULO) $(KSEEDERO) $(KVERTEXMO) $(KFASTTRACKO) $(KVERTEXFO) $(KONLINETRACKO) $(MILLEALIGNO)
SLIBS         = $(KTRACKERSO)
PROGRAMS      = $(KSEEDER) $(KTRACKERMUL) $(KVERTEX) $(KVERTEXM) $(MILLEALIGN) $(KFASTTRACK) $(KVERTEXF) $(KONLINETRACK)

all:            $(PROGRAMS) $(SLIBS)

.SUFFIXES: .cxx .o

$(MILLEPEDEO): $(MILLEPEDES) 
	$(FORTRAN) $(FFLAGS) -c $^ -o $@

$(KTRACKERSO):  $(CLASSOBJS) $(TRKEXTOBJS)
	$(LD) $(SOFLAGS) $(LDFLAGS) $^ -o $@
	@echo "$@ done."

$(KSEEDER):   $(KSEEDERO) $(CLASSOBJS) $(TRKEXTOBJS)
	$(LD) $(LDFLAGS) $^ -o $@ 
	@echo "$@ done."

$(KTRACKERMUL):   $(KTRACKERMULO) $(CLASSOBJS) $(TRKEXTOBJS)
	$(LD) $(LDFLAGS) $^ -o $@ 
	@echo "$@ done."

$(KFASTTRACK):   $(KFASTTRACKO) $(CLASSOBJS) $(TRKEXTOBJS)
	$(LD) $(LDFLAGS) $^ -o $@ 
	@echo "$@ done."

$(KONLINETRACK):   $(KONLINETRACKO) $(CLASSOBJS) $(TRKEXTOBJS)
	$(LD) $(LDFLAGS) $^ -o $@ 
	@echo "$@ done."

$(KVERTEX):   $(KVERTEXO) $(CLASSOBJS) $(TRKEXTOBJS)
	$(LD) $(LDFLAGS) $^ -o $@ 
	@echo "$@ done."

$(KVERTEXM):   $(KVERTEXMO) $(CLASSOBJS) $(TRKEXTOBJS)
	$(LD) $(LDFLAGS) $^ -o $@ 
	@echo "$@ done."

$(KVERTEXF):   $(KVERTEXFO) $(CLASSOBJS) $(TRKEXTOBJS)
	$(LD) $(LDFLAGS) $^ -o $@ 
	@echo "$@ done."

$(MILLEALIGN):   $(MILLEALIGNO) $(CLASSOBJS) $(TRKEXTOBJS)
	$(LD) $(LDFLAGS) $^ -o $@ 
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
