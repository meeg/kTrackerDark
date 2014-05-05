ROOTCONFIG   := root-config
ROOTCINT     := rootcint
ROOTCFLAGS   := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS  := $(shell $(ROOTCONFIG) --ldflags)
ROOTGLIBS    := $(shell $(ROOTCONFIG) --libs)

G4CONFIG     := geant4-config
G4CFLAGS     := $(shell $(G4CONFIG) --cflags) 
G4LDFLAGS    := $(shell $(G4CONFIG) --libs)

MYSQLCONIG   := mysql_config
MYSQLCFLAGS  := $(shell $(MYSQLCONIG) --include)
MYSQLLDFLAGS := $(shell $(MYSQLCONIG) --libs)

CXX           = g++ -g
CXXFLAGS      = -O3 -Wall -fPIC
LD            = g++ -g
LDFLAGS       = -O3
SOFLAGS       = -shared

CXXFLAGS     += $(ROOTCFLAGS)
LDFLAGS      += $(ROOTLDFLAGS) $(ROOTGLIBS)

#CXXFLAGS     += -I$(CLHEP_INCLUDE_DIR)
#LDFLAGS      += -L$(CLHEP_LIB_DIR) -l$(CLHEP_LIB)

CXXFLAGS     += $(G4CFLAGS)
LDFLAGS      += $(G4LDFLAGS)

CXXFLAGS     += $(MYSQLCFLAGS)
LDFLAGS      += $(MYSQLLDFLAGS)

LDFLAGS      += -L$(KTRACKER_LIB)/ -lkTrackerSvc

TRKEXTO       = TrackExtrapolator.o DetectorConstruction.o Field.o TabulatedField3D.o Settings.o GenericSD.o MCHit.o TPhysicsList.o 
TRKEXTSO      = $(KTRACKER_LIB)/libTrkExt.so


TESTO         = TrkExt.o
TRKEXT        = TrkExt

SLIBS         = $(TRKEXTSO)
OBJS          = $(TRKEXTO) $(TESTO)
PROGRAMS      = $(TRKEXT)

all:            $(TRKEXT)

.SUFFIXES: .cc .o

$(TRKEXTSO):    $(TRKEXTO)
	$(LD) $^ -o $@ $(SOFLAGS) $(LDFLAGS) 
	@echo "$@ done."
#	@cp $@ $(KTRACKER_LIB)

$(TRKEXT):      $(TRKEXTSO) $(TESTO)
	$(LD) $^ -o $@ $(LDFLAGS)
	@echo "$@ done."

.cc.o:
	$(CXX) $(CXXFLAGS) -c $<

.PHONY: clean

clean:
	@echo "Cleanning everything ... "
	@rm $(PROGRAMS) $(OBJS) $(SLIBS)
