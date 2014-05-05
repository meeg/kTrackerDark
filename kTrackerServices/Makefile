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

CXXFLAGS     += $(G4CFLAGS)
LDFLAGS      += $(G4LDFLAGS)

CXXFLAGS     += $(MYSQLCFLAGS)
LDFLAGS      += $(MYSQLLDFLAGS)

JOBOPTSSVCO   = JobOptsSvc.o

KTRACKERSVCSO    = $(KTRACKER_LIB)/libkTrackerSvc.so

OBJS    = $(JOBOPTSSVCO)
SRC     = JobOptsSvc.cxx
HEADERS = JobOptsSvc.h

all:           $(KTRACKERSVCSO)

.SUFFIXES: .cxx .o

$(OBJS): %.o: %.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(KTRACKERSVCSO):  $(OBJS)
	$(LD) $^ -o $@  $(SOFLAGS) $(LDFLAGS) 
	@echo "$@ done."

.PHONY: clean

clean:
	@echo "Cleanning everything ... "
	@rm $(OBJS) $(SLIBS) $(KTRACKERSVCSO)
