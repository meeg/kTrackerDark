CXX = g++
CXXFLAGS = -g -Wall -fPIC
ROOTFLAGS = `root-config --cflags --glibs`

SERVICES = $(KTRACKER_ROOT)/kTrackerServices

INCLUDE += -I$(KTRACKER_ROOT)/ 

LDLIBS += -L$(KTRACKER_ROOT) -lkTracker
LDLIBS += -L$(SERVICES) -lkTrackerSvc


BINARIES = TestJobOpts.o
TARGETS = TestJobOpts

#--- if using 'make all' ---#
all : $(TARGETS)

#--- if making individual targets ---#

TestJobOpts.o : TestJobOpts.cxx $(SERVICES)/JobOptsSvc.o
	$(CXX) $(INCLUDE) $(CXXFLAGS) $(ROOTFLAGS) -o $*.o $(LDLIBS) -c $*.cxx #compile
	$(CXX) $(INCLUDE) $(CXXFLAGS) $(ROOTFLAGS) $(LDLIBS) -o $* $*.o        #link

clean:
	rm -f $(BINARIES) $(TARGETS)

