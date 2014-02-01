XX = g++
CXXFLAGS = -g -Wall -fPIC
ROOTFLAGS = `root-config --cflags --glibs`
INCLUDE += -I$(KTRACKER_ROOT)/ 
LDLIBS += -L$(KTRACKER_ROOT) -lkTracker

BINARIES = TestJobOpts.o
TARGETS = TestJobOpts

#--- if using 'make all' ---#
all : $(TARGETS)

#--- if making individual targets ---#

TestJobOpts.o : TestJobOpts.cxx ../JobOptsSvc.o
	$(CXX) $(INCLUDE) $(CXXFLAGS) $(ROOTFLAGS) -o $*.o $(LDLIBS) -c $*.cxx #compile
	$(CXX) $(INCLUDE) $(CXXFLAGS) $(ROOTFLAGS) $(LDLIBS) -o $* $*.o        #link

clean:
	rm -f $(BINARIES) $(TARGETS)

