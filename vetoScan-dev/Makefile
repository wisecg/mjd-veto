# Generic MJD makefile (MGSW+ROOT)
# 60% of the time, it works all the time.
# C. Wiseman

# executable name
PROGRAM = vetoScan

# this is a list of the .cc files
SOURCES = $(PROGRAM).cc $(wildcard ./code/*.cc)

# this is a list of .o files created by the makefile
OBJECTS = $(SOURCES:.cc=.o)

# ROOT5 (osx and linux)
# CXX = g++
# CXXFLAGS = -c -Wall -fPIC -O2 -g
# LD = g++
# LDFLAGS = -O  $(shell root-config --ldflags)

#ROOT6 (osx: clang, linux: g++)
UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
    CXX = g++
    CXXFLAGS = -O -Wall -c -fPIC $(shell root-config --cflags)
    LD = g++
    LDFLAGS  = -O -m64
endif
ifeq ($(UNAME), Darwin)
    CXX = clang++
    CXXFLAGS = -O -Wall -c -fPIC $(shell root-config --cflags)
    LD = clang++
    LDFLAGS  = -O  $(shell root-config --ldflags)
endif

MGDOINCLUDE = -I$(MGDODIR)/Base -I$(MGDODIR)/Root -I$(MGDODIR)/Transforms -I$(MGDODIR)/Majorana -I$(MGDODIR)/Tabree -I$(MGDODIR)/MJDB
MGDOLIB = -L$(MGDODIR)/lib -lMGDOBase -lMGDOMajorana -lMGDORoot  -lMGDOMJDB #-lMGDOTransforms #-lMGDOTabree
TAMINCLUDE = -I$(MGDODIR)/tam/inc
TAMLIB = -L$(MGDODIR)/tam/lib -lTAM
GATINCLUDE= -I$(GATDIR)/BaseClasses -I$(GATDIR)/DCProcs -I$(GATDIR)/MGOutputMCRunProcessing -I$(GATDIR)/MGTEventProcessing -I$(GATDIR)/buildtools
GATLIB = -L$(GATDIR)/lib -lGATBaseClasses -lGATDCProcs -lGATMGOutputMCRunProcessing -lGATMGTEventProcessing #-lGATbuildtools
CLHEPINCLUDE = -I$(CLHEP_INCLUDE_DIR)
ROOTLIB= $(shell root-config --libs)
ROOTINCLUDE = -I$(ROOTSYS)/include
ALLINC= -I. $(ROOTINCLUDE) $(MGDOINCLUDE) $(GATINCLUDE) $(CLHEPINCLUDE) $(TAMINCLUDE)
ALLLIB= $(ROOTLIB) $(MGDOLIB) $(GATLIB) $(TAMLIB)

#####################

all: $(SOURCES) $(PROGRAM)
	@echo sources: $(SOURCES)  executable is: $(PROGRAM)

# $< is the first term in the dependencies list
# -o creates an object file
# $@ refers to the target (left side of the : symbol)

.cc.o:
#	@echo creating object file.
#	@echo target list: $(OBJECTS)  source list: $(SOURCES)
	$(CXX) $(CXXFLAGS) $(ALLINC) $< -o $@

$(PROGRAM): $(OBJECTS)
	@echo creating executable ...
	$(LD) $(LDFLAGS) $(ALLLIB) $(OBJECTS) -o $@

.PHONY: clean

clean:
	find . -name "*.o" -type f -delete
	rm $(PROGRAM)