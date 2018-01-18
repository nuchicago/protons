BINDIR := bin
LIBDIR := lib

override TOPDIR := $(PWD)/

LIB := libanalysisutils.a

export BINDIR LIBDIR TOPDIR

LIBSRCS := $(wildcard src/Selection/*.cpp) $(wildcard src/Utilities/*.cpp) 
LIBOBJS := $(foreach v,$(notdir $(LIBSRCS:.cpp=.o)),$(LIBDIR)/$v)

vpath %.cpp src/Selection src/Utilities

CXX := g++
CPPFLAGS := $(foreach path,$(shell root-config --incdir),-I$(path))
CXXFLAGS := $(filter-out -I%,$(shell root-config --cflags)) -I$(BHEP)/include -I$(ROOTSYS)/include \
            -I$(GSL_INC) -Isrc/Selection -Isrc/Utilities 

all: bin

$(LIBDIR)/$(LIB): $(LIBOBJS)

%.a:
	$(AR) $(ARFLAGS) $(@) $(^) 

$(LIBDIR)/%.o : %.cpp 
	$(CXX) -g -Wall -std=c++11 $(CPPFLAGS) -c $(CXXFLAGS) -o $(@) $(<)

echo_%:
	@echo "$(subst echo_,,$@)=\"$($(subst echo_,,$@))\""
	@echo "origin $(subst echo_,,$@) returns $(origin $(subst echo_,,$@))"


BIN_SUBDIRS = src/ProtonAnalyzerMC \
	      src/ProtonXsec	


lib: $(LIBDIR)/$(LIB)

bin: $(foreach dir,$(BIN_SUBDIRS),$(dir)/Makefile) lib
	$(foreach dir,$(BIN_SUBDIRS),$(MAKE) -C $(dir) bin;)

clean:
	-$(RM) $(LIBDIR)/$(LIB) $(LIBOBJS)
	-$(RM) $(BINDIR)/*.exe
	-$(RM) src/Selection/*.o
	-$(RM) src/Utilities/*.o
	$(foreach dir,$(BIN_SUBDIRS),$(RM) $(dir)/*.o;)

