override BINDIR := $(TOPDIR)$(BINDIR)

override LIBDIR := $(TOPDIR)$(LIBDIR)

PROG = $(BINDIR)/ProtonXsec.exe

LIBS  = -L$(ROOTSYS)/lib -lCore -lHist -lGraf -lGraf3d -lMatrix -lTree -lGpad -lRIO -lPostscript -ldl \
	-L$(LIBDIR) -lanalysisutils 

INCLUDES = -I$(ROOTSYS)/include


SRCS = ProtonXsec.cpp 
OBJS = ProtonXsec.o 


all : bin

bin: $(PROG)

clean :
	rm -f $(OBJS) $(PROG)

$(PROG) : $(OBJS)
	$(CXX) -g -Wall -std=c++11 $(OBJS) $(LIBS) -o $(PROG)

%.o : %.cpp
	$(CXX) -g -Wall -std=c++1y -c $(INCLUDES) $<
	#$(CXX) -g -Wall -std=c++11 -c $(INCLUDES) $<
