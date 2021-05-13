ifeq ($(OS),Windows_NT)
    $(error Operation system not support)
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux) # cluster in school
        ifeq (, $(shell which icpc 2>/dev/null))
            CC = g++
            COMPFLG = -I$(INCLUDE)
            # condor
            ifneq (, $(shell which condor_submit 2>/dev/null))
                LINKERFLAGS = -static
            else
                LINKERFLAGS = 
            endif
        else
            CC = icpc
            COMPFLG = -I$(INCLUDE)
            LINKERFLAGS = 
        endif
    endif
    ifeq ($(UNAME_S),Darwin)
        CC = clang++
        COMPFLG = 
        LINKERFLAGS = 
    endif
endif


CPPFLAGS = -std=c++11 -march=haswell -O3 -Wall -Wno-unused-variable -Wno-overloaded-virtual $(COMPFLG)
LINKERFLAGS = 

SRCDIR = DumpConfig
SRCDIR_TMAT = TMat
BUILDDIR = build/$(SRCDIR)
BUILDDIR_TMAT = build/$(SRCDIR_TMAT)
BINDIR = ../bin
BINNAME = $(BINDIR)/dumpconfig

SOURCES = utilityb.cpp dumpconfig.cpp dumpreduce.cpp main.cpp
SOURCES_TMAT = read_input.cpp tmatreduce.cpp tmat.cpp tmatphys.cpp utility.cpp

OBJECTS := $(addprefix $(BUILDDIR)/,$(SOURCES:.cpp=.o)) 
OBJECTS_TMAT = $(addprefix $(BUILDDIR_TMAT)/,$(SOURCES_TMAT:.cpp=.o))
DEPENDS = $(OBJECTS:%.o=%.d)

.PHONY: clean

all: $(BINNAME)

$(BINNAME): $(OBJECTS) $(OBJECTS_TMAT) | $(BINDIR)
	$(CC) $(CPPFLAGS) -o $(BINNAME) $(OBJECTS) $(OBJECTS_TMAT) $(LINKERFLAGS)

-include $(DEPENDS)

$(OBJECTS): $(BUILDDIR)/%.o: $(SRCDIR)/%.cpp | $(BUILDDIR)
	$(CC) $(CPPFLAGS) $(COMPFLAGS) -MMD -c $< -o $@

$(OBJECTS_TMAT): $(BUILDDIR_TMAT)/%.o: $(SRCDIR_TMAT)/%.cpp | $(BUILDDIR_TMAT)
	$(CC) $(CPPFLAGS) $(COMPFLAGS) -MMD -c $< -o $@

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(BUILDDIR_TMAT):
	mkdir -p $(BUILDDIR_TMAT)
	
$(BINDIR):
	mkdir -p $(BINDIR)

clean:
	rm -rf $(BUILDDIR) $(BINNAME) $(DEPENDS)
