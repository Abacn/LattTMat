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
            COMPFLG_OMP = -fopenmp
            LINKERFLAGS_OMP = -lgomp 
        else
            CC = icpc
            COMPFLG = -I$(INCLUDE)
            LINKERFLAGS = 
            COMPFLG_OMP = -qopenmp
            LINKERFLAGS_OMP = -qopenmp
        endif
    endif
    ifeq ($(UNAME_S),Darwin)
        CC = clang++
        COMPFLG = 
        LINKERFLAGS = 
        COMPFLG_OMP = -Xpreprocessor -fopenmp
        LINKERFLAGS_OMP = -lpthread -lomp
    endif
endif


CPPFLAGS = -std=c++11 -march=haswell -O3 -Wall -Wno-unused-variable -Wno-overloaded-virtual $(COMPFLG)

SRCDIR = TMat
BUILDDIR = build/single/$(SRCDIR)
BINDIR = ../bin
BINNAME = $(BINDIR)/tmat
BUILDDIR_OMP = build/omp/$(SRCDIR)
BINNAME_OMP = $(BINDIR)/tmat_omp

SOURCES = utility.cpp tmat.cpp tmatreduce.cpp tmat3d.cpp tmatphys.cpp pfunc.cpp read_input.cpp tests.cpp main.cpp
OBJECTS := $(addprefix $(BUILDDIR)/,$(SOURCES:.cpp=.o))
DEPENDS = $(OBJECTS:%.o=%.d)
OBJECTS_OMP := $(addprefix $(BUILDDIR_OMP)/,$(SOURCES:.cpp=.o))
DEPENDS_OMP = $(OBJECTS_OMP:%.o=%.d)

.PHONY: clean

all: $(BINNAME) $(BINNAME_OMP)

$(BINNAME): $(OBJECTS) | $(BINDIR)
	$(CC) $(CPPFLAGS) -o $(BINNAME) $(OBJECTS) $(LINKERFLAGS)

$(BINNAME_OMP): $(OBJECTS_OMP) | $(BINDIR)
	$(CC) $(CPPFLAGS) -o $(BINNAME_OMP) $(OBJECTS_OMP) $(LINKERFLAGS) $(LINKERFLAGS_OMP)

-include $(DEPENDS)

$(OBJECTS): $(BUILDDIR)/%.o: $(SRCDIR)/%.cpp | $(BUILDDIR)
	$(CC) $(CPPFLAGS) -MMD -c $< -o $@

$(OBJECTS_OMP): $(BUILDDIR_OMP)/%.o: $(SRCDIR)/%.cpp | $(BUILDDIR_OMP)
	$(CC) $(CPPFLAGS) $(COMPFLG_OMP) -MMD -c $< -o $@

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(BUILDDIR_OMP):
	mkdir -p $(BUILDDIR_OMP)

$(BINDIR):
	mkdir -p $(BINDIR)

clean:
	rm -rf $(BUILDDIR) $(BUILDDIR_OMP) $(BINNAME) $(BINNAME_OMP) $(DEPENDS)
