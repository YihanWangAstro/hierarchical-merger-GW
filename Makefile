# Three builds of the same pipeline (include/bhs/), differing only in the integrator:
#   bin/scatter_spacehub  SpaceHub AR-chain (needs the SpaceHub sources): the one to use
#   bin/scatter_bs        Bulirsch-Stoer extrapolation   } for comparison
#   bin/scatter_leapfrog  time-transformed leapfrog      }

CXX      ?= c++
CXXFLAGS ?= -O3 -std=c++17 -pthread -Wall -Wextra
# On a cluster: make CXX=g++ LDFLAGS="-static-libstdc++ -static-libgcc" keeps the programs
# from needing the compiler's module when they run.
LDFLAGS  ?=
SPACEHUB ?= $(HOME)/Repositories/SpaceHub/src

HEADERS := $(wildcard include/bhs/*.hpp)

all: bin/scatter_bs bin/scatter_leapfrog bin/scatter_spacehub

bin/scatter_%: src/scatter_%.cpp $(HEADERS)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -Iinclude $< -o $@ $(LDFLAGS)

# SpaceHub is header-only; -isystem keeps its warnings out of the build log.
bin/scatter_spacehub: src/scatter_spacehub.cpp $(HEADERS)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -Iinclude -isystem $(SPACEHUB) $< -o $@ $(LDFLAGS)

clean:
	rm -rf bin

.PHONY: all clean
