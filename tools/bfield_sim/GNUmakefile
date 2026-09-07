#

SRCS = $(wildcard src/*.cc)
TARGETS = $(basename $(subst src,bin,$(SRCS)))

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)

all: $(TARGETS)

bin/calc_bfield: src/calc_bfield.cc include/calc_bfield.h include/special_coil.h
	$(info ===== Compiling $* =====)
	$(CXX) $(ROOTFLAGS) -I include src/calc_bfield.cc $(ROOTLIBS) -o bin/calc_bfield

bin/%: src/%.cc
	$(info ===== Compiling $* =====)
	$(CXX) $(ROOTFLAGS) -I include src/$*.cc $(ROOTLIBS) -o bin/$*

clean: 
	$(RM) bin/* 
