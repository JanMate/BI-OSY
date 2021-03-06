CC        = gcc
CFLAGS    = -Wall -pedantic -ansi 
LIBS      = -pthread  
CXX       = g++
CXXFLAGS  = -Wall -pedantic

TARGETS   = sample.cpp test.inc

# -fsanitize=address

all: $(TARGETS)

sample: sample.cpp test.inc
	$(CXX) -std=c++11 $(CXXFLAGS) -g -o $@ $< $(LIBS)

clean:
	\rm -f $(TARGETS) *~ core

