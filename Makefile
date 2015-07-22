
# Makefile for Branching Process Simulation Project
# homebrew c++ compiler (in /usr/local/Cellar/gcc/4.9.1)
CC=g++-4.9

# compiler flags
CFLAGS=-W -Wno-long-long -pedantic -Wno-variadic-macros -std=c++11 -O3

INCLUDE=

all: branching

branching: cauloprocess.cpp branching.h
	$(CC) $(CFLAGS) $(INCLUDE) cauloprocess.cpp -o branching

clean: 
	rm -rf *.o branching *~
