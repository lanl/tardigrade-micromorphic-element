#Compiler option
CC=g++

#Standard option
STD=-std=gnu++11

#Compiler flags
CFLAGS=-I. -I ../..

#Include Eigen Library
EIGEN = -I C:/usr/share/cpp/eigen3.3.4

#Debugging flag
DBG = -ggdb

all: test_finite_difference
#Terminate after N errors
ERRORFLG=-fmax-errors=5

test_finite_difference: test_finite_difference.cpp ../../finite_difference.h
	$(CC) $(STD) -o $@ test_finite_difference.cpp $(CFLAGS) $(EIGEN) $(ERRORFLG) $(DBG)
clean:
	rm *o test_finite_difference