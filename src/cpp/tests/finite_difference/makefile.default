#Compiler option
CC=COMPILER_COMMAND

#Standard option
STD=-std=gnu++11

#Compiler flags
CFLAGS=-I. -I ../..

#Include Eigen Library
EIGEN = -I EIGEN_LOCATION

#Debugging flag
DBG = -ggdb

all: test_finite_difference
#Terminate after N errors
ERRORFLG=-fmax-errors=5

test_finite_difference: test_finite_difference.cpp ../../finite_difference.h
	$(CC) $(STD) -o $@ test_finite_difference.cpp $(CFLAGS) $(EIGEN) $(ERRORFLG) $(DBG)
clean:
	rm *o test_finite_difference
