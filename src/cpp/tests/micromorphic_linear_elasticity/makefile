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

all: test_micromorphic_linear_elasticity
#Terminate after N errors
ERRORFLG=-fmax-errors=5

test_micromorphic_linear_elasticity: test_micromorphic_linear_elasticity.o micro_material.o tensor.o
	$(CC) $(STD) -o $@ test_micromorphic_linear_elasticity.o micro_material.o tensor.o $(CFLAGS) $(EIGEN) $(ERRORFLG) $(DBG)

test_micromorphic_linear_elasticity.o: test_micromorphic_linear_elasticity.cpp ../../tensor.h ../../finite_difference.h ../../micromorphic_linear_elasticity.h
	$(CC) $(STD) -o $@ -c test_micromorphic_linear_elasticity.cpp $(CFLAGS) $(EIGEN) $(ERRORFLG) $(DBG)

micro_material.o: ../../tensor.h ../../micromorphic_linear_elasticity.h ../../micromorphic_linear_elasticity.cpp
	$(CC) $(STD) -o $@ -c ../../micromorphic_linear_elasticity.cpp $(CFLAGS) $(EIGEN) $(ERRORFLG) $(DBG)

tensor.o: ../../tensor.h ../../tensor.cpp
	$(CC) $(STD) -o $@ -c ../../tensor.cpp $(CFLAGS) $(EIGEN) $(ERRORFLG) $(DBG)

clean:
	rm *o test_micromorphic_linear_elasticity