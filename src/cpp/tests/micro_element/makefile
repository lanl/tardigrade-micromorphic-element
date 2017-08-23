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

all: test_micro_element
#Terminate after N errors
ERRORFLG=-fmax-errors=5

test_micro_element: test_micro_element.o micro_element.o tensor.o micro_material.o
	$(CC) $(STD) -o $@ test_micro_element.o micro_element.o micro_material.o tensor.o $(CFLAGS) $(EIGEN) $(ERRORFLG) $(DBG)

test_micro_element.o: test_micro_element.cpp ../../micro_element.h ../../tensor.h ../../finite_difference.h
	$(CC) $(STD) -o $@ -c test_micro_element.cpp $(CFLAGS) $(EIGEN) $(ERRORFLG) $(DBG)

micro_element.o: ../../micro_element.h ../../tensor.h ../../micro_element.cpp ../../micromorphic_linear_elasticity.h
	$(CC) $(STD) -o $@ -c ../../micro_element.cpp $(CFLAGS) $(EIGEN) $(ERRORFLG) $(DBG)

micro_material.o: ../../tensor.h ../../micromorphic_linear_elasticity.h ../../micromorphic_linear_elasticity.cpp
	$(CC) $(STD) -o $@ -c ../../micromorphic_linear_elasticity.cpp $(CFLAGS) $(EIGEN) $(ERRORFLG) $(DBG)

tensor.o: ../../tensor.h ../../tensor.cpp
	$(CC) $(STD) -o $@ -c ../../tensor.cpp $(CFLAGS) $(EIGEN) $(ERRORFLG) $(DBG)

clean:
	rm *o test_micro_element