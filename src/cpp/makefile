#Compiler option
CC=g++

#Standard option
STD=-std=gnu++11

#Compiler flags
CFLAGS=-I.

#Include Eigen Library
EIGEN = -I C:/usr/share/cpp/eigen3.3.4

#Debugging flag
DBG = -ggdb

all: driver
#Terminate after N errors
ERRORFLG=-fmax-errors=5

driver: driver.o micro_element.o tensor.o micro_material.o newton_krylov.o
	$(CC) $(STD) -o $@ driver.o micro_element.o micro_material.o tensor.o newton_krylov.o $(CFLAGS) $(EIGEN) $(ERRORFLG) $(DBG)

driver.o: driver.cpp micro_element.h tensor.h
	$(CC) $(STD) -o $@ -c driver.cpp $(CFLAGS) $(EIGEN) $(ERRORFLG) $(DBG)

micro_element.o: micro_element.h tensor.h micro_element.cpp micromorphic_linear_elasticity.h
	$(CC) $(STD) -o $@ -c micro_element.cpp $(CFLAGS) $(EIGEN) $(ERRORFLG) $(DBG)

micro_material.o: tensor.h micromorphic_linear_elasticity.h micromorphic_linear_elasticity.cpp
	$(CC) $(STD) -o $@ -c micromorphic_linear_elasticity.cpp $(CFLAGS) $(EIGEN) $(ERRORFLG) $(DBG)

tensor.o: tensor.h tensor.cpp
	$(CC) $(STD) -o $@ -c tensor.cpp $(CFLAGS) $(EIGEN) $(ERRORFLG) $(DBG)

newton_krylov.o: newton_krylov.h newton_krylov.cpp
	$(CC) $(STD) -o $@ -c newton_krylov.cpp $(CFLAGS) $(EIGEN) $(ERRORFLG) $(DBG)

clean:
	rm *o test_micro_element