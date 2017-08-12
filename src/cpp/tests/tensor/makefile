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

all: test_tensor
#Terminate after N errors
ERRORFLG=-fmax-errors=5

test_tensor: test_tensor.o tensor.o
	$(CC) $(STD) -o $@ test_tensor.o tensor.o $(CFLAGS) $(EIGEN) $(ERRORFLG) $(DBG)

test_tensor.o: ../../tensor.h test_tensor.cpp
	$(CC) $(STD) -o $@ -c test_tensor.cpp $(CFLAGS) $(EIGEN) $(ERRORFLG) $(DBG)

tensor.o: ../../tensor.h ../../tensor.cpp
	$(CC) $(STD) -o $@ -c ../../tensor.cpp $(CFLAGS) $(EIGEN) $(ERRORFLG) $(DBG)

clean:
	rm *o test_tensor