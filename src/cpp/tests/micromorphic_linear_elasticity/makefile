# Micromorphic linear elasticity library
# 
# Author: Nathan A. Miller (LANL / CU Boulder)
# Email:  nathanm@lanl.gov
# Date:   February 21, 2020
#
# This is the makefile for test_micromorphic_linear_elasticity.cpp it tests the 
# micromorphic linear elastic constitutive model.

include ../../../../config.mk

# List of executables
EXECUTABLE=test_micromorphic_linear_elasticity
TESTOBJECT=test_micromorphic_linear_elasticity.o
TESTSOURCE=test_micromorphic_linear_elasticity.cpp

OBJECT=micromorphic_linear_elasticity.o
SOURCE=../../micromorphic_linear_elasticity.cpp
HEADER=../../micromorphic_linear_elasticity.h

MICROMORPHICTOOLSOBJECT=micromorphic_tools.o
CONSTITUTIVEOBJECT=constitutive_tools.o
ERROROBJECT=error_tools.o
MICROMORPHICMATERIALLIBOBJECT=micromorphic_material_library.o

OBJECTS = $(TESTOBJECT)
OBJECTS += $(OBJECT) $(ERROROBJECT) $(CONSTITUTIVEOBJECT) $(MICROMORPHICTOOLSOBJECT)
OBJECTS += $(MICROMORPHICMATERIALLIBOBJECT)

INC+=-I../..
LOCAL_CFLAGS=-ggdb

#Set the mode due to the method variable
DBG =
ifeq ($(METHOD),dbg)
  DBG = -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -g
endif
ifeq ($(METHOD),oprof)
  DBG = -pg
endif

# Makefile rules
all: $(EXECUTABLE)

$(EXECUTABLE): $(TESTOBJECT) $(OBJECT) $(ERROROBJECT) $(CONSTITUTIVEOBJECT) $(MICROMORPHICTOOLSOBJECT) $(MICROMORPHICMATERIALLIBOBJECT)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) $^ -o $@ $(CFLAGS) -lstdc++ -std=c++11

$(TESTOBJECT): $(TESTSOURCE) $(HEADER) $(VECTORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) $< -c -o $@ $(CFLAGS) -lstdc++ -std=c++11

$(OBJECT): $(SOURCE) $(HEADER) $(MICROMORPHICTOOLSHEADER) $(VECTORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=c++11

$(MICROMORPHICTOOLSOBJECT): $(MICROMORPHICTOOLSSOURCE) $(MICROMORPHICTOOLSHEADER) $(CONSITUTIVEHEADER) $(VECTORHEADER) $(ERRORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=c++11

$(ERROROBJECT): $(ERRORSOURCE) $(ERRORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=c++11

$(CONSTITUTIVEOBJECT): $(CONSTITUTIVESOURCE) $(CONSTITUTIVEHEADER) $(VECTORHEADER) $(ERRORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=c++11

$(MICROMORPHICMATERIALLIBOBJECT): $(MICROMORPHICMATERIALLIBSOURCE) $(MICROMORPHICMATERIALLIBHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=c++11

clean:
	rm -f $(EXECUTABLE)
	rm $(OBJECTS)

.PHONY: all clean
