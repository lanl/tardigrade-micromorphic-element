/*!=======================================================
  |                                                     |
  |              test_micro_element.cpp                 |
  |                                                     |
  -------------------------------------------------------
  | The unit test file for micro_element.h/cpp.         |
  | This file tests the classes and functions defined   |
  | in micro_element.h/cpp.                             |
  |                                                     |
  | Generated files:                                    |
  |    Results.tex:  A LaTeX file which contains the    |
  |                  results as they will be included   |
  |                  in the generated report.           |
  =======================================================
  | Dependencies:                                       |
  | Eigen:  An implementation of various matrix         |
  |         commands. The implementation of the data    |
  |         matrix uses such a matrix.                  |
  | tensor: An implementation of a tensor which stores  |
  |         the data in a 2D matrix but allows access   |
  |         through n dimensional index notation.       |
  =======================================================*/
  
#include <functional>
#include <iostream>
#include <fstream>
#include <numeric>
#include <vector>
#include <Eigen/Dense>
#include <tensor.h>
#include <micro_element.h>
#include <ctime>
#include <stdlib.h>

int test_constructors(std::ofstream &results){
    /*!===========================
    |    test_constructors    |
    ===========================
    
    Run tests on the constructors to ensure that 
    they result in the expected behavior.
    
    */
    
    //Seed the random number generator
    srand (1);
    
    //!Form the required vectors
    std::vector< double > reference_coords = {0,0,0,1,0,0,1,1,0,0,1,0,0,0,1,1,0,1,1,1,1,0,1,1};
    std::vector< double > U;
    std::vector< double > dU;
    U.resize(96);
    dU.resize(96);
    for(int i=0; i<96; i++){
        U[i]  = (rand()%100-50)/1000.;
        dU[i] = (rand()%100-50)/10000.;
    }
    
    //!Test the empty constructor
    micro_element::Hex8 A = micro_element::Hex8();
    
    //!Test the constructor with the reference node locations
    micro_element::Hex8 B = micro_element::Hex8(reference_coords);
    
    //!Test the constructor with the reference node locations and the degree of freedom vectors
    micro_element::Hex8 C = micro_element::Hex8(reference_coords,U,dU);
}

int main(){
    /*!==========================
    |         main            |
    ===========================
    
    The main loop which runs the tests defined in the 
    accompanying functions. Each function should output
    the function name followed by & followed by True or 
    False if the test passes or fails respectively.*/
    
    std::ofstream results;
    //Open the results file
    results.open ("results.tex");
        
    //!Run the test functions
    
    
    //Close the results file
    results.close();
}

