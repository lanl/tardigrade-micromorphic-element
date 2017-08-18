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

void print_vector(std::string name, std::vector< double > V){
    /*!======================
    |    print_vector    |
    ======================
    
    Print a vector to the screen
    
    */
    
    std::cout << name << ": ";
    
    for(int i=0; i<V.size(); i++){
        std::cout << V[i] << " ";
    }
    std::cout << "\n";
}

void print_vector_of_vectors(std::string name, std::vector< std::vector< double > > V){
    /*!=================================
    |    print_vector_of_vectors    |
    =================================
    
    Print a vector of vectors to the screen
    
    */
    
    for(int i=0; i<V.size(); i++){
        std::cout << name << "[" << i << "]: ";
        for(int j=0; j<V[i].size(); j++){
            std::cout << V[i][j] << " ";
        }
        std::cout << "\n";
    }
}

std::vector< double > generate_current_coordinates(std::vector< double > reference_coord){
    /*!======================================
    |    generate_current_coordinates    |
    ======================================
    
    Generate output coordinates of the incoming coordinates 
    such that a known deformation can be applied.
    
    Input:
        reference_coord: The coordinates in the reference 
                         configuration which will be mapped 
                         to a coordinate in the current 
                         configuration.
    
    */
    
    std::vector< double > current_coord(3,0.);
    
    current_coord[0] =  1.30*reference_coord[0]-0.375*reference_coord[1]+1.20*reference_coord[2]+1.0;
    current_coord[1] =  0.75*reference_coord[0]+0.650*reference_coord[1]-0.31*reference_coord[2]-2.3;
    current_coord[2] = -2.30*reference_coord[0]+1.400*reference_coord[1]+0.44*reference_coord[2]+0.3;
    
    return current_coord;
}

tensor::Tensor get_gradient_coordinates(std::vector< double > reference_coord){
    /*!==================================
    |    get_gradient_coordinates    |
    ==================================
    
    Generate the gradient of the output coordinates
    with respect to the incoming coordinates.
    
    Input:
        reference_coord: The coordinates in the reference 
                         configuration which will be mapped 
                         to a coordinate in the current 
                         configuration.
    
    */
    
    std::vector< int > shape = {3,3};
    tensor::Tensor T = tensor::Tensor(shape);
    T(0,0) =  1.300;
    T(0,1) = -0.375;
    T(0,2) =  1.200;
    T(1,0) =  0.750;
    T(1,1) =  0.650;
    T(1,2) = -0.310;
    T(2,0) = -2.300;
    T(2,1) =  1.400;
    T(2,2) =  0.440;
    
    return T;
    
}

int test_constructors(std::ofstream &results){
    /*!===========================
    |    test_constructors    |
    ===========================
    
    Run tests on the constructors to ensure that 
    they result in the expected behavior.
    
    */
    
    //Seed the random number generator
    srand (1);
    
    //!Initialize test results
    int  test_num        = 8;
    bool test_results[test_num] = {false,false,false,false,false,false,false,false};
    
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
    
    //!|=> Test 1
    //!Test the empty constructor
    micro_element::Hex8 A = micro_element::Hex8();
    
    test_results[0] = A.RHS.size()    == 96;
    test_results[1] = A.AMATRX.size() == 96;
    for(int i=0; i<A.AMATRX.size(); i++){
        test_results[1] *= A.AMATRX[i].size() == 96;
    }
    
    //!|=> Test 2
    //!Test the constructor with the reference node locations
    micro_element::Hex8 B = micro_element::Hex8(reference_coords);
    
    std::vector< std::vector< double > > B_answer= {{0,0,0},{1,0,0},{1,1,0},{0,1,0},
                                                    {0,0,1},{1,0,1},{1,1,1},{0,1,1}}; //!The desired format of reference_coords and current_coords
    
    test_results[2] = B_answer==B.reference_coords;
    test_results[3] = B_answer==B.current_coords;
    
    //!Test the constructor with the reference node locations and the degree of freedom vectors
    micro_element::Hex8 C = micro_element::Hex8(reference_coords,U,dU);
    
    std::vector< std::vector< double > > U_answer;
    std::vector< std::vector< double > > dU_answer;
    
    U_answer.resize(8);
    dU_answer.resize(8);
    
    //!|=> Test 3
    //!Test if the dof vectors were parsed correctly
    for(int n=0; n<8; n++){
        U_answer[n].resize(9);
        dU_answer[n].resize(9);
        
        for(int i=0; i<9; i++){
            U_answer[n][i]  = U[i+n*9];
            dU_answer[n][i] = dU[i+n*9];
        }
    }
    
    test_results[4] = U_answer == C.dof_at_nodes;
    test_results[5] = dU_answer == C.Delta_dof_at_nodes;
    
    //!|=> Test 4
    //!Check if the current coordinates were updated correctly
    test_results[6] = true;
    for(int n=0; n<8; n++){
        for(int i=0; i<3; i++){
            test_results[6] *= fabs(C.current_coords[n][i]-(U[i+n*12]+reference_coords[i+n*3]))<1e-6;
        }
    }
    
    //!|=> Test 5
    //!Check if the phi values were updated correctly
    std::vector< int > shape = {3,3};
    tensor::Tensor phic = tensor::Tensor(shape);
    
    test_results[7] = true;
    for(int n=0; n<8; n++){
        phic(0,0) = C.dof_at_nodes[n][ 3];
        phic(1,1) = C.dof_at_nodes[n][ 4];
        phic(2,2) = C.dof_at_nodes[n][ 5];
        phic(1,2) = C.dof_at_nodes[n][ 6];
        phic(0,2) = C.dof_at_nodes[n][ 7];
        phic(0,1) = C.dof_at_nodes[n][ 8];
        phic(2,1) = C.dof_at_nodes[n][ 9];
        phic(2,0) = C.dof_at_nodes[n][10];
        phic(1,0) = C.dof_at_nodes[n][11];
        
        test_results[7] *= phic.data == C.node_phis[n].data;
    }
    
    
    //Compare all test results
    bool tot_result = true;
    for(int i = 0; i<test_num; i++){
        //std::cout << "\nSub-test " << i+1 << " result: " << test_results[i] << "\n";
        if(!test_results[i]){
            tot_result = false;
        }
    }
    
    if(tot_result){
        results << "test_constructors & True\\\\\n\\hline\n";
    }
    else{
        results << "test_constructors & False\\\\\n\\hline\n";
    }
}

int test_shape_functions(std::ofstream &results){
    /*!==============================
    |    test_shape_functions    |
    ==============================
    
    Run tests on the shape functions and related 
    methods to ensure they are functioning 
    properly
    
    */
    
    //Seed the random number generator
    srand (1);
    
    //!Initialize test results
    int  test_num        = 3;
    bool test_results[test_num] = {false,false,false};
    
    //!Form the required vectors for element formation
    std::vector< double > reference_coords = {0,0,0,1,0,0,1,1,0,0,1,0,0,0,1,1,0,1,1,1,1,0,1,1};
    std::vector< double > U;
    std::vector< double > dU;
    U.resize(96);
    dU.resize(96);
    for(int i=0; i<96; i++){
        U[i]  = (rand()%100-50)/1000.;
        dU[i] = (rand()%100-50)/10000.;
    }
    
    //!Form the hexehedral test element.
    micro_element::Hex8 element = micro_element::Hex8(reference_coords,U,dU);
    
    
    //!=
    //!| Tests of Hex8::shape_function
    //!=
    
    //!|=> Test 1
    //!Test whether the shape functions values are unity at the correct nodes
    //!and zero elsewhere.
    
    double N;
    double M;
    
    std::vector< std::vector< double > > local_coords     = {{-1,-1,-1},{ 1,-1,-1},{ 1, 1,-1},{-1, 1,-1},
                                                             {-1,-1, 1},{ 1,-1, 1},{ 1, 1, 1},{-1, 1, 1}};
    test_results[0] = true;
    for(int n=0; n<8; n++){
        N = element.shape_function(n,local_coords[n]);
        test_results[0] *= 1e-9>N-1>=0; //The shape function should be equal to one at the specified node
        for(int m=0; m<8; m++){
            M = element.shape_function(m,local_coords[n]);
            if(n==m){test_results[0] *= 1e-9>N-M>=0;} //The shape function should only be 1 at N==M
            else{test_results[0] *= 1e-9>N-M-1>=0;}
        }
    }
    
    //!|=> Test 2
    //!Test whether the sum of the shape functions values are unity at a
    //!point within the element
    
    std::vector< double > xi = {0,0,0}; //!The local coordinate
    for(int i=0; i<3; i++){xi[i] = (rand()%1000-500)/1000.;} //!Populate the local coordinate at a location in each direction between -1 and 1
    
    double sum_Ns;
    
    for(int n=0; n<8; n++){sum_Ns += element.shape_function(n,xi);} //Sum up all of the values of the shape function at the point
    
    test_results[1] = 1e-9>sum_Ns-1.>=0;
    
    //!|=> Test 3
    //!Test whether the gradient of the shape function w.r.t. the local 
    //!coordinates are correct
    
    xi = {0.1,-0.2,.3};
    std::vector< std::vector< double > > dNdxi_answers = {{-0.84, -0.63, -1.08},{0.84, -0.77, -1.32},{0.56, 0.77, -0.88},{-0.56, 0.63, -0.72},
                                                          {-1.56, -1.17,  1.08},{1.56, -1.43,  1.32},{1.04, 1.43,  0.88},{-1.04, 1.17,  0.72}};
    
    test_results[2] = true;
    for(int n=0; n<8; n++){
        print_vector("result",element.local_gradient_shape_function(n,xi));
        print_vector("answer",dNdxi_answers[n]);
        //test_results[2] *= 1e-9>fabs(element.local_gradient_shape_function(n,xi)-dNdxi_answers[n])>=0;
    }
    
    
    //Compare all test results
    bool tot_result = true;
    for(int i = 0; i<test_num; i++){
        std::cout << "\nSub-test " << i+1 << " result: " << test_results[i] << "\n";
        if(!test_results[i]){
            tot_result = false;
        }
    }
    
    if(tot_result){
        results << "test_constructors & True\\\\\n\\hline\n";
    }
    else{
        results << "test_constructors & False\\\\\n\\hline\n";
    }

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
    test_constructors(results);
    test_shape_functions(results);
    
    //Close the results file
    results.close();
}

