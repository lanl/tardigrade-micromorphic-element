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
        U_answer[n].resize(12);
        dU_answer[n].resize(12);
        
        for(int i=0; i<12; i++){
            U_answer[n][i]  = U[i+n*12];
            dU_answer[n][i] = dU[i+n*12];
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
    
    return 1;
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
    int  test_num        = 6;
    bool test_results[test_num] = {false,false,false,false,false,false};
    
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
    //!| Tests of Hex8::set_shape_function and Hex8::set_shape_functions
    //!=
    
    //!|=> Test 1
    //!Test whether the shape functions values are unity at the correct nodes
    //!and zero elsewhere.
    
    element.points     = {{-1,-1,-1},{ 1,-1,-1},{ 1, 1,-1},{-1, 1,-1},
                          {-1,-1, 1},{ 1,-1, 1},{ 1, 1, 1},{-1, 1, 1}}; //Set the gauss points of the element to the nodes (done for testing purposes)
    test_results[0] = true;
    for(int n=0; n<8; n++){
        element.set_gpt_num(n);        //Set the, "gauss point," of the element to the current node
        element.set_shape_functions(); //Compute all of the values of the shape functions
        test_results[0] *= (1e-9>(element.get_N(n)-1) && (element.get_N(n)-1))>=0; //The shape function should be equal to one at the specified node
        for(int m=0; m<8; m++){
            if(n==m){test_results[0] *= ((1e-9>(element.get_N(n)-element.get_N(m))) && ((element.get_N(n)-element.get_N(m))>=0));} //The shape function should only be 1 at N==M
            else{test_results[0] *= 1e-9>element.get_N(n)-element.get_N(m)-1>=0;}
        }
    }
    
    
    
    //!|=> Test 2
    //!Test whether the sum of the shape functions values are unity at a
    //!point within the element
    
    element.points[0] = {0,0,0}; //!The local coordinate (set as a false gauss point)
    for(int i=0; i<3; i++){element.points[0][i] = (rand()%1000-500)/1000.;} //!Populate the local coordinate at a location in each direction between -1 and 1
    element.set_gpt_num(0);

    double sum_Ns;
    
    element.set_shape_functions(); //Compute all of the nodal shape functions at the given point

    for(int n=0; n<8; n++){sum_Ns += element.get_N(n);} //Sum up all of the values of the shape function at the point
    
    test_results[1] = ((1e-9>(sum_Ns-1.)) && ((sum_Ns-1.)>=0));
    
    //!|=> Test 3
    //!Test whether the gradient of the shape function w.r.t. the local 
    //!coordinates are correct at the center of the element
    
    element.points[0] = {0.,0.,0.};
    std::vector< std::vector< double > > dNdxi_answers = {{ -0.125, -0.125, -0.125},{  0.125, -0.125, -0.125},{  0.125,  0.125, -0.125},{-0.125,  0.125, -0.125},
                                                          { -0.125, -0.125,  0.125},{  0.125, -0.125,  0.125},{  0.125,  0.125,  0.125},{-0.125,  0.125,  0.125}};
    
    element.set_gpt_num(0);
    std::cout << "Setting the local gradient of the shape functions\n";
    element.set_local_gradient_shape_functions();
    std::cout << "Comparing the results vs. the answers\n";
  
    test_results[2] = true;
    double temp_diff;
    for(int n=0; n<8; n++){
        for(int i=0; i<3; i++){
            temp_diff        = fabs(element.get_dNdxi(n)[i]-dNdxi_answers[n][i]);
            test_results[2] *= (1e-9>temp_diff);
        }
    }
    
    //!|=> Test 4
    //!Test whether the gradient of the shape function w.r.t. the local
    //!coordinates are correct at a location off the center of the 
    //!element.
    
    element.points[0] = {0.3,-0.6,0.7};
    element.set_gpt_num(0);
    dNdxi_answers = {{-0.06   , -0.02625, -0.14},{ 0.06   , -0.04875, -0.26},{ 0.015  ,  0.04875, -0.065},{-0.015  ,  0.02625, -0.035},
                     {-0.34   , -0.14875,  0.14},{ 0.34   , -0.27625,  0.26},{ 0.085  ,  0.27625,  0.065},{-0.085  ,  0.14875,  0.035}};
    
    element.set_local_gradient_shape_functions();
  
    test_results[3] = true;
    bool temp;
    for(int n=0; n<8; n++){
        for(int i=0; i<3; i++){
            temp_diff = fabs(element.get_dNdxi(n)[i]-dNdxi_answers[n][i]);
            test_results[3] *= ((1e-9>temp_diff) && (temp_diff>=0));
        }
    }
    
    //!|=> Test 5
    //!Test whether the jacobian is computed correctly for the reference coordinates
    
    tensor::Tensor J_answer({3,3}); //!The answer jacobian.
    tensor::Tensor J_result = element.get_jacobian(0);
	
    for(int n=0; n<8; n++){
        J_answer += micro_element::vector_dyadic_product(dNdxi_answers[n],element.reference_coords[n]); //Compute the expected value of the jacobian
    }
    
    test_results[4] = J_result.data.isApprox(J_answer.data); //Test the results
	
    //!|=> Test 6
    //!Test whether the jacobian is computed correctly for the current coordinates
    
    J_answer = tensor::Tensor({3,3});
    J_result = element.get_jacobian(1);
    
    for(int n=0; n<8; n++){
        J_answer += micro_element::vector_dyadic_product(dNdxi_answers[n],element.current_coords[n]);
    }
    
    test_results[5] = J_result.data.isApprox(J_answer.data);
    
    //!|=> Test 7
    //!Test whether the gradient of the shape function w.r.t. the 
    //!current coordinates are correct at a location in the 
    //!element.
    
    //xi = {0.52,-.42,.73};
    
    
    
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
    
    return 1;

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

