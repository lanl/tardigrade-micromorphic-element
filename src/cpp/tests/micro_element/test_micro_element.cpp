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
#include <finite_difference.h>
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

tensor::Tensor23 get_gradient_coordinates(std::vector< double > reference_coord){
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
    
    tensor::Tensor23 T;
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
    tensor::Tensor23 phic({3,3});
    
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
    int  test_num        = 8;
    bool test_results[test_num] = {false,false,false,false,false,false,false,false};
    
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
    //std::cout << "Setting the local gradient of the shape functions\n";
    element.set_local_gradient_shape_functions();
    //std::cout << "Comparing the results vs. the answers\n";
  
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
    //!at a location off the center of the element.
    
    tensor::Tensor23 J_answer({3,3}); //!The answer jacobian.
    tensor::Tensor23 J_result = element.get_jacobian(0);
    for(int n=0; n<8; n++){
        //J_answer += micro_element::vector_dyadic_product(dNdxi_answers[n],element.reference_coords[n]); //Compute the expected value of the jacobian
        J_answer += micro_element::vector_dyadic_product(element.reference_coords[n],dNdxi_answers[n]);
    }
    test_results[4] = J_result.data.isApprox(J_answer.data); //Test the results
    
    //!|=> Test 6
    //!Test whether the jacobian is computed correctly for the current coordinates
    
    J_answer.data.setZero();
    J_result = element.get_jacobian(1);
    
    for(int n=0; n<8; n++){
        //J_answer += micro_element::vector_dyadic_product(dNdxi_answers[n],element.current_coords[n]);
        J_answer += micro_element::vector_dyadic_product(element.current_coords[n],dNdxi_answers[n]);
    }
    
    test_results[5] = J_result.data.isApprox(J_answer.data);
    
    //!|=> Test 7
    //!Test whether the gradient of the shape function with respect to the 
    //!reference coordinates is computed correctly.
    
    std::vector< std::vector< double > > dNdX_answer; //Initialize the expected answer
    dNdX_answer.resize(8);
    tensor::Tensor23 Jtemp = element.get_jacobian(0); //Get dXdxi
    Jtemp = Jtemp.inverse();                          //Get dxidX
    
    for(int n=0; n<8; n++){
        dNdX_answer[n].resize(3); //Resize the gradient to be in three dimensions
        for(int i=0; i<3; i++){dNdX_answer[n][i]=0;} //Zero out the gradient
        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                dNdX_answer[n][i] += Jtemp(j,i)*dNdxi_answers[n][j]; // Compute dxi_j dX_i dN dxi_j 
            }
        }
    }
    
    element.set_global_gradient_shape_functions(0);
    
    //Check the answer to the result
    test_results[6] = true;
    for(int n=0; n<8; n++){
        //print_vector("result",element.get_dNdx(0,n));
        //print_vector("answer",dNdX_answer[n]);
        for(int i=0; i<3; i++){
            test_results[6] *= 1e-9>fabs(element.get_dNdx(0,n)[i]-dNdX_answer[n][i]);
        }
    }
    
    //!|=> Test 8
    //!Test whether the gradient of the shape function with respect to the 
    //!current coordinates is computed correctly.
    
    std::vector< std::vector< double > > dNdx_answer; //Initialize the expected answer
    dNdx_answer.resize(8);
    Jtemp = element.get_jacobian(1); //Get dxdxi
    Jtemp = Jtemp.inverse();         //Get dxidx
    
    for(int n=0; n<8; n++){
        dNdx_answer[n].resize(3); //Resize the gradient to be in three dimensions
        for(int i=0; i<3; i++){dNdx_answer[n][i]=0;} //Zero out the gradient
        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                dNdx_answer[n][i] += Jtemp(j,i)*dNdxi_answers[n][j]; // Compute dxi_j dx_i dN dxi_j 
            }
        }
    }
    
    element.set_global_gradient_shape_functions(1);
    
    //Check the answer to the result
    test_results[7] = true;
    for(int n=0; n<8; n++){
        //print_vector("result",element.get_dNdx(1,n));
        //print_vector("answer",dNdx_answer[n]);
        for(int i=0; i<3; i++){
            test_results[7] *= 1e-9>fabs(element.get_dNdx(1,n)[i]-dNdx_answer[n][i]);
        }
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
        results << "test_shape_functions & True\\\\\n\\hline\n";
    }
    else{
        results << "test_shape_functions & False\\\\\n\\hline\n";
    }
    
    return 1;

}

std::vector< double > test_deformation(std::vector< double > reference_position){
    /*!==========================
    |    test_deformation    |
    ==========================
    
    Compute a test deformation 
    to show that the deformation gradient 
    and microdisplacement are being computed 
    correctly.
    
    */
    
    double X = reference_position[0];
    double Y = reference_position[1];
    double Z = reference_position[2];
    
    std::vector< double > U;
    U.resize(12,0);
    
    //!Set the displacements
    U[ 0] =  0.32*X-0.14*Y+0.61*Z;
    U[ 1] = -0.50*X+0.24*Y-0.38*Z;
    U[ 2] = -0.22*X+0.47*Y+0.62*Z;
    U[ 3] = -1.10*X+0.04*Y+2.30*Z; //phi_11
    U[ 4] = -0.74*X+1.22*Y+2.22*Z; //phi_22
    U[ 5] = -2.24*X+5.51*Y+1.11*Z; //phi_33
    U[ 6] = -5.75*X+2.26*Y+7.66*Z; //phi_23
    U[ 7] = -6.22*X+8.63*Y+2.72*Z; //phi_13
    U[ 8] = -2.76*X+3.37*Y+3.93*Z; //phi_12
    U[ 9] = -6.32*X+6.73*Y+7.22*Z; //phi_32
    U[10] = -3.83*X+4.29*Y+1.51*Z; //phi_31
    U[11] = -9.18*X+3.61*Y+9.08*Z; //phi_21
    
    return U;
    
}
    
std::vector< std::vector< double > > compute_gradients(std::vector< double > reference_position){
    /*!===========================
    |    compute_gradients    |
    ===========================
    
    Compute the gradients of the test 
    deformation numerically so that the 
    results are consistent.
    
    */
    
    finite_difference::FiniteDifference FD(test_deformation,2,reference_position,1e-6);
    std::vector< std::vector< double > > gradient = FD.numeric_gradient();
    return gradient;
    
}

std::vector< double > parse_F(std::vector< double > U_in){
    /*!=================
    |    parse_F    |
    =================
    
    Parse the computed deformation gradient 
    into a form which can be read by the 
    gradient computation.
    
    */
    
    std::vector< double > reference_coords = {0,0,0,1,0,0,1,1,0,0,1,0,0.1,-0.2,1,1.1,-0.2,1.1,1.1,0.8,1.1,0.1,0.8,1}; //!The reference coordinates.
    
    micro_element::Hex8 test_element(reference_coords,U_in,U_in);  //Note: dU is a copy of U. This shouldn't matter.
    test_element.set_gpt_num(0); //Use the first gauss point
    test_element.update_shape_function_values();
    test_element.set_fundamental_measures();                       //Set the fundamental deformation measures
    tensor::Tensor23 Ftmp = test_element.get_F();
        
    //Parse the current value of the deformation gradient
    std::vector< double > out_F(9,0);
    out_F[0] = Ftmp(0,0);
    out_F[1] = Ftmp(0,1);
    out_F[2] = Ftmp(0,2);
    out_F[3] = Ftmp(1,0);
    out_F[4] = Ftmp(1,1);
    out_F[5] = Ftmp(1,2);
    out_F[6] = Ftmp(2,0);
    out_F[7] = Ftmp(2,1);
    out_F[8] = Ftmp(2,2);
        
    return out_F;
}

std::vector< std::vector< double > > compute_gradient_F(std::vector< double > U){
    /*!============================
    |    compute_gradient_F    |
    ============================
    
    Compute a numeric gradient of the 
    deformation gradient with respect to 
    the degree of freedom vector.
    
    */
    
    //Define a lambda function to parse the deformation gradient
    
    finite_difference::FiniteDifference FD(parse_F,2,U,1e-6);
    std::vector< std::vector< double > > gradient = FD.numeric_gradient();
    return gradient;
}

std::vector< double > parse_chi(std::vector< double > U_in){
    /*!===================
    |    parse_chi    |
    ===================
    
    Parse the computed micro-displacement
    tensor into a form which can be 
    read by the gradient computation.
    
    */
    
    std::vector< double > reference_coords = {0,0,0,1,0,0,1,1,0,0,1,0,0.1,-0.2,1,1.1,-0.2,1.1,1.1,0.8,1.1,0.1,0.8,1}; //!The reference coordinates.
    
    micro_element::Hex8 test_element(reference_coords,U_in,U_in);  //Note: dU is a copy of U. This shouldn't matter.
    test_element.set_gpt_num(0); //Use the first gauss point
    test_element.update_shape_function_values();
    test_element.set_fundamental_measures();                       //Set the fundamental deformation measures
    tensor::Tensor23 chitmp = test_element.get_chi();
        
    //Parse the current value of the deformation gradient
    std::vector< double > out_chi(9,0);
    out_chi[0] = chitmp(0,0);
    out_chi[1] = chitmp(0,1);
    out_chi[2] = chitmp(0,2);
    out_chi[3] = chitmp(1,0);
    out_chi[4] = chitmp(1,1);
    out_chi[5] = chitmp(1,2);
    out_chi[6] = chitmp(2,0);
    out_chi[7] = chitmp(2,1);
    out_chi[8] = chitmp(2,2);
        
    return out_chi;
}

std::vector< std::vector< double > > compute_gradient_chi(std::vector< double > U){
    /*!==============================
    |    compute_gradient_chi    |
    ==============================
    
    Compute a numeric gradient of the 
    micro-displacement tensor with 
    respect to the degree of freedom 
    vector.
    
    */
    
    //Define a lambda function to parse the deformation gradient
    
    finite_difference::FiniteDifference FD(parse_chi,2,U,1e-6);
    std::vector< std::vector< double > > gradient = FD.numeric_gradient();
    return gradient;
}

std::vector< double > parse_grad_chi(std::vector< double > U_in){
    /*!========================
    |    parse_grad_chi    |
    ========================
    
    Parse the computed gradient of 
    the micro-displacement tensor 
    into a form which can be read 
    by the gradient computation.
    
    */
    
    std::vector< double > reference_coords = {0,0,0,1,0,0,1,1,0,0,1,0,0.1,-0.2,1,1.1,-0.2,1.1,1.1,0.8,1.1,0.1,0.8,1}; //!The reference coordinates.
    
    micro_element::Hex8 test_element(reference_coords,U_in,U_in);  //Note: dU is a copy of U. This shouldn't matter.
    test_element.set_gpt_num(0); //Use the first gauss point
    test_element.update_shape_function_values();
    test_element.set_fundamental_measures();                       //Set the fundamental deformation measures
    tensor::Tensor33 grad_chitmp = test_element.get_grad_chi();
        
    //Parse the current value of the deformation gradient
    std::vector< double > out_grad_chi(27,0);
    out_grad_chi[ 0] = grad_chitmp(0,0,0);
    out_grad_chi[ 1] = grad_chitmp(0,1,0);
    out_grad_chi[ 2] = grad_chitmp(0,2,0);
    out_grad_chi[ 3] = grad_chitmp(1,0,0);
    out_grad_chi[ 4] = grad_chitmp(1,1,0);
    out_grad_chi[ 5] = grad_chitmp(1,2,0);
    out_grad_chi[ 6] = grad_chitmp(2,0,0);
    out_grad_chi[ 7] = grad_chitmp(2,1,0);
    out_grad_chi[ 8] = grad_chitmp(2,2,0);
    out_grad_chi[ 9] = grad_chitmp(0,0,1);
    out_grad_chi[10] = grad_chitmp(0,1,1);
    out_grad_chi[11] = grad_chitmp(0,2,1);
    out_grad_chi[12] = grad_chitmp(1,0,1);
    out_grad_chi[13] = grad_chitmp(1,1,1);
    out_grad_chi[14] = grad_chitmp(1,2,1);
    out_grad_chi[15] = grad_chitmp(2,0,1);
    out_grad_chi[16] = grad_chitmp(2,1,1);
    out_grad_chi[17] = grad_chitmp(2,2,1);
    out_grad_chi[18] = grad_chitmp(0,0,2);
    out_grad_chi[19] = grad_chitmp(0,1,2);
    out_grad_chi[20] = grad_chitmp(0,2,2);
    out_grad_chi[21] = grad_chitmp(1,0,2);
    out_grad_chi[22] = grad_chitmp(1,1,2);
    out_grad_chi[23] = grad_chitmp(1,2,2);
    out_grad_chi[24] = grad_chitmp(2,0,2);
    out_grad_chi[25] = grad_chitmp(2,1,2);
    out_grad_chi[26] = grad_chitmp(2,2,2);
        
    return out_grad_chi;
}

std::vector< std::vector< double > > compute_gradient_grad_chi(std::vector< double > U){
    /*!===================================
    |    compute_gradient_grad_chi    |
    ===================================
    
    Compute a numeric gradient of the 
    gradient of the micro-displacement 
    tensor with respect to the degree 
    of freedom vector.
    
    */
    
    //Define a lambda function to parse the deformation gradient
    
    finite_difference::FiniteDifference FD(parse_grad_chi,2,U,1e-6);
    std::vector< std::vector< double > > gradient = FD.numeric_gradient();
    return gradient;
}

std::vector< double > parse_C(std::vector< double > U_in){
    /*!=================
    |    parse_C    |
    =================
    
    Parse the computed right Cauchy-Green 
    deformation tensor into a form which 
    can be read by the gradient computation.
    
    */
    
    std::vector< double > reference_coords = {0,0,0,1,0,0,1,1,0,0,1,0,0.1,-0.2,1,1.1,-0.2,1.1,1.1,0.8,1.1,0.1,0.8,1}; //!The reference coordinates.
    
    micro_element::Hex8 test_element(reference_coords,U_in,U_in);  //Note: dU is a copy of U. This shouldn't matter.
    test_element.set_gpt_num(0); //Use the first gauss point
    test_element.update_shape_function_values();
    test_element.set_fundamental_measures();                       //Set the fundamental deformation measures
    test_element.set_deformation_measures();                       //Set the deformation measures
    tensor::Tensor23 Ctmp = test_element.get_C();
        
    //Parse the current value of the deformation gradient
    std::vector< double > out_C(9,0);
    out_C[0] = Ctmp(0,0);
    out_C[1] = Ctmp(0,1);
    out_C[2] = Ctmp(0,2);
    out_C[3] = Ctmp(1,0);
    out_C[4] = Ctmp(1,1);
    out_C[5] = Ctmp(1,2);
    out_C[6] = Ctmp(2,0);
    out_C[7] = Ctmp(2,1);
    out_C[8] = Ctmp(2,2);
        
    return out_C;
}

std::vector< std::vector< double > > compute_gradient_C(std::vector< double > U){
    /*!============================
    |    compute_gradient_C    |
    ============================
    
    Compute a numeric gradient of the 
    right Cauchy-Green deformation tensor  
    the degree of freedom vector.
    
    */
    
    //Define a lambda function to parse the deformation gradient
    
    finite_difference::FiniteDifference FD(parse_C,2,U,1e-6);
    std::vector< std::vector< double > > gradient = FD.numeric_gradient();
    return gradient;
}

std::vector< double > parse_Psi(std::vector< double > U_in){
    /*!===================
    |    parse_Psi    |
    ===================
    
    Parse the computed micro-deformation
    tensor into a form which can be read 
    by the gradient computation.
    
    */
    
    std::vector< double > reference_coords = {0,0,0,1,0,0,1,1,0,0,1,0,0.1,-0.2,1,1.1,-0.2,1.1,1.1,0.8,1.1,0.1,0.8,1}; //!The reference coordinates.
    
    micro_element::Hex8 test_element(reference_coords,U_in,U_in);  //Note: dU is a copy of U. This shouldn't matter.
    test_element.set_gpt_num(0); //Use the first gauss point
    test_element.update_shape_function_values();
    test_element.set_fundamental_measures();                       //Set the fundamental deformation measures
    test_element.set_deformation_measures();                       //Set the deformation measures
    tensor::Tensor23 Psitmp = test_element.get_Psi();
        
    //Parse the current value of the deformation gradient
    std::vector< double > out_Psi(9,0);
    out_Psi[0] = Psitmp(0,0);
    out_Psi[1] = Psitmp(0,1);
    out_Psi[2] = Psitmp(0,2);
    out_Psi[3] = Psitmp(1,0);
    out_Psi[4] = Psitmp(1,1);
    out_Psi[5] = Psitmp(1,2);
    out_Psi[6] = Psitmp(2,0);
    out_Psi[7] = Psitmp(2,1);
    out_Psi[8] = Psitmp(2,2);
        
    return out_Psi;
}

std::vector< std::vector< double > > compute_gradient_Psi(std::vector< double > U){
    /*!==============================
    |    compute_gradient_Psi    |
    ==============================
    
    Compute a numeric gradient of the 
    micro-deformation tensor with 
    respect to the degree of freedom 
    vector.
    
    */
    
    //Define a lambda function to parse the deformation gradient
    
    finite_difference::FiniteDifference FD(parse_Psi,2,U,1e-6);
    std::vector< std::vector< double > > gradient = FD.numeric_gradient();
    return gradient;
}

std::vector< double > parse_Gamma(std::vector< double > U_in){
    /*!=====================
    |    parse_Gamma    |
    =====================
    
    Parse the computed value of 
    gamma into a form which can 
    be read by the gradient 
    computation.
    
    */
    
    std::vector< double > reference_coords = {0,0,0,1,0,0,1,1,0,0,1,0,0.1,-0.2,1,1.1,-0.2,1.1,1.1,0.8,1.1,0.1,0.8,1}; //!The reference coordinates.
    
    micro_element::Hex8 test_element(reference_coords,U_in,U_in);  //Note: dU is a copy of U. This shouldn't matter.
    test_element.set_gpt_num(0); //Use the first gauss point
    test_element.update_shape_function_values();
    test_element.set_fundamental_measures();                       //Set the fundamental deformation measures
    test_element.set_deformation_measures();                       //Set the deformation measures.
    tensor::Tensor33 Gamma = test_element.get_Gamma();
        
    //Parse the current value of the deformation gradient
    std::vector< double > out_Gamma(27,0);
    out_Gamma[ 0] = Gamma(0,0,0);
    out_Gamma[ 1] = Gamma(0,1,0);
    out_Gamma[ 2] = Gamma(0,2,0);
    out_Gamma[ 3] = Gamma(1,0,0);
    out_Gamma[ 4] = Gamma(1,1,0);
    out_Gamma[ 5] = Gamma(1,2,0);
    out_Gamma[ 6] = Gamma(2,0,0);
    out_Gamma[ 7] = Gamma(2,1,0);
    out_Gamma[ 8] = Gamma(2,2,0);
    out_Gamma[ 9] = Gamma(0,0,1);
    out_Gamma[10] = Gamma(0,1,1);
    out_Gamma[11] = Gamma(0,2,1);
    out_Gamma[12] = Gamma(1,0,1);
    out_Gamma[13] = Gamma(1,1,1);
    out_Gamma[14] = Gamma(1,2,1);
    out_Gamma[15] = Gamma(2,0,1);
    out_Gamma[16] = Gamma(2,1,1);
    out_Gamma[17] = Gamma(2,2,1);
    out_Gamma[18] = Gamma(0,0,2);
    out_Gamma[19] = Gamma(0,1,2);
    out_Gamma[20] = Gamma(0,2,2);
    out_Gamma[21] = Gamma(1,0,2);
    out_Gamma[22] = Gamma(1,1,2);
    out_Gamma[23] = Gamma(1,2,2);
    out_Gamma[24] = Gamma(2,0,2);
    out_Gamma[25] = Gamma(2,1,2);
    out_Gamma[26] = Gamma(2,2,2);
        
    return out_Gamma;
}

std::vector< std::vector< double > > compute_gradient_Gamma(std::vector< double > U){
    /*!================================
    |    compute_gradient_Gamma    |
    ================================
    
    Compute a numeric gradient of the 
    gradient of Gamma with respect to 
    the degree of freedom vector.
    
    */
    
    //Define a lambda function to parse the deformation gradient
    
    finite_difference::FiniteDifference FD(parse_Gamma,2,U,1e-6);
    std::vector< std::vector< double > > gradient = FD.numeric_gradient();
    return gradient;
}
    
int test_fundamental_measures(std::ofstream &results){
    /*!===================================
    |    test_fundamental_measures    |
    ===================================
    
    Run tests on the fundamental deformation measures
    and related methods to ensure they are functioning 
    properly
    
    */
    
    //Seed the random number generator
    srand (1);
    
    //!Initialize test results
    int  test_num        = 6;
    bool test_results[test_num] = {false,false,false,false,false,false};
    
    //!Form the required vectors for element formation
    std::vector< double > reference_coords = {0,0,0,1,0,0,1,1,0,0,1,0,0.1,-0.2,1,1.1,-0.2,1.1,1.1,0.8,1.1,0.1,0.8,1};
    std::vector< double > Unode;
    std::vector< double > Xnode;
    std::vector< double > U;
    std::vector< double > dU;
    Xnode.resize(3);
    Unode.resize(96);
    U.resize(96);
    dU.resize(96);
    int inc = 0;
    for(int n=0; n<8; n++){
        Xnode[0] = reference_coords[0+n*3]; //Get the position of the current node
        Xnode[1] = reference_coords[1+n*3];
        Xnode[2] = reference_coords[2+n*3];
        
        Unode = test_deformation(Xnode);    //Update the deformation
        
        for(int i=0; i<12; i++){
            U[inc]  = Unode[i];   //Assign the deformation
            dU[inc] = 0.1*Unode[i];  //Assign the change in deformation (1/10 of the deformation)
            inc++;
        }
    }
    
    //!Form the hexehedral test element.
    micro_element::Hex8 element = micro_element::Hex8(reference_coords,U,dU);
    
    //!Set the gauss point
    element.set_gpt_num(0); //Use the first gauss point since the gradients should be constant
    
    //!Compute the shape function values
    element.update_shape_function_values();
    
    //!Set the fundamental deformation measures
    element.set_fundamental_measures();
    
    //!Compare the computed values to the expected result
    tensor::Tensor23 F_answer({3,3});
    tensor::Tensor23 chi_answer({3,3});
    tensor::Tensor33 grad_chi_answer({3,3,3});
    
    //!Populate the expected gradients
    std::vector< std::vector< double > > gradients = compute_gradients(element.points[0]); //Compute the numeric gradients
    
    for(int i=0; i<3; i++){
        //!Populate the deformation gradient
        F_answer(0,i)          = gradients[i][ 0];
        F_answer(1,i)          = gradients[i][ 1];
        F_answer(2,i)          = gradients[i][ 2];
        //!Populate the gradient of chi
        grad_chi_answer(0,0,i) = gradients[i][ 3];
        grad_chi_answer(1,1,i) = gradients[i][ 4];
        grad_chi_answer(2,2,i) = gradients[i][ 5];
        grad_chi_answer(1,2,i) = gradients[i][ 6];
        grad_chi_answer(0,2,i) = gradients[i][ 7];
        grad_chi_answer(0,1,i) = gradients[i][ 8];
        grad_chi_answer(2,1,i) = gradients[i][ 9];
        grad_chi_answer(2,0,i) = gradients[i][10];
        grad_chi_answer(1,0,i) = gradients[i][11];
    }
    
    //Because we are specifying the deformation, and 
    //not the actual current coordinates, we have to 
    //add 1 to the diagonal terms
    F_answer(0,0) += 1;
    F_answer(1,1) += 1;
    F_answer(2,2) += 1;
    
    //!Populate the expected chi
    chi_answer(0,0) = 1.;
    chi_answer(1,1) = 1.;
    chi_answer(2,2) = 1.;
    for(int n=0; n<8; n++){
        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                chi_answer(i,j) += element.get_N(n)*element.node_phis[n](i,j);
            }
        }
    }
    
    test_results[0] = F_answer.data.isApprox(element.get_F().data,1e-9);
    test_results[1] = chi_answer.data.isApprox(element.get_chi().data);
    test_results[2] = grad_chi_answer.data.isApprox(element.get_grad_chi().data,1e-9);
    
    //!Compare tangents
    element.set_fundamental_tangents(); //Compute the tangents for the element
    
    //!Compare the deformation gradient tangent
    std::vector< std::vector< double > > gradient = compute_gradient_F(U);
    
    //Compare the numeric and analytic tangents
    test_results[3] = true;
    
    tensor::BaseTensor<3,288> dFdU_result = element.get_dFdU();
    for(int K=0; K<96; K++){
        
        for(int I=0; I<3; I++){
            for(int J=0; J<3; J++){
                test_results[3] *= 1e-9>fabs(gradient[K][3*I+J] - dFdU_result(I,J,K));
                //std::cout << "answer: " << gradient[K][3*I+J] << "\nresult: " << dFdU_result(I,J,K) << "\n";
                if(!test_results[3]){break;}
            }
            if(!test_results[3]){break;}
        }
        if(!test_results[3]){break;}
    }
    
    //!Compare the micro-displacement tangent
    gradient = compute_gradient_chi(U);
    
    //print_vector_of_vectors("gradient",gradient);
    
    test_results[4] = true;
    
    tensor::BaseTensor<3,288> dchidU_result = element.get_dchidU();
    
    //std::cout << "dchidU:\n" << dchidU_result.data << "\n";
    
    for(int K=0; K<96; K++){
        
        for(int I=0; I<3; I++){
            for(int J=0; J<3; J++){
                test_results[4] *= 1e-9>fabs(gradient[K][3*I+J] - dchidU_result(I,J,K));
                //std::cout << "answer: " << gradient[K][3*I+J] << "\nresult: " << dchidU_result(I,J,K) << "\n";
                if(!test_results[4]){break;}
            }
            if(!test_results[4]){break;}
        }
        if(!test_results[4]){break;}
    }
    
    //!Compare the micro-displacement gradient tangent
    gradient = compute_gradient_grad_chi(U);
    
    //print_vector_of_vectors("gradient",gradient);
    
    test_results[5] = true;
    
    tensor::BaseTensor<9,288> dgrad_chidU_result = element.get_dgrad_chidU();
    
    //std::cout << "dgrad_chidU_result:\n" << dgrad_chidU_result.data << "\n";
    
    for(int I=0; I<3; I++){
        for(int J=0; J<3; J++){
            for(int K=0; K<3; K++){
                for(int L=0; L<96; L++){
                    test_results[5] *= 1e-9>fabs(gradient[L][9*K+3*I+J] - dgrad_chidU_result(I,J,K,L));
                    //std::cout << "answer: " << gradient[L][9*K+3*I+J] << "\nresult: " << dgrad_chidU_result(I,J,K,L) << "\n";
                    if(!test_results[5]){break;}
                }
                if(!test_results[5]){break;}
            }
            if(!test_results[5]){break;}
        }
        if(!test_results[5]){break;}
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
        results << "test_fundamental_measures & True\\\\\n\\hline\n";
    }
    else{
        results << "test_fundamental_measures & False\\\\\n\\hline\n";
    }
    
    return 1;
    
}

int test_deformation_measures(std::ofstream &results){
    /*!===================================
    |    test_deformation_measures    |
    ===================================
    
    Run tests on the deformation measures and 
    related methods to ensure they are functioning 
    properly
    
    */
    
    //Seed the random number generator
    srand (1);
    
    //!Initialize test results
    int  test_num        = 6;
    bool test_results[test_num] = {false,false,false,false,false,false};
    
    //!Initialize the floating point parameters
    std::vector< double > fparams(18,0.);
    
    for(int i=0; i<18; i++){
        fparams[i] = 0.1*(i+1);
    }
    
    //!Form the required vectors for element formation
    std::vector< double > reference_coords = {0,0,0,1,0,0,1,1,0,0,1,0,0.1,-0.2,1,1.1,-0.2,1.1,1.1,0.8,1.1,0.1,0.8,1};
    std::vector< double > Unode;
    std::vector< double > Xnode;
    std::vector< double > U;
    std::vector< double > dU;
    Xnode.resize(3);
    Unode.resize(96);
    U.resize(96);
    dU.resize(96);
    int inc = 0;
    for(int n=0; n<8; n++){
        Xnode[0] = reference_coords[0+n*3]; //Get the position of the current node
        Xnode[1] = reference_coords[1+n*3];
        Xnode[2] = reference_coords[2+n*3];
        
        Unode = test_deformation(Xnode);    //Update the deformation
        
        for(int i=0; i<12; i++){
            U[inc]  = Unode[i];   //Assign the deformation
            dU[inc] = 0.1*Unode[i];  //Assign the change in deformation (1/10 of the deformation)
            inc++;
        }
    }
    
    //!Form the hexehedral test element.
    micro_element::Hex8 element = micro_element::Hex8(reference_coords,U,dU,fparams);
    
    //!Set the gauss point
    element.set_gpt_num(0); //Use the first gauss point since the gradients should be constant
    
    //!Compute the shape function values
    element.update_shape_function_values();
    
    //!Set the fundamental deformation measures
    element.set_fundamental_measures();
    
    //!Set the deformation measures
    element.set_deformation_measures();
    
    //!Compare the computed values to the expected result
    tensor::Tensor23 C_answer({3,3});
    tensor::Tensor23 Psi_answer({3,3});
    tensor::Tensor33 Gamma_answer({3,3,3});
    
    for(int I=0; I<3; I++){
        for(int J=0; J<3; J++){
            
            for(int i=0; i<3; i++){
                C_answer(I,J)   += element.get_F()(i,I)*element.get_F()(i,J);
                Psi_answer(I,J) += element.get_F()(i,I)*element.get_chi()(i,J);
            }
            
        }
    }
    
    for(int I=0; I<3; I++){
        for(int J=0; J<3; J++){
            for(int K=0; K<3; K++){
                for(int i=0; i<3; i++){
                    Gamma_answer(I,J,K)   += element.get_F()(i,I)*element.get_grad_chi()(i,J,K);
                }
            }            
        }
    }
    
    test_results[0] = C_answer.data.isApprox(element.get_C().data);
    test_results[1] = Psi_answer.data.isApprox(element.get_Psi().data);
    test_results[2] = Gamma_answer.data.isApprox(element.get_Gamma().data);
    
    //!Compare tangents
    
    element.set_fundamental_tangents();
    element.set_deformation_tangents();
    
    //!Compare the deformation gradient tangent
    std::vector< std::vector< double > > gradient = compute_gradient_C(U);
    
    //print_vector_of_vectors("gradient",gradient);
    
    //Compare the numeric and analytic tangents
    test_results[3] = true;
    
    tensor::BaseTensor<3,288> dCdU_result = element.get_dCdU();
    
    //std::cout << "dCdU:\n" << dCdU_result.data << "\n";
    
    for(int K=0; K<96; K++){
        
        for(int I=0; I<3; I++){
            for(int J=0; J<3; J++){
                test_results[3] *= 1e-9>fabs(gradient[K][3*I+J] - dCdU_result(I,J,K));
                //std::cout << "answer: " << gradient[K][3*I+J] << "\nresult: " << dFdU_result(I,J,K) << "\n";
                if(!test_results[3]){break;}
            }
            if(!test_results[3]){break;}
        }
        if(!test_results[3]){break;}
    }
    
    //!Compare the micro-deformation tangent
    gradient = compute_gradient_Psi(U);
    
    //print_vector_of_vectors("gradient",gradient);
    
    test_results[4] = true;
    
    tensor::BaseTensor<3,288> dPsidU_result = element.get_dPsidU();
    
    //std::cout << "dPsidU:\n" << dPsidU_result.data << "\n";
    
    for(int K=0; K<96; K++){
        
        for(int I=0; I<3; I++){
            for(int J=0; J<3; J++){
                test_results[4] *= 1e-9>fabs(gradient[K][3*I+J] - dPsidU_result(I,J,K));
                //std::cout << "answer: " << gradient[K][3*I+J] << "\nresult: " << dchidU_result(I,J,K) << "\n";
                if(!test_results[4]){break;}
            }
            if(!test_results[4]){break;}
        }
        if(!test_results[4]){break;}
    }
    
    //!Compare the gradient of Gamma tangent
    gradient = compute_gradient_Gamma(U);
    
    //print_vector_of_vectors("gradient",gradient);
    
    test_results[5] = true;
    
    tensor::BaseTensor<9,288> dGammadU_result = element.get_dGammadU();
    
    //std::cout << "dGammadU_result:\n" << dGammadU_result.data << "\n";
    
    for(int I=0; I<3; I++){
        for(int J=0; J<3; J++){
            for(int K=0; K<3; K++){
                for(int L=0; L<96; L++){
                    test_results[5] *= 1e-8>fabs(gradient[L][9*K+3*I+J] - dGammadU_result(I,J,K,L));
                    //std::cout << "answer: " << gradient[L][9*K+3*I+J] << "\nresult: " << dGammadU_result(I,J,K,L) << "\n";
                    //std::cout << I << J << K << L << "\n";
                    if(!test_results[5]){break;}
                }
                if(!test_results[5]){break;}
            }
            if(!test_results[5]){break;}
        }
        if(!test_results[5]){break;}
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
        results << "test_deformation_measures & True\\\\\n\\hline\n";
    }
    else{
        results << "test_deformation_measures & False\\\\\n\\hline\n";
    }
    
    return 1;
}

int test_balance_of_linear_momentum(std::ofstream &results){
    /*!=========================================
    |    test_balance_of_linear_momentum    |
    =========================================
    
    Run tests on the formation of the residual of the 
    balance of linear momentum.
    
    */
    
    //Seed the random number generator
    srand (1);
    
    //!Initialize test results
    int  test_num        = 1;
    bool test_results[test_num] = {false};
    
    //!Initialize the floating point parameters
    std::vector< double > fparams(18,0.);
    
    for(int i=0; i<18; i++){
        fparams[i] = 0.1*(i+1);
    }
    
    //!Form the required vectors for element formation
    std::vector< double > reference_coords = {0,0,0,1,0,0,1,1,0,0,1,0,0.1,-0.2,1,1.1,-0.2,1.1,1.1,0.8,1.1,0.1,0.8,1};
    std::vector< double > Unode;
    std::vector< double > Xnode;
    std::vector< double > U;
    std::vector< double > dU;
    Xnode.resize(3);
    Unode.resize(96);
    U.resize(96);
    dU.resize(96);
    int inc = 0;
    for(int n=0; n<8; n++){
        Xnode[0] = reference_coords[0+n*3]; //Get the position of the current node
        Xnode[1] = reference_coords[1+n*3];
        Xnode[2] = reference_coords[2+n*3];
        
        Unode = test_deformation(Xnode);    //Update the deformation
        
        for(int i=0; i<12; i++){
            U[inc]  = Unode[i];   //Assign the deformation
            dU[inc] = 0.1*Unode[i];  //Assign the change in deformation (1/10 of the deformation)
            inc++;
        }
    }
    
    //!Form the hexehedral test element.
    micro_element::Hex8 element = micro_element::Hex8(reference_coords,U,dU,fparams);
    
    //!Set the gauss point
    element.set_gpt_num(0); //Use the first gauss point since the gradients should be constant
    
    //!Compute the shape function values
    element.update_shape_function_values();
    
    //!Set the fundamental deformation measures
    element.set_fundamental_measures();
    
    //!Set the deformation measures
    element.set_deformation_measures();
    
    //!Set the stresses at the gauss point
    element.set_stresses();
    
    //!Set the residual vector due to the given gauss point
    element.add_internal_nodal_forces();
    
    //!Compute the expected residual vector
    std::vector< double > RHS(96,0.);
    
    //!Compute the internal forces
    for(int n=0; n<8; n++){
        for(int j=0; j<3; j++){
            for(int I=0; I<3; I++){
                for(int J=0; J<3; J++){
                    RHS[j+12*n] += -element.get_dNdx(0,n)[I]*element.PK2[0](I,J)*element.get_F()(j,J)*element.get_Jhatdet(0)*element.weights[0];
                }
            }
        }
    }
    
    //!Compare the expected results to the element results
    test_results[0] = true;
    for(int i=0; i<96; i++){
        test_results[0] *= 1e-9>fabs(element.RHS[i]-RHS[i]);
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
        results << "test_balance_of_linear_momentum & True\\\\\n\\hline\n";
    }
    else{
        results << "test_balance_of_linear_momentum & False\\\\\n\\hline\n";
    }
    
    return 1;
}

int test_balance_of_first_moment_of_momentum(std::ofstream &results){
    /*!==================================================
    |    test_balance_of_first_moment_of_momentum    |
    ==================================================
    
    Run tests on the formation of the residual of the 
    first moment of momentum.
    
    */
    
    //Seed the random number generator
    srand (1);
    
    //!Initialize test results
    int  test_num        = 1;
    bool test_results[test_num] = {false};
    
    //!Initialize the floating point parameters
    std::vector< double > fparams(18,0.);
    
    for(int i=0; i<18; i++){
        fparams[i] = 0.1*(i+1);
    }
    
    //!Form the required vectors for element formation
    std::vector< double > reference_coords = {0,0,0,1,0,0,1,1,0,0,1,0,0.1,-0.2,1,1.1,-0.2,1.1,1.1,0.8,1.1,0.1,0.8,1};
    std::vector< double > Unode;
    std::vector< double > Xnode;
    std::vector< double > U;
    std::vector< double > dU;
    Xnode.resize(3);
    Unode.resize(96);
    U.resize(96);
    dU.resize(96);
    int inc = 0;
    for(int n=0; n<8; n++){
        Xnode[0] = reference_coords[0+n*3]; //Get the position of the current node
        Xnode[1] = reference_coords[1+n*3];
        Xnode[2] = reference_coords[2+n*3];
        
        Unode = test_deformation(Xnode);    //Update the deformation
        
        for(int i=0; i<12; i++){
            U[inc]  = Unode[i];   //Assign the deformation
            dU[inc] = 0.1*Unode[i];  //Assign the change in deformation (1/10 of the deformation)
            inc++;
        }
    }
    
    //!Form the hexehedral test element.
    micro_element::Hex8 element = micro_element::Hex8(reference_coords,U,dU,fparams);
    
    //!Set the gauss point
    element.set_gpt_num(0); //Use the first gauss point since the gradients should be constant
    
    //!Compute the shape function values
    element.update_shape_function_values();
    
    //!Set the fundamental deformation measures
    element.set_fundamental_measures();
    
    //!Set the deformation measures
    element.set_deformation_measures();
    
    //!Set the stresses at the gauss point
    element.set_stresses();
    
    //!Set the residual vector due to the given gauss point
    element.add_internal_moments();
    
    //!Compute the expected residual vector
    std::vector< double > RHS(96,0.);
    
    //!Define the internal stress balance
    tensor::Tensor23 mu_int({3,3});
    
    //!Compute the internal stresses
    for(int n=0; n<8; n++){
        mu_int.data.setZero();
        
        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                for(int I=0; I<3; I++){
                    for(int J=0; J<3; J++){
                        mu_int(i,j) += -element.get_N(n)*element.get_F()(i,I)*(element.SIGMA[0](I,J) - element.PK2[0](I,J))*element.get_F()(j,J)*element.get_Jhatdet(0)*element.weights[0];
                    }
                }
                for(int I=0; I<3; I++){
                    for(int J=0; J<3; J++){
                        for(int K=0; K<3; K++){
                            mu_int(i,j) += -element.get_dNdx(0,n)[K]*element.get_F()(j,J)*element.get_chi()(i,I)*element.M[0](K,J,I)*element.get_Jhatdet(0)*element.weights[0];
                        }
                    }
                }
            }
        }
        
        RHS[0+n*12+3] = mu_int(0,0);
        RHS[1+n*12+3] = mu_int(1,1);
        RHS[2+n*12+3] = mu_int(2,2);
        RHS[3+n*12+3] = mu_int(1,2);
        RHS[4+n*12+3] = mu_int(0,2);
        RHS[5+n*12+3] = mu_int(0,1);
        RHS[6+n*12+3] = mu_int(2,1);
        RHS[7+n*12+3] = mu_int(2,0);
        RHS[8+n*12+3] = mu_int(1,0);
        
    }
    
    //!Compare the expected results to the element results
    test_results[0] = true;
    for(int i=0; i<96; i++){
        test_results[0] *= 1e-9>fabs(element.RHS[i]-RHS[i]);
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
        results << "test_balance_of_linear_momentum & True\\\\\n\\hline\n";
    }
    else{
        results << "test_balance_of_linear_momentum & False\\\\\n\\hline\n";
    }
    
    return 1;
}

int test_integrate_element(std::ofstream &results){
    /*!================================
    |    test_integrate_element    |
    ================================
    
    Run tests on the integration of the 
    finite element.
    
    */
    
    //Seed the random number generator
    srand (1);
    
    //!Initialize test results
    int  test_num        = 1;
    bool test_results[test_num] = {false};
    
    //!Initialize the floating point parameters
    std::vector< double > fparams(18,0.);
    
    for(int i=0; i<18; i++){
        fparams[i] = 0.1*(i+1);
    }
    
    //!Form the required vectors for element formation
    std::vector< double > reference_coords = {0,0,0,1,0,0,1,1,0,0,1,0,0.1,-0.2,1,1.1,-0.2,1.1,1.1,0.8,1.1,0.1,0.8,1};
    std::vector< double > Unode;
    std::vector< double > Xnode;
    std::vector< double > U;
    std::vector< double > dU;
    Xnode.resize(3);
    Unode.resize(96);
    U.resize(96);
    dU.resize(96);
    int inc = 0;
    for(int n=0; n<8; n++){
        Xnode[0] = reference_coords[0+n*3]; //Get the position of the current node
        Xnode[1] = reference_coords[1+n*3];
        Xnode[2] = reference_coords[2+n*3];
        
        Unode = test_deformation(Xnode);    //Update the deformation
        
        for(int i=0; i<12; i++){
            U[inc]  = Unode[i];   //Assign the deformation
            dU[inc] = 0.1*Unode[i];  //Assign the change in deformation (1/10 of the deformation)
            inc++;
        }
    }
    
    //std::clock_t start;
    //double duration;

    
    
    //!Form the hexehedral test element.
    micro_element::Hex8 element = micro_element::Hex8(reference_coords,U,dU,fparams);
    
    //start = std::clock();
    
    //!Integrate the element
    element.integrate_element();
    
    //duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    //std::cout<<"printf: "<< duration <<'\n';
    
    //!Compute the expected residual vector
    std::vector< double > RHS(96,0.);
    
    //!Define the internal stress balance
    tensor::Tensor23 mu_int({3,3});
    
    //!Integrate the element
    for(int gpt_num=0; gpt_num<8; gpt_num++){
        element.set_gpt_num(gpt_num);
        element.update_gauss_point();
        
        //!Compute the internal stresses
        for(int n=0; n<8; n++){
            
            //!Add the internal force residual
            for(int j=0; j<3; j++){
                for(int I=0; I<3; I++){
                    for(int J=0; J<3; J++){
                        RHS[j+12*n] += -element.get_dNdx(0,n)[I]*element.PK2[gpt_num](I,J)*element.get_F()(j,J)*element.get_Jhatdet(0)*element.weights[gpt_num];
                    }
                }
            }
            
            //!Add the internal stress residual
            mu_int.data.setZero();
        
            for(int i=0; i<3; i++){
                for(int j=0; j<3; j++){
                    for(int I=0; I<3; I++){
                        for(int J=0; J<3; J++){
                            mu_int(i,j) += -element.get_N(n)*element.get_F()(i,I)*(element.SIGMA[gpt_num](I,J) - element.PK2[gpt_num](I,J))*element.get_F()(j,J)*element.get_Jhatdet(0)*element.weights[gpt_num];
                        }
                    }
                    for(int I=0; I<3; I++){
                        for(int J=0; J<3; J++){
                            for(int K=0; K<3; K++){
                                mu_int(i,j) += -element.get_dNdx(0,n)[K]*element.get_F()(j,J)*element.get_chi()(i,I)*element.M[gpt_num](K,J,I)*element.get_Jhatdet(0)*element.weights[gpt_num];
                            }
                        }
                    }
                }
            }
            
            RHS[0+n*12+3] += mu_int(0,0);
            RHS[1+n*12+3] += mu_int(1,1);
            RHS[2+n*12+3] += mu_int(2,2);
            RHS[3+n*12+3] += mu_int(1,2);
            RHS[4+n*12+3] += mu_int(0,2);
            RHS[5+n*12+3] += mu_int(0,1);
            RHS[6+n*12+3] += mu_int(2,1);
            RHS[7+n*12+3] += mu_int(2,0);
            RHS[8+n*12+3] += mu_int(1,0);
        
        }
    }
    
    //!Compare the expected results to the element results
    test_results[0] = true;
    for(int i=0; i<96; i++){
        test_results[0] *= 1e-9>fabs(element.RHS[i]-RHS[i]);
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
        results << "test_balance_of_linear_momentum & True\\\\\n\\hline\n";
    }
    else{
        results << "test_balance_of_linear_momentum & False\\\\\n\\hline\n";
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
    test_fundamental_measures(results);
    test_deformation_measures(results);
    test_balance_of_linear_momentum(results);
    test_balance_of_first_moment_of_momentum(results);
    test_integrate_element(results);
    
    //Close the results file
    results.close();
}

