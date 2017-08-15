/*!=======================================================
  |                                                     |
  |                 test_tensor.cpp                     |
  |                                                     |
  -------------------------------------------------------
  | The unit test file for tensor.h/cpp.                |
  | This file tests the classes and functions defined   |
  | in tensor.h/cpp.                                    |
  |                                                     |
  | Generated files:                                    |
  |    Results.tex:  A LaTeX file which contains the    |
  |                  results as they will be included   |
  |                  in the generated report.           |
  =======================================================
  | Dependencies:                                       |
  | Eigen: An implementation of various matrix          |
  |        commands. The implementation of the data     |
  |        matrix uses such a matrix.                   |
  =======================================================*/
  
#include <functional>
#include <iostream>
#include <fstream>
#include <numeric>
#include <vector>
#include <Eigen/Dense>
#include <tensor.h>
#include <ctime>

int test_tensor_functionality(std::ofstream &results){
    /*!====================================
    |       test_tensor_functionality   |
    =====================================
    
    A test of some of the basic tensor functionality. 
    This should show that tensors of different orders 
    can be generated and manipulated.*/

    int  test_num        = 6;
    bool test_results[test_num] = {false,false,false,false,false,false};
    
    //Define a test vector
    std::vector< int > v_shape;                 //Vector shape vector
    v_shape.resize(1);                          //Resizing the shape vector to have only one index
    v_shape[0] = 6;                             //Setting the vector size
    tensor::Tensor V = tensor::Tensor(v_shape); //Form the vector
    
    Eigen::MatrixXd V_compare(6,1); //Initialize the comparison matrix
    
    V_compare << 1,2,3,4,5,6;  //Set the initial values
    
    double inc = 1; //Set an initial increment vector
    
    //!Compare expected storage pattern vs. actual for a vector
    for(int i=0; i<v_shape[0]; i++){//Iterate through the vector setting the required values
        V(i) = inc;                 //Set the i'th value of V equal to inc
        inc++;                      //Increment inc
    }
    
    test_results[0] = V_compare.isApprox(V.data); //Check the results of the test
    
    //!Setting index test for a vector
    V_compare(3) = -1; //Set the 4th index value to -1
    V(3)         = -1; //Do the same for the tensor
    
    test_results[1] = V_compare.isApprox(V.data); //Check the results of the test
    
    //Define a test matrix
    std::vector < int > m_shape; //Matrix shape vector
    m_shape.resize(2);           //Resizing the shape vector to having two indices
    
    m_shape[0] = 3;              //Create a 3x3 matrix
    m_shape[1] = 3;
    
    tensor::Tensor M = tensor::Tensor(m_shape); //Initialize the matrix
    
    Eigen::MatrixXd M_compare(3,3); //Initialize the comparison matrix
    M_compare << 1,2,3,4,5,6,7,8,9;  //Set the initial values
    
    
    //!Compare expected storage pattern to actual for a 2nd order tensor
    inc = 1; //Set an initial increment vector
    
    for(int i=0; i<m_shape[0]; i++){
        for(int j=0; j<m_shape[1]; j++){
            M(i,j) = inc;
            inc++;
        }
    }
    
    test_results[2] =  M_compare.isApprox(M.data); //Check the results of the test
    
    //!Setting index test for a 2nd order tensor
    M_compare(2,1) = -8;
    M(2,1)         = -8;
    
    test_results[3] =  M_compare.isApprox(M.data); //Check the results of the test
    
    //Define a test 3rd order tensor
    std::vector < int > t_shape;                //Initialize the the tensor shape vector
    t_shape.resize(3);                          //Resize the shape vector of the tensor
    
    t_shape[0] = 4;                             //Initialize the tensor size
    t_shape[1] = 3;
    t_shape[2] = 7;
    
    tensor::Tensor T = tensor::Tensor(t_shape); //Initialize the tensor
    
    Eigen::MatrixXd T_compare(4,21); //Initialize the comparison tensor
    T_compare <<  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,
                 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42,
                 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,
                 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84;
    
    //!Compare expected storage pattern to actual for a 3rd order tensor
    inc = 1; //Set an initial increment vector
    
    for(int i=0; i<t_shape[0]; i++){
        for(int j=0; j<t_shape[1]; j++){
            for(int k=0; k<t_shape[2]; k++){
                T(i,j,k) = inc;
                inc++;
            }
        }
    }
    
    test_results[4] =  T_compare.isApprox(T.data); //Check the results of the test
    
    //!Setting index test for a 3nd order tensor
    T_compare(2,13) = -8;
    T(2,1,6)        = -8;
    
    test_results[5] =  T_compare.isApprox(T.data); //Check the results of the test
    
    //std::cout << "\nT_compare:\n" << T_compare << "\n" << "T:\n" << T.data << "\n";
    
    //Compare all test results
    bool tot_result = true;
    for(int i = 0; i<test_num; i++){
        //std::cout << "\nSub-test " << i+1 << " result: " << test_results[i] << "\n";
        if(!test_results[i]){
            tot_result = false;
        }
    }
    
    if(tot_result){
        results << "test_tensor_functionality & True\\\\\n\\hline\n";
    }
    else{
        results << "test_tensor_functionality & False\\\\\n\\hline\n";
    }
    
    return 1;
    
}

int test_eye(std::ofstream &results){
    /*!========================
    |       test_eye       |
    ========================
    
    A test of the generation of the second
    order identity tensor. Note that with 
    different storage schemes this tensor 
    may not always be the same in terms of 
    the way it is stored.*/
    
    //Initialize the results
    int  test_num        = 2;
    bool test_results[test_num] = {false,false};
    
    //Compute the identity tensor
    tensor::Tensor I = tensor::eye();
    
    //Define a test matrix
    std::vector < int > m_shape; //Matrix shape vector
    m_shape.resize(2);           //Resizing the shape vector to having two indices
    
    m_shape[0] = 3;              //Create a 3x3 tensor
    m_shape[1] = 3;
    
    tensor::Tensor M = tensor::Tensor(m_shape); //Initialize the matrix
    M.data << 1,2,3,4,5,6,7,8,9;  //Set the initial values
    
    //!Run multiplication tests
    tensor::Tensor R1 = tensor::Tensor(m_shape); //Initialize the first results matrix
    tensor::Tensor R2 = tensor::Tensor(m_shape); //Initialize the second results matrix
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            for(int k=0; k<3; k++){
                R1(i,j) += I(i,k)*M(k,j);
                R2(i,j) += M(i,k)*I(k,j);
            }
        }
    }
    
    //Check the results of the tests
    test_results[0] = R1.data.isApprox(M.data);
    test_results[1] = R2.data.isApprox(M.data);
    
    //Compare all test results
    bool tot_result = true;
    for(int i = 0; i<test_num; i++){
        //std::cout << "\nSub-test " << i+1 << " result: " << test_results[i] << "\n";
        if(!test_results[i]){
            tot_result = false;
        }
    }
    
    if(tot_result){
        results << "test_eye & True\\\\\n\\hline\n";
    }
    else{
        results << "test_eye & False\\\\\n\\hline\n";
    }
    
}

int test_FOT_eye(std::ofstream &results){
    /*!========================
    |     test_FOT_eye     |
    ========================
    
    A test of the generation of the fourth
    order identity tensor. Note that with 
    different storage schemes this tensor 
    may not always be the same in terms of 
    the way it is stored.*/
    
    //Initialize the results
    int  test_num        = 2;
    bool test_results[test_num] = {false,false};
    
    //Compute the identity tensor
    tensor::Tensor FOTI = tensor::FOT_eye();
    
    //Define a test matrix
    std::vector < int > m_shape; //Matrix shape vector
    m_shape.resize(4);           //Resizing the shape vector to having two indices
    
    m_shape[0] = 3;              //Create a 3x3x3x3 tensor
    m_shape[1] = 3;
    m_shape[2] = 3;
    m_shape[3] = 3;
    
    tensor::Tensor FOT = tensor::Tensor(m_shape); //Initialize the matrix
    FOT.data << 2,4,3,5,1,3,4,6,1,
                4,5,2,7,2,5,3,5,7,
                6,2,8,3,6,2,6,4,5,
                7,2,9,1,5,3,7,3,5,
                9,8,9,4,1,2,4,3,5,
                4,2,5,6,2,7,4,5,3,
                6,6,5,2,4,1,5,4,8,
                9,1,3,2,2,4,5,6,7,
                5,2,1,2,1,1,4,5,6;
    
    //!Run multiplication tests
    tensor::Tensor R1 = tensor::Tensor(m_shape); //Initialize the first results matrix
    tensor::Tensor R2 = tensor::Tensor(m_shape); //Initialize the second results matrix
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            for(int k=0; k<3; k++){
                for(int l=0; l<3; l++){
                    for(int m=0; m<3; m++){
                        for(int n=0; n<3; n++){
                            R1(i,j,k,l) += FOTI(i,j,m,n)*FOT(m,n,k,l);
                            R2(i,j,k,l) += FOT(i,j,m,n)*FOTI(m,n,k,l);
                        }
                    }
                }
            }
        }
    }
    
    //Check the results of the tests
    test_results[0] = R1.data.isApprox(FOT.data);
    test_results[1] = R2.data.isApprox(FOT.data);
    
    //Compare all test results
    bool tot_result = true;
    for(int i = 0; i<test_num; i++){
        //std::cout << "\nSub-test " << i+1 << " result: " << test_results[i] << "\n";
        if(!test_results[i]){
            tot_result = false;
        }
    }
    
    if(tot_result){
        results << "test_FOT_eye & True\\\\\n\\hline\n";
    }
    else{
        results << "test_FOT_eye & False\\\\\n\\hline\n";
    }
    
}

int test_inverse(std::ofstream &results){
    /*!========================
    |     test_inverse     |
    ========================
    
    A test of the inverse of the tensor computation. 
    Several different tensor formulations are compared 
    and examined to make sure that the product of the 
    inverse and the original tensor is the identity 
    tensor.*/
    
    int  test_num               = 2;
    bool test_results[test_num] = {false,false};
    
    /*!Test the inverse of a second order tensor*/
    std::vector< int > m_shape; //Initialize the shape vector
    m_shape.resize(2);          //Resize the shape vector to a second order tensor
    
    m_shape[0] = 3;             //Set the tensor dimensions
    m_shape[1] = 3;
    
    tensor::Tensor T = tensor::Tensor(m_shape); //Initialize and populate the tensor
    T.data << 2,4,3,5,1,3,4,6,1;
    
    tensor::Tensor Tinv = T.inverse(); //Invert the tensor
    
    tensor::Tensor product = tensor::Tensor(m_shape); //Compute the product of the tensor and its inverse
    
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            for(int k=0; k<3; k++){
                product(i,j) += T(i,k)*Tinv(k,j);
            }
        }
    }
    
    tensor::Tensor I = tensor::eye(); //Get the second order identity tensor
    
    test_results[0] = I.data.isApprox(product.data);
    
    /*!Test the inverse of a fourth order tensor*/
    std::vector< int > fot_shape; //Initialize the shape vector
    fot_shape.resize(4);          //Resize the shape vector to a fourth order tensor
    
    fot_shape[0] = 3;             //Set the tensor dimensions
    fot_shape[1] = 3;
    fot_shape[2] = 3;
    fot_shape[3] = 3;
    
    tensor::Tensor FOT = tensor::Tensor(fot_shape); //Initialize and populate the tensor
    
    FOT.data << 2,4,3,5,1,3,4,6,1,
                4,5,2,7,2,5,3,5,7,
                6,2,8,3,6,2,6,4,5,
                7,2,9,1,5,3,7,3,5,
                9,8,9,4,1,2,4,3,5,
                4,2,5,6,2,7,4,5,3,
                6,6,5,2,4,1,5,4,8,
                9,1,3,2,2,4,5,6,7,
                5,2,1,2,1,1,4,5,6;
    
    tensor::Tensor FOTinv = FOT.inverse(); //Invert the tensor
    
    tensor::Tensor productFOT = tensor::Tensor(fot_shape); //Compute the product of the tensor and its inverse
    
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            for(int k=0; k<3; k++){
                for(int l=0; l<3; l++){
                    for(int m=0; m<3; m++){
                        for(int n=0; n<3; n++){
                            productFOT(i,j,k,l) += FOT(i,j,m,n)*FOTinv(m,n,k,l);
                        }
                    }
                }
            }
        }
    }
    tensor::Tensor FOTI = tensor::FOT_eye();
    test_results[1] = FOTI.data.isApprox(productFOT.data);
    
    //Compare all test results
    bool tot_result = true;
    for(int i = 0; i<test_num; i++){
        //std::cout << "\nSub-test " << i+1 << " result: " << test_results[i] << "\n";
        if(!test_results[i]){
            tot_result = false;
        }
    }
    
    if(tot_result){
        results << "test_inverse & True\\\\\n\\hline\n";
    }
    else{
        results << "test_inverse & False\\\\\n\\hline\n";
    }
    return 1;
}

int test_det(std::ofstream &results){
    /*!========================
    |       test_det       |
    ========================
    
    A test of the determinant of the tensor computation. 
    Tensors with known determinants are compared to 
    those computed.*/
    
    int  test_num               = 2;
    bool test_results[test_num] = {false,false};
    
    /*!Test the determinant of a second order tensor*/
    std::vector< int > m_shape; //Initialize the shape vector
    m_shape.resize(2);          //Resize the shape vector to a second order tensor
    
    m_shape[0] = 3;             //Set the tensor dimensions
    m_shape[1] = 3;
    
    tensor::Tensor T = tensor::Tensor(m_shape); //Initialize and populate the tensor
    T.data << 2,4,3,5,1,3,4,6,1;
    
    double Tdet = T.det(); //Compute the determinant
    
    if(fabs(Tdet-72.)<1e-6){test_results[0] = true;}
    
    /*!Test the determinant of a fourth order tensor*/
    std::vector< int > fot_shape; //Initialize the shape vector
    fot_shape.resize(4);          //Resize the shape vector to a fourth order tensor
    
    fot_shape[0] = 3;             //Set the tensor dimensions
    fot_shape[1] = 3;
    fot_shape[2] = 3;
    fot_shape[3] = 3;
    
    tensor::Tensor FOT = tensor::Tensor(fot_shape); //Initialize and populate the tensor
    
    FOT.data << 2,4,3,5,1,3,4,6,1,
                4,5,2,7,2,5,3,5,7,
                6,2,8,3,6,2,6,4,5,
                7,2,9,1,5,3,7,3,5,
                9,8,9,4,1,2,4,3,5,
                4,2,5,6,2,7,4,5,3,
                6,6,5,2,4,1,5,4,8,
                9,1,3,2,2,4,5,6,7,
                5,2,1,2,1,1,4,5,6;
    
    double FOTdet = FOT.det(); //Compute the determinant
    
    tensor::Tensor productFOT = tensor::Tensor(fot_shape); //Compute the product of the tensor and its inverse
    
    if(fabs(FOTdet-226209.)<1e-6){test_results[1] = true;}
    
    //Compare all test results
    bool tot_result = true;
    for(int i = 0; i<test_num; i++){
        //std::cout << "\nSub-test " << i+1 << " result: " << test_results[i] << "\n";
        if(!test_results[i]){
            tot_result = false;
        }
    }
    
    if(tot_result){
        results << "test_inverse & True\\\\\n\\hline\n";
    }
    else{
        results << "test_inverse & False\\\\\n\\hline\n";
    }
    return 1;
}

int test_operators(std::ofstream &results){
    /*!==============================
    |       test_operators       |
    ==============================
    
    A test of the operators on the tensor 
    object. These include, addition, subtraction, 
    +=, and -=*/
    
    //Define a test matrix
    std::vector < int > m_shape1; //Matrix shape vector
    m_shape1.resize(2);           //Resizing the shape vector to having two indices
    
    m_shape1[0] = 3;              //Create a 3x3 tensor
    m_shape1[1] = 3;
    
    tensor::Tensor M1 = tensor::Tensor(m_shape1); //Initialize the matrix
    M1.data << 1,2,3,4,5,6,7,8,9;  //Set the initial values
    
    //Define a second test matrix
    std::vector < int > m2_shape; //Matrix shape vector
    m2_shape.resize(2);           //Resizing the shape vector to having two indices
    
    m2_shape[0] = 3;              //Create a 3x3 tensor
    m2_shape[1] = 3;
    
    tensor::Tensor M2 = tensor::Tensor(m2_shape); //Initialize the matrix
    M2.data << 6,2,2,7,2,1,3,9,0;  //Set the initial values
    
    //Compute the values
    tensor::Tensor Tsum =  M1+M2;
    tensor::Tensor Tsub =  M1-M2;
    tensor::Tensor Tneg = -M1;
    M1 += M2;
    M2 -= M2;
    
    //Set the answers
    Eigen::MatrixXd MA1 = Eigen::MatrixXd(3,3) << 1+6,2+2,3+2,4+7,5+2,6+1,7+3,8+9,9+0;
    Eigen::MatrixXd MA2 = Eigen::MatrixXd(3,3) << 1-6,2-2,3-2,4-7,5-2,6-1,7-3,8-9,9-0;
    Eigen::MatrixXd MA3 = Eigen::MatrixXd(3,3) << -1,-2,-3,-4,-5,-6,-7,-8,-9;
    Eigen::MatrixXd MA4 = Eigen::MatrixXd(3,3) << 1+6,2+2,3+2,4+7,5+2,6+1,7+3,8+9,9+0;
    Eigen::MatrixXd MA5 = Eigen::Zeros(3,3);
    
    //Compare the results to the solutions
    int test_num = 5;
    std::array<bool, test_num> test_results = {false,false,false,false,false};
    
    test_results[0] = MA1.isApprox(Tsum.data);
    test_results[1] = MA2.isApprox(Tsub.data);
    test_results[2] = MA3.isApprox(Tneg.data);
    test_results[3] = MA4.isApprox(M1.data);
    test_results[4] = MA5.isApprox(M2.data);
    
    //Compare all test results
    bool tot_result = true;
    for(int i = 0; i<test_num; i++){
        std::cout << "\nSub-test " << i+1 << " result: " << test_results[i] << "\n";
        if(!test_results[i]){
            tot_result = false;
        }
    }
    
    if(tot_result){
        results << "test_inverse & True\\\\\n\\hline\n";
    }
    else{
        results << "test_inverse & False\\\\\n\\hline\n";
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
    test_tensor_functionality(results);
    test_inverse(results);
    test_eye(results);
    test_FOT_eye(results);
    test_det(results);
    
    //Close the results file
    results.close();
}

