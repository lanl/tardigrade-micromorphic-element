/*!=======================================================
  |                                                     |
  |            test_finite_difference.cpp               |
  |                                                     |
  -------------------------------------------------------
  | The unit test file for finite_difference.h          |
  =======================================================*/
  
#include <functional>
#include <iostream>
#include <fstream>
#include <numeric>
#include <vector>
#include <Eigen/Dense>
#include <finite_difference.h>
#include <ctime>
#include <stdlib.h>

void print_matrix(std::vector< std::vector< double > > M){
    /*!======================
    |    print_matrix    |
    ======================
    
    Print a matrix out
    
    */
    for(int i=0; i<M.size(); i++){
        for(int j=0; j<M[i].size(); j++){
            std::cout << M[i][j] << " ";
        }
        std::cout << "\n";
    }
    return;
}

bool compare_vectors(const std::vector<double> &V1, const std::vector<double> &V2){
    /*!=========================
    |    compare_vectors    |
    =========================
    
    Compare two vectors and determine if they are 
    equal.
    
    */
    
    //!Check vector's size
    if(V1.size()!=V2.size()){return false;}
    
    //!Compare the values
    else{
        for(int i=0; i<V1.size(); i++){
            if(!(1e-9>fabs(V1[i]-V2[i]))){return false;}
        }
        return true;
    }
}

bool compare_vectors_of_vectors(const std::vector< std::vector< double > > &V1, const std::vector< std::vector< double > > &V2){
    /*!====================================
    |    compare_vectors_of_vectors    |
    ====================================
    
    Compare two vectors of vectors and determine if they are equal.
    
    */
    
    bool result = true;     //!The overall result
    
    //!Check vector's size
    if(V1.size()!=V2.size()){std::cout << "\nHI!\n"; return false;}
    
    //Compare each column against the value in the other column
    for(int i=0; i<V1.size(); i++){result *= compare_vectors(V1[i],V2[i]); if(result==false){return false;}}
    
    return result;
}

std::vector<double> test_fxn_1(std::vector<double> x){
    /*! test_fxn_1
    
    A test function used in test_finite_difference
    
    */
    
    return {x[0]*x[0]};
}

std::vector<double> answer_fxn_1(std::vector<double> x){
    /*! answer_fxn_1
    
    A solution function used in test_finite_difference 
        
    */
    
    return {2.*x[0]};
}

std::vector<double> test_fxn_2(std::vector<double> x){
    /*! test_fxn_2
    
    A test function used in test_finite_difference 
    and test_numeric_gradient
    
    */
    
    return {x[0]*x[0],x[1]-1};
}

std::vector< std::vector< double > > answer_fxn_2(std::vector<double> x){
    /*! answer_fxn_2
    
    A solution function used in test_finite_difference
    and test_numeric_gradient
    
    */
    
    return {{2.*x[0], 0},{0., 1.}};
}

std::vector< double > test_fxn_3(std::vector< double > x){
    /*! test_fxn_3
    
    A test function used in test_numeric_gradient
    
    */
    
    return {x[1]*x[0]*x[0],x[1]+x[2]*x[0], x[0]};
}

std::vector< std::vector< double > > answer_fxn_3(std::vector< double > x){
    /*! answer_fxn_3
    
    A solution function used in test_numeric_gradient
    
    */
    
    return {{2*x[0]*x[1], x[2], 1.},
            {  x[0]*x[0],   1., 0.},
            {         0., x[0], 0.}};
}

std::vector< double > test_fxn_4(std::vector< double > x){
    /*! test_fxn_4
    
    A test function used in test_numeric_gradient
    
    */
    
    return {x[1]*x[0]*x[0],x[1]+x[2]*x[0]};
}

std::vector< std::vector< double > > answer_fxn_4(std::vector< double > x){
    /*! answer_fxn_4
    
    A solution function used in test_numeric_gradient
    
    */
    
    return {{2*x[0]*x[1], x[2]},
            {  x[0]*x[0],   1.},
            {         0., x[0]}};
}

int test_finite_difference(std::ofstream &results){
    /*!================================
    |    test_finite_difference    |
    ================================
    
    Test the finite difference method in the
    class FiniteDifference
    
    */
    
    //!Initialize test results
    int  test_num        = 3;
    bool test_results[test_num] = {false,false,false};
    
    //!Initialize the points about which to compute the finite differences
    std::vector< double > x01(1,2.6);
    std::vector< double > x02 = {2.3,-5.7};
    
    //!Initialize the finite difference
    finite_difference::FiniteDifference FD1 = finite_difference::FiniteDifference(test_fxn_1, 2, x01, 1e-6);
    finite_difference::FiniteDifference FD2 = finite_difference::FiniteDifference(test_fxn_2, 2, x02, 1e-6);
    
    //!Initialize the perturbation vectors
    std::vector< double > hvec1   = {1e-6};
    std::vector< double > hvec2a  = {1e-6,  0.};
    std::vector< double > hvec2b  = {  0.,1e-6};
    
    //!Compute the results of the finite difference
    std::vector< double > result1  = FD1.finite_difference(hvec1);
    std::vector< double > result2a = FD2.finite_difference(hvec2a);
    std::vector< double > result2b = FD2.finite_difference(hvec2b);
    
    //!Compute the answers
    std::vector< double > answer1 = answer_fxn_1(x01);
    std::vector< std::vector< double > > answer2 = answer_fxn_2(x02);
    
    test_results[0] = 1e-9>fabs(answer1[0]-result1[0]);
    test_results[1] = 1e-9>fabs(answer2[0][0]-result2a[0]);
    test_results[2] = 1e-9>fabs(answer2[1][1]-result2b[1]);
    
    //Compare all test results
    bool tot_result = true;
    for(int i = 0; i<test_num; i++){
        //std::cout << "\nSub-test " << i+1 << " result: " << test_results[i] << "\n";
        if(!test_results[i]){
            tot_result = false;
        }
    }
    
    if(tot_result){
        results << "test_finite_difference & True\\\\\n\\hline\n";
    }
    else{
        results << "test_finite_difference & False\\\\\n\\hline\n";
    }
    
    return 1;
}

int test_numeric_gradient(std::ofstream &results){
    /*!===============================
    |    test_numeric_gradient    |
    ===============================
    
    Test the numeric gradient method in the
    class FiniteDifference
    
    */
    
    //!Initialize test results
    int  test_num        = 3;
    bool test_results[test_num] = {false,false,false};
    
    //!Initialize the points about which to compute the finite differences
    std::vector< double > x01 = {2.3,-5.7};
    std::vector< double > x02 = {1.4,-2.0,3.4};
    std::vector< double > x03 = {1.4,-2.0,3.4};
    
    //!Initialize the finite difference
    finite_difference::FiniteDifference FD1 = finite_difference::FiniteDifference(test_fxn_2, 2, x01, 1e-6);
    finite_difference::FiniteDifference FD2 = finite_difference::FiniteDifference(test_fxn_3, 2, x02, 1e-6);
    finite_difference::FiniteDifference FD3 = finite_difference::FiniteDifference(test_fxn_4, 2, x03, 1e-6);
    
    //!Compute the results of the finite difference
    std::vector< std::vector< double > > result1  = FD1.numeric_gradient();
    std::vector< std::vector< double > > result2  = FD2.numeric_gradient();
    std::vector< std::vector< double > > result3  = FD3.numeric_gradient();
    
    //!Compute the answers
    std::vector< std::vector< double > > answer1 = answer_fxn_2(x01);
    std::vector< std::vector< double > > answer2 = answer_fxn_3(x02);
    std::vector< std::vector< double > > answer3 = answer_fxn_4(x03);
    
    //!Compare the resulting gradients
    test_results[0] = compare_vectors_of_vectors(answer1,result1);
    test_results[1] = compare_vectors_of_vectors(answer2,result2);
    test_results[2] = compare_vectors_of_vectors(answer3,result3);
    
    //Compare all test results
    bool tot_result = true;
    for(int i = 0; i<test_num; i++){
        //std::cout << "\nSub-test " << i+1 << " result: " << test_results[i] << "\n";
        if(!test_results[i]){
            tot_result = false;
        }
    }
    
    if(tot_result){
        results << "test_numeric_gradient & True\\\\\n\\hline\n";
    }
    else{
        results << "test_numeric_gradient & False\\\\\n\\hline\n";
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
    test_finite_difference(results);
    test_numeric_gradient(results);
    
    //Close the results file
    results.close();
}