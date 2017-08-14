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
//#include <cstdarg>
#include <Eigen/Dense>
#include <tensor.h>

int test_tensor_functionality(std::ofstream &results){
    /*!====================================
    |       test_tensor_functionality   |
    =====================================
    
    A test of some of the basic tensor functionality. 
    This should show that tensors of different orders 
    can be generated and manipulated.*/

    std::vector< int > v_shape;
    v_shape.resize(1);
    v_shape[0] = 6;
    
    std::cout << "Forming tensor\n";
    tensor::Tensor V = tensor::Tensor(v_shape);
    
    std::cout << "V is of size " << V.data.rows() << "x" << V.data.cols() << "\n";
    
    std::cout << "V.data:\n" << V.data << "\n";
    
    V(2) = -7;
    
    std::cout << "V.data:\n" << V.data << "\n";
    
    std::vector < int > m_shape;
    m_shape.resize(2);
    
    m_shape[0] = 3;
    m_shape[1] = 3;
    
    tensor::Tensor M = tensor::Tensor(m_shape);
    
    std::cout << "M is of size " << M.data.rows() << "x" << M.data.cols() << "\n";
    
    std::cout << "M.data:\n" << M.data << "\n";
    
    M(2,1) = -8;
    
    std::cout << "M.data:\n" << M.data << "\n";
    
    std::vector < int > t_shape;
    t_shape.resize(3);
    
    t_shape[0] = 3;
    t_shape[1] = 3;
    t_shape[2] = 7;
    
    tensor::Tensor T = tensor::Tensor(t_shape);
}

int main(){
    /*!==========================
    |         main            |
    ===========================
    
    The main loop which runs the tests defined in the 
    accompanying functions. Each function should output
    the function name followed by & followed by True or 
    False if the test passes or fails respectively.*/
    
    std::cout << "Compiled successfully";
    
    std::ofstream results;
    //Open the results file
    results.open ("results.tex");
    results << "Writing this to a file.\n";
        
    //!Run the test functions
    test_tensor_functionality(results);
    results.close();
}

