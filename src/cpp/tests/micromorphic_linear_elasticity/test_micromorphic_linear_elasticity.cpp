/*!=======================================================
  |                                                     |
  |      test_micromorphic_linear_elasticity.cpp        |
  |                                                     |
  -------------------------------------------------------
  | The unit test file for                              |
  | micromorphic_linear_elasticity.h/cpp. This file     |
  | tests the classes and functions defined in          |
  | micromorphic_linear_elasticity.h/cpp.               |
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
#include <micromorphic_linear_elasticity.h>
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

void test_stiffness_tensors(std::ofstream &results){
    /*!================================
    |    test_stiffness_tensors    |
    ================================
    
    Test the computation of the different 
    stiffness tensors to make sure they 
    are consistent with the expected value.
    
    */
    
    //!Initialize test results
    int  test_num        = 4;
    bool test_results[test_num] = {false,false,false,false};
    
    //!Initialize the floating point parameters
    std::vector< double > fparams(18,0.);
    
    for(int i=0; i<18; i++){
        fparams[i] = i+1;
    }
    
    //!Initialize the common tensors
    tensor::Tensor I = tensor::eye();
    
    //!initialize the tensor answers
    tensor::Tensor A_answer({3,3,3,3});
    tensor::Tensor B_answer({3,3,3,3});
    tensor::Tensor C_answer({3,3,3,3,3,3});
    tensor::Tensor D_answer({3,3,3,3});
    
    //!Compute the expected tensors
    tensor::Tensor A_result = micro_material::generate_A_stiffness(fparams);
    tensor::Tensor B_result = micro_material::generate_B_stiffness(fparams);
    tensor::Tensor C_result = micro_material::generate_C_stiffness(fparams);
    tensor::Tensor D_result = micro_material::generate_D_stiffness(fparams);
    
    //!Extract the values of fparams
    double lambda = fparams[ 0];             //!lambda micromorphic material parameter
    double mu     = fparams[ 1];             //!mu micromorphic material parameter
    double eta    = fparams[ 2];             //!eta micromorphic material parameter
    double tau    = fparams[ 3];             //!tau micromorphic material parameter
    double kappa  = fparams[ 4];             //!kappa micromorphic material parameter
    double nu     = fparams[ 5];             //!nu micromorphic material parameter
    double sigma  = fparams[ 6];             //!sigma micromorphic material parameter
    double tau1   = fparams[ 7];             //!tau1  micromorphic material parameter
    double tau2   = fparams[ 8];             //!tau2  micromorphic material parameter
    double tau3   = fparams[ 9];             //!tau3  micromorphic material parameter
    double tau4   = fparams[10];             //!tau4  micromorphic material parameter
    double tau5   = fparams[11];             //!tau5  micromorphic material parameter
    double tau6   = fparams[12];             //!tau6  micromorphic material parameter
    double tau7   = fparams[13];             //!tau7  micromorphic material parameter
    double tau8   = fparams[14];             //!tau8  micromorphic material parameter
    double tau9   = fparams[15];             //!tau9  micromorphic material parameter
    double tau10  = fparams[16];             //!tau10 micromorphic material parameter
    double tau11  = fparams[17];             //!tau11 micromorphic material parameter
    
    //!Compute the expected value of the A stiffness tensor
    for(int K=0; K<3; K++){
        for(int L=0; L<3; L++){
            for(int M=0; M<3; M++){
                for(int N=0; N<3; N++){
                    A_answer(K,L,M,N) = lambda*I(K,L)*I(M,N)+mu*(I(K,M)*I(L,N)+I(K,N)*I(L,M));
                }
            }
        }
    }
    
    //!Compute the expected value of the B stiffness tensor
    for(int K=0; K<3; K++){
        for(int L=0; L<3; L++){
            for(int M=0; M<3; M++){
                for(int N=0; N<3; N++){
                    B_answer(K,L,M,N) = (eta-tau)*I(K,L)*I(M,N) + kappa*I(K,M)*I(L,N) + nu*I(K,N)*I(L,M)
                                        - sigma*(I(K,M)*I(L,N) + I(K,N)*I(L,M));
                }
            }
        }
    }
    
    //!Compute the expected value of the C stiffness tensor
    for(int K=0; K<3; K++){
        for(int L=0; L<3; L++){
            for(int M=0; M<3; M++){
                for(int N=0; N<3; N++){
                    for(int P=0; P<3; P++){
                        for(int Q=0; Q<3; Q++){
                            C_answer(K,L,M,N,P,Q) =    tau1*(I(K,L)*I(M,N)*I(P,Q) + I(K,Q)*I(L,M)*I(N,P))
                                                    +  tau2*(I(K,L)*I(M,P)*I(N,Q) + I(K,M)*I(L,Q)*I(N,P))
                                                    +  tau3*I(K,L)*I(M,Q)*I(N,P)  + tau4*I(K,N)*I(L,M)*I(P,Q)
                                                    +  tau5*(I(K,M)*I(L,N)*I(P,Q) + I(K,P)*I(L,M)*I(N,Q))
                                                    +  tau6*I(K,M)*I(L,P)*I(N,Q)  + tau7*(I(K,N)*I(L,P)*I(M,Q))
                                                    +  tau8*(I(K,P)*I(L,Q)*I(M,N) + I(K,Q)*I(L,N)*I(M,P)) + tau9*I(K,N)*I(L,Q)*I(M,P)
                                                    + tau10*I(K,P)*I(L,N)*I(M,Q)  + tau11*I(K,Q)*I(L,P)*I(M,N);
                        }
                    }
                }
            }
        }
    }
    
    //!Compute the expected value of the D stiffness tensor
    for(int K=0; K<3; K++){
        for(int L=0; L<3; L++){
            for(int M=0; M<3; M++){
                for(int N=0; N<3; N++){
                    D_answer(K,L,M,N) = tau*I(K,L)*I(M,N) + sigma*(I(K,M)*I(L,N) + I(K,N)*I(L,M));
                }
            }
        }
    }
    
    //!Compare the answers to the results
    test_results[0] = A_answer.data.isApprox(A_result.data);
    test_results[1] = B_answer.data.isApprox(B_result.data);
    test_results[2] = C_answer.data.isApprox(C_result.data);
    test_results[3] = D_answer.data.isApprox(D_result.data);
    
    //Compare all test results
    bool tot_result = true;
    for(int i = 0; i<test_num; i++){
        //std::cout << "\nSub-test " << i+1 << " result: " << test_results[i] << "\n";
        if(!test_results[i]){
            tot_result = false;
        }
    }
    
    if(tot_result){
        results << "test_stiffness_tensors & True\\\\\n\\hline\n";
    }
    else{
        results << "test_stiffness_tensors & False\\\\\n\\hline\n";
    }
    
    return;
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
    test_stiffness_tensors(results);
    
    
    //Close the results file
    results.close();
}

