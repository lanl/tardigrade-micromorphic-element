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
  | finite_difference: An implementation of the         |
  |                    computation of the numeric       |
  |                    gradient via finite differences. |
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
        fparams[i] = 0.1*(i+1);
    }
    
    //!Initialize the common tensors
    tensor::Tensor23 I = tensor::eye();
    
    //!initialize the tensor answers
    tensor::Tensor43 A_answer({3,3,3,3});
    tensor::Tensor43 B_answer({3,3,3,3});
    tensor::Tensor63 C_answer({3,3,3,3,3,3});
    tensor::Tensor43 D_answer({3,3,3,3});
    
    //!initialize the tensor results
    tensor::Tensor43 A_result({3,3,3,3});
    tensor::Tensor43 B_result({3,3,3,3});
    tensor::Tensor63 C_result({3,3,3,3,3,3});
    tensor::Tensor43 D_result({3,3,3,3});
    
    //!Compute the expected tensors
    micro_material::generate_A_stiffness(fparams,A_result);
    micro_material::generate_B_stiffness(fparams,B_result);
    micro_material::generate_C_stiffness(fparams,C_result);
    micro_material::generate_D_stiffness(fparams,D_result);
    
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

void populate_deformation_measures(tensor::Tensor23& C, tensor::Tensor23& Psi, tensor::Tensor33& Gamma){
    /*!=======================================
    |    populate_deformation_measures    |
    =======================================
    
    Populate the deformation measures in a consistent 
    way.
    
    */
    
    //!Initialize the fundamental deformation measures
    tensor::Tensor F({3,3});
    tensor::Tensor chi({3,3});
    tensor::Tensor grad_chi({3,3,3});
    
    //!Populate the gradients
    std::vector< double > position = {1.2, 2.4, -0.42};
    std::vector< std::vector< double > > gradients = compute_gradients(position); //Compute the numeric gradients
    
    for(int i=0; i<3; i++){
        //!Populate the deformation gradient
        F(0,i)          = gradients[i][ 0];
        F(1,i)          = gradients[i][ 1];
        F(2,i)          = gradients[i][ 2];
        //!Populate the gradient of chi
        grad_chi(0,0,i) = gradients[i][ 3];
        grad_chi(1,1,i) = gradients[i][ 4];
        grad_chi(2,2,i) = gradients[i][ 5];
        grad_chi(1,2,i) = gradients[i][ 6];
        grad_chi(0,2,i) = gradients[i][ 7];
        grad_chi(0,1,i) = gradients[i][ 8];
        grad_chi(2,1,i) = gradients[i][ 9];
        grad_chi(2,0,i) = gradients[i][10];
        grad_chi(1,0,i) = gradients[i][11];
    }
    
    //Because we are specifying the deformation, and 
    //not the actual current coordinates, we have to 
    //add 1 to the diagonal terms
    F(0,0) += 1;
    F(1,1) += 1;
    F(2,2) += 1;
    
    //!Populate the chi
    std::vector< double > U = test_deformation(position);
    chi(0,0) = 1.+U[ 3];
    chi(1,1) = 1.+U[ 4];
    chi(2,2) = 1.+U[ 5];
    chi(1,2) =    U[ 6];
    chi(0,2) =    U[ 7];
    chi(0,1) =    U[ 8];
    chi(2,1) =    U[ 9];
    chi(2,0) =    U[10];
    chi(1,0) =    U[11];
    
    for(int I=0; I<3; I++){
        for(int J=0; J<3; J++){
            for(int i=0; i<3; i++){
                C(I,J)   += F(i,I)*F(i,J);
                Psi(I,J) += F(i,I)*chi(i,J);
            }
            for(int K=0; K<3; K++){
                for(int i=0; i<3; i++){
                    Gamma(I,J,K) += F(i,I)*grad_chi(i,J,K);
                }
            }
        }
    }
    return;
}


std::vector<double> dStressddef_parser(std::vector<double> C0){
    /*!=======================
    |    dPK2dC_parser    |
    =======================
    
    A function which parses the PK2
    stress being varied by the 
    right Cauchy-Green deformation 
    tensor.
    
    Input:
    
        C0:            Initial point to compute the 
                       gradient around.
    
    */
    
    std::vector<double> parsed_stress; //!The stress after being parsed
    parsed_stress.resize(9);
    
    //!Populate derived deformation measures
    tensor::Tensor23 C({3,3});          //!The right Cauchy-Green deformation tensor (will be discarded)
    tensor::Tensor23 Psi({3,3});        //!The micro-deformation measure
    tensor::Tensor33 Gamma({3,3,3});    //!The higher order micro-deformation measure
    
    populate_deformation_measures(C,Psi,Gamma);
    
    //!The initialization of the result stress tensors
    tensor::Tensor23 PK2_result({3,3});
    tensor::Tensor23 SIGMA_result({3,3});
    tensor::Tensor33 M_result({3,3,3});
    
    //!Reset C
    C(0,0) = C0[0];
    C(0,1) = C0[1];
    C(0,2) = C0[2];
    C(1,0) = C0[3];
    C(1,1) = C0[4];
    C(1,2) = C0[5];
    C(2,0) = C0[6];
    C(2,1) = C0[7];
    C(2,2) = C0[8];
    
    //!Compute and assign the stresses to the result measures
    micro_material::get_stress(fparams, {}, C0, Psi, Gamma, PK2_result, SIGMA_result, M_result);
    
    parsed_stress[0] = PK2_result(0,0);
    parsed_stress[1] = PK2_result(0,1);
    parsed_stress[2] = PK2_result(0,2);
    parsed_stress[3] = PK2_result(1,0);
    parsed_stress[4] = PK2_result(1,1);
    parsed_stress[5] = PK2_result(1,2);
    parsed_stress[6] = PK2_result(2,0);
    parsed_stress[7] = PK2_result(2,1);
    parsed_stress[8] = PK2_result(2,2);
    
    return parsed_stress;
}

void test_get_stress(std::ofstream &results){
    /*!=========================
    |    test_get_stress    |
    =========================
    
    A test for the computation of the 
    stress tensors. They are compared to 
    another formulation of the tensors 
    to ensure that they are consistent.
    
    */
    
    //!Initialize test results
    int  test_num        = 3;
    bool test_results[test_num] = {false,false,false};
    
    //!Seed the random number generator
    srand (1);
    
    //!Initialize the common tensors
    tensor::Tensor23 ITEN = tensor::eye();
    
    //!Initialize the floating point parameters
    std::vector< double > fparams(18,0.);
    
    for(int i=0; i<18; i++){
        fparams[i] = i+1;
    }
    
    //!Populate derived deformation measures
    tensor::Tensor23 C({3,3});          //!The right Cauchy-Green deformation tensor
    tensor::Tensor23 Psi({3,3});        //!The micro-deformation measure
    tensor::Tensor33 Gamma({3,3,3});    //!The higher order micro-deformation measure
    
    populate_deformation_measures(C,Psi,Gamma);
    
    tensor::Tensor23 Cinv = C.inverse();
    
    //!Compute the strain measures
    tensor::Tensor23 macro_E = 0.5*(C-ITEN);
    tensor::Tensor23 micro_E = Psi-ITEN;
    
    //!Initialize the stiffness tensors
    tensor::Tensor43 A_stiffness({3,3,3,3});
    tensor::Tensor43 B_stiffness({3,3,3,3});
    tensor::Tensor63 C_stiffness({3,3,3,3,3,3});
    tensor::Tensor43 D_stiffness({3,3,3,3});
        
    //!Compute the stiffness tensors
    micro_material::generate_A_stiffness(fparams,A_stiffness);
    micro_material::generate_B_stiffness(fparams,B_stiffness);
    micro_material::generate_C_stiffness(fparams,C_stiffness);
    micro_material::generate_D_stiffness(fparams,D_stiffness);
    
    tensor::Tensor23 PK2_answer({3,3});   //!The expected second Piola-Kirchhoff Stress
    tensor::Tensor23 SIGMA_answer({3,3}); //!The symmetric micro-stress
    tensor::Tensor33 M_answer({3,3,3});   //!The expected higher order stress
    
    //!Compute the answer stress tensors
    for(int I=0; I<3; I++){
        for(int J=0; J<3; J++){
            for(int K=0; K<3; K++){
                for(int L=0; L<3; L++){
                    PK2_answer(I,J)   += A_stiffness(I,J,K,L)*macro_E(K,L) + D_stiffness(I,J,K,L)*micro_E(K,L);
                    SIGMA_answer(I,J) += A_stiffness(I,J,K,L)*macro_E(K,L) + D_stiffness(I,J,K,L)*micro_E(K,L);
                }
            }
            
            for(int K=0; K<3; K++){
                for(int L=0; L<3; L++){
                    for(int Q=0; Q<3; Q++){
                        for(int R=0; R<3; R++){
                            PK2_answer(I,J)   +=   (B_stiffness(I,Q,K,L)*micro_E(K,L) + D_stiffness(I,Q,K,L)*macro_E(K,L))*(micro_E(R,Q) + ITEN(R,Q))*Cinv(J,R);
                            SIGMA_answer(I,J) +=   (B_stiffness(I,Q,K,L)*micro_E(K,L) + D_stiffness(I,Q,K,L)*macro_E(K,L))*(micro_E(R,Q) + ITEN(R,Q))*Cinv(J,R)
                                                 + (B_stiffness(J,Q,K,L)*micro_E(K,L) + D_stiffness(J,Q,K,L)*macro_E(K,L))*(micro_E(R,Q) + ITEN(R,Q))*Cinv(I,R);
                        }
                    }
                }
            }
            
            for(int Q=0; Q<3; Q++){
                for(int R=0; R<3; R++){
                    for(int L=0; L<3; L++){
                        for(int M=0; M<3; M++){
                            for(int N=0; N<3; N++){
                                for(int S=0; S<3; S++){
                                    PK2_answer(I,J)   +=   C_stiffness(I,Q,R,L,M,N)*Gamma(L,M,N)*Cinv(S,J)*Gamma(S,Q,R);
                                    SIGMA_answer(I,J) +=   C_stiffness(I,Q,R,L,M,N)*Gamma(L,M,N)*Cinv(S,J)*Gamma(S,Q,R)
                                                         + C_stiffness(J,Q,R,L,M,N)*Gamma(L,M,N)*Cinv(S,I)*Gamma(S,Q,R);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    for(int I=0; I<3; I++){
        for(int J=0; J<3; J++){
            for(int K=0; K<3; K++){
                for(int L=0; L<3; L++){
                    for(int M=0; M<3; M++){
                        for(int N=0; N<3; N++){
                            M_answer(I,J,K) += C_stiffness(I,J,K,L,M,N)*Gamma(L,M,N);
                        }
                    }
                }
            }
        }
    }
    
    //!The initialization of the result stress tensors
    tensor::Tensor23 PK2_result({3,3});
    tensor::Tensor23 SIGMA_result({3,3});
    tensor::Tensor33 M_result({3,3,3});
    
    //!Compute and assign the stresses to the result measures
    micro_material::get_stress(fparams, {}, C, Psi, Gamma, PK2_result, SIGMA_result, M_result);
    
    test_results[0] = PK2_answer.data.isApprox(PK2_result.data);
    test_results[1] = SIGMA_answer.data.isApprox(SIGMA_result.data);
    test_results[2] = M_answer.data.isApprox(M_result.data);
    
    //Compare all test results
    bool tot_result = true;
    for(int i = 0; i<test_num; i++){
        //std::cout << "\nSub-test " << i+1 << " result: " << test_results[i] << "\n";
        if(!test_results[i]){
            tot_result = false;
        }
    }
    
    if(tot_result){
        results << "test_get_stress & True\\\\\n\\hline\n";
    }
    else{
        results << "test_get_stress & False\\\\\n\\hline\n";
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
    test_get_stress(results);
    
    //Close the results file
    results.close();
}

