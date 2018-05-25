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

template< int n_m, int m_m>
bool compare_vectors_matrix(std::vector<std::vector< double > > V, Eigen::Matrix<double,n_m,m_m> M, double tol=1e-5){
    /*!================================
    |    compare_vectors_matrix    |
    ================================
    
    Compare a vector of equal length vectors to a matrix
    
    */
    
    bool result = true;
    double abs_error;   //!The absolute error
    double rel_error;   //!The relative error
    
    if(V.size() != M.rows()){std::cout << "Error: vector of vectors and matrix do not have the same number of rows.\n";return false;}
    
    for(int i=0; i<M.rows(); i++){
        if(V[i].size()!=M.cols()){std::cout << "Error: vector of vectors and matrix do not have the same number of columns in row "<<i <<".";return false;}
        
        for(int j=0; j<M.cols(); j++){
            //std::cout << "M("<<i<<","<<j<<"): "<<M(i,j)<<" V["<<i<<"]["<<j<<"]: " << V[i][j] <<"\n";
            abs_error = fabs(M(i,j)-V[i][j]);
            rel_error = fabs(M(i,j)-V[i][j])/std::max(fabs(M(i,j)),fabs(V[i][j]));
            //std::cout << "abs error: " << abs_error << "\n";
            //std::cout << "rel error: " << rel_error << "\n";
            if((tol<abs_error)&&(tol<rel_error)){return false;}
        }
    }
    
    return true;
}

template< int n_m, int m_m>
bool compare_vectors_matrix_transpose(std::vector<std::vector< double > > V, Eigen::Matrix<double,n_m,m_m> M, double tol=1e-5){
    /*!================================
    |    compare_vectors_matrix    |
    ================================
    
    Compare a vector of equal length vectors to a matrix
    
    */
    
    bool result = true;
    double abs_error;   //!The absolute error
    double rel_error;   //!The relative error
    
    if(V.size() != M.cols()){std::cout << "Error: vector of vectors and matrix do not have the same number of rows.\n";return false;}
    
    for(int i=0; i<M.cols(); i++){
        if(V[i].size()!=M.rows()){std::cout << "Error: vector of vectors and matrix do not have the same number of columns in row "<<i <<".";return false;}
        
        for(int j=0; j<M.rows(); j++){
            //std::cout << "M("<<j<<","<<i<<"): "<<M(j,i)<<" V["<<i<<"]["<<j<<"]: " << V[i][j] <<"\n";
            abs_error = fabs(M(j,i)-V[i][j]);
            rel_error = fabs(M(j,i)-V[i][j])/std::max(fabs(M(j,i)),fabs(V[i][j]));
            //std::cout << "abs error: " << abs_error << "\n";
            //std::cout << "rel error: " << rel_error << "\n";
            if((tol<abs_error)&&(tol<rel_error)){return false;}
            
        }
    }
    return true;
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
    std::vector<bool> test_results(test_num, false);
    
    //!Initialize the floating point parameters
    double fpointer[19];
    Vector fparams = Vector_Xd_Map(fpointer,19,1);
    
    fparams[0] = 1000.;
    
    for(int i=1; i<19; i++){
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
    double rho    = fparams[ 0];             //!The density
    double lambda = fparams[ 1];             //!lambda micromorphic material parameter
    double mu     = fparams[ 2];             //!mu micromorphic material parameter
    double eta    = fparams[ 3];             //!eta micromorphic material parameter
    double tau    = fparams[ 4];             //!tau micromorphic material parameter
    double kappa  = fparams[ 5];             //!kappa micromorphic material parameter
    double nu     = fparams[ 6];             //!nu micromorphic material parameter
    double sigma  = fparams[ 7];             //!sigma micromorphic material parameter
    double tau1   = fparams[ 8];             //!tau1  micromorphic material parameter
    double tau2   = fparams[ 9];             //!tau2  micromorphic material parameter
    double tau3   = fparams[10];             //!tau3  micromorphic material parameter
    double tau4   = fparams[11];             //!tau4  micromorphic material parameter
    double tau5   = fparams[12];             //!tau5  micromorphic material parameter
    double tau6   = fparams[13];             //!tau6  micromorphic material parameter
    double tau7   = fparams[14];             //!tau7  micromorphic material parameter
    double tau8   = fparams[15];             //!tau8  micromorphic material parameter
    double tau9   = fparams[16];             //!tau9  micromorphic material parameter
    double tau10  = fparams[17];             //!tau10 micromorphic material parameter
    double tau11  = fparams[18];             //!tau11 micromorphic material parameter
    
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


std::vector<double> dPK2dC_parser(std::vector<double> C_in){
    /*!=======================
    |    dPK2dC_parser    |
    =======================
    
    A function which parses the PK2
    stress being varied by the 
    right Cauchy-Green deformation 
    tensor.
    
    Input:
    
        C_in:    The perturbed value of C.
    
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
    C(0,0) = C_in[0];
    C(0,1) = C_in[1];
    C(0,2) = C_in[2];
    C(1,0) = C_in[3];
    C(1,1) = C_in[4];
    C(1,2) = C_in[5];
    C(2,0) = C_in[6];
    C(2,1) = C_in[7];
    C(2,2) = C_in[8];
    
    //!Set the floating point parameters
    //!Initialize the floating point parameters
    double fpointer[19];
    Vector fparams = Vector_Xd_Map(fpointer,19,1);
    
    fparams[0] = 1000.;
    
    for(int i=1; i<19; i++){
        fparams[i] = i+1;
    }
    
    int ipointer[1];
    Vectori iparams = Vector_Xi_Map(ipointer,1,1);
    
    //!Compute and assign the stresses to the result measures
    micro_material::get_stress(fparams, iparams, C, Psi, Gamma, PK2_result, SIGMA_result, M_result);
    
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

std::vector<double> dSIGMAdC_parser(std::vector<double> C_in){
    /*!=========================
    |    dSIGMAdC_parser    |
    =========================
    
    A function which parses the symmetric
    stress being varied by the 
    right Cauchy-Green deformation 
    tensor.
    
    Input:
    
        C_in:    The perturbed value of C.
    
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
    
    //!Reset Psi
    C(0,0) = C_in[0];
    C(0,1) = C_in[1];
    C(0,2) = C_in[2];
    C(1,0) = C_in[3];
    C(1,1) = C_in[4];
    C(1,2) = C_in[5];
    C(2,0) = C_in[6];
    C(2,1) = C_in[7];
    C(2,2) = C_in[8];
    
    //!Set the floating point parameters
    double fpointer[19];
    Vector fparams = Vector_Xd_Map(fpointer,19,1);
    
    fparams[0] = 1000.;
    
    for(int i=1; i<19; i++){
        fparams[i] = i+1;
    }
    
    int ipointer[1];
    Vectori iparams = Vector_Xi_Map(ipointer,1,1);
    
    //!Compute and assign the stresses to the result measures
    micro_material::get_stress(fparams, iparams, C, Psi, Gamma, PK2_result, SIGMA_result, M_result);
    
    parsed_stress[0] = SIGMA_result(0,0);
    parsed_stress[1] = SIGMA_result(0,1);
    parsed_stress[2] = SIGMA_result(0,2);
    parsed_stress[3] = SIGMA_result(1,0);
    parsed_stress[4] = SIGMA_result(1,1);
    parsed_stress[5] = SIGMA_result(1,2);
    parsed_stress[6] = SIGMA_result(2,0);
    parsed_stress[7] = SIGMA_result(2,1);
    parsed_stress[8] = SIGMA_result(2,2);
    
    return parsed_stress;
}

std::vector<double> dPK2dPsi_parser(std::vector<double> Psi_in){
    /*!=========================
    |    dPK2dPsi_parser    |
    =========================
    
    A function which parses the PK2
    stress being varied by the 
    micro-deformation tensor Psi.
    
    Input:
    
        Psi_in:    The perturbed value of Psi.
    
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
    
    //!Reset Psi
    Psi(0,0) = Psi_in[0];
    Psi(0,1) = Psi_in[1];
    Psi(0,2) = Psi_in[2];
    Psi(1,0) = Psi_in[3];
    Psi(1,1) = Psi_in[4];
    Psi(1,2) = Psi_in[5];
    Psi(2,0) = Psi_in[6];
    Psi(2,1) = Psi_in[7];
    Psi(2,2) = Psi_in[8];
    
    //!Set the floating point parameters
    double fpointer[19];
    Vector fparams = Vector_Xd_Map(fpointer,19,1);
    
    fparams[0] = 1000.;
    
    for(int i=1; i<19; i++){
        fparams[i] = i+1;
    }
    
    int ipointer[1];
    Vectori iparams = Vector_Xi_Map(ipointer,1,1);
    
    //!Compute and assign the stresses to the result measures
    micro_material::get_stress(fparams, iparams, C, Psi, Gamma, PK2_result, SIGMA_result, M_result);
    
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

std::vector<double> dSIGMAdPsi_parser(std::vector<double> Psi_in){
    /*!===========================
    |    dSIGMAdPsi_parser    |
    ===========================
    
    A function which parses the symmetric
    stress being varied by the 
    micro-deformation tensor Psi.
    
    Input:
    
        Psi_in:    The perturbed value of Psi.
    
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
    Psi(0,0) = Psi_in[0];
    Psi(0,1) = Psi_in[1];
    Psi(0,2) = Psi_in[2];
    Psi(1,0) = Psi_in[3];
    Psi(1,1) = Psi_in[4];
    Psi(1,2) = Psi_in[5];
    Psi(2,0) = Psi_in[6];
    Psi(2,1) = Psi_in[7];
    Psi(2,2) = Psi_in[8];
    
    //!Set the floating point parameters
    double fpointer[19];
    Vector fparams = Vector_Xd_Map(fpointer,19,1);
    
    fparams[0] = 1000.;
    
    for(int i=1; i<19; i++){
        fparams[i] = i+1;
    }
    
    int ipointer[1];
    Vectori iparams = Vector_Xi_Map(ipointer,1,1);
    
    //!Compute and assign the stresses to the result measures
    micro_material::get_stress(fparams, iparams, C, Psi, Gamma, PK2_result, SIGMA_result, M_result);
    
    parsed_stress[0] = SIGMA_result(0,0);
    parsed_stress[1] = SIGMA_result(0,1);
    parsed_stress[2] = SIGMA_result(0,2);
    parsed_stress[3] = SIGMA_result(1,0);
    parsed_stress[4] = SIGMA_result(1,1);
    parsed_stress[5] = SIGMA_result(1,2);
    parsed_stress[6] = SIGMA_result(2,0);
    parsed_stress[7] = SIGMA_result(2,1);
    parsed_stress[8] = SIGMA_result(2,2);
    
    return parsed_stress;
}

std::vector<double> dPK2dGamma_parser(std::vector<double> Gamma_in){
    /*!===========================
    |    dPK2dGamma_parser    |
    ===========================
    
    A function which parses the PK2
    stress being varied by the 
    micro-deformation gradient tensor Gamma.
    
    Input:
    
        Psi_in:    The perturbed value of Psi.
    
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
    
    //!Reset Gamma
    int temp_indx = 0;
    for(int K=0; K<3; K++){
        for(int I=0; I<3; I++){
            for(int J=0; J<3; J++){
                Gamma(I,J,K) = Gamma_in[temp_indx];
                temp_indx++;
            }
        }
    }
    
    //!Set the floating point parameters
    double fpointer[19];
    Vector fparams = Vector_Xd_Map(fpointer,19,1);
    
    fparams[0] = 1000.;
    
    for(int i=1; i<19; i++){
        fparams[i] = i+1;
    }
    
    int ipointer[1];
    Vectori iparams = Vector_Xi_Map(ipointer,1,1);
    
    //!Compute and assign the stresses to the result measures
    micro_material::get_stress(fparams, iparams, C, Psi, Gamma, PK2_result, SIGMA_result, M_result);
    
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

std::vector<double> dSIGMAdGamma_parser(std::vector<double> Gamma_in){
    /*!=========================
    |    dSIGMAdGamma_parser    |
    =========================
    
    A function which parses the symmetric
    stress being varied by the 
    micro-deformation gradient tensor Gamma.
    
    Input:
    
        Psi_in:    The perturbed value of Psi.
    
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
    
    //!Reset Gamma
    int temp_indx = 0;
    for(int K=0; K<3; K++){
        for(int I=0; I<3; I++){
            for(int J=0; J<3; J++){
                Gamma(I,J,K) = Gamma_in[temp_indx];
                temp_indx++;
            }
        }
    }
    
    //!Set the floating point parameters
    double fpointer[19];
    Vector fparams = Vector_Xd_Map(fpointer,19,1);
    
    fparams[0] = 1000.;
    
    for(int i=1; i<19; i++){
        fparams[i] = i+1;
    }
    
    int ipointer[1];
    Vectori iparams = Vector_Xi_Map(ipointer,1,1);
    
    //!Compute and assign the stresses to the result measures
    micro_material::get_stress(fparams, iparams, C, Psi, Gamma, PK2_result, SIGMA_result, M_result);
    
    parsed_stress[0] = SIGMA_result(0,0);
    parsed_stress[1] = SIGMA_result(0,1);
    parsed_stress[2] = SIGMA_result(0,2);
    parsed_stress[3] = SIGMA_result(1,0);
    parsed_stress[4] = SIGMA_result(1,1);
    parsed_stress[5] = SIGMA_result(1,2);
    parsed_stress[6] = SIGMA_result(2,0);
    parsed_stress[7] = SIGMA_result(2,1);
    parsed_stress[8] = SIGMA_result(2,2);
    
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
    int  test_num        = 10;
    std::vector<bool> test_results(test_num,false);
    
    //!Seed the random number generator
    srand (1);
    
    //!Initialize the common tensors
    tensor::Tensor23 ITEN = tensor::eye();
    
    //!Initialize the floating point parameters
    double fpointer[19];
    Vector fparams = Vector_Xd_Map(fpointer,19,1);
    
    fparams[0] = 1000.;
    
    for(int i=1; i<19; i++){
        fparams[i] = i+1;
    }
    
    int ipointer[1];
    
    Vectori iparams = Vector_Xi_Map(ipointer,1,1);
    
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
    micro_material::get_stress(fparams, iparams, C, Psi, Gamma, PK2_result, SIGMA_result, M_result);
    
    test_results[0] = PK2_answer.data.isApprox(PK2_result.data);
    test_results[1] = SIGMA_answer.data.isApprox(SIGMA_result.data);
    test_results[2] = M_answer.data.isApprox(M_result.data);
    
    
    //!Compute the tangents and compare them to the numeric results
    
    //!Initialize the tangents
    tensor::Tensor43 dPK2dC = tensor::Tensor43({3,3,3,3});
    tensor::Tensor43 dPK2dPsi = tensor::Tensor43({3,3,3,3});
    tensor::Tensor53 dPK2dGamma = tensor::Tensor53({3,3,3,3,3});
    
    tensor::Tensor43 dSIGMAdC = tensor::Tensor43({3,3,3,3});
    tensor::Tensor43 dSIGMAdPsi = tensor::Tensor43({3,3,3,3});
    tensor::Tensor53 dSIGMAdGamma = tensor::Tensor53({3,3,3,3,3});
    
    tensor::Tensor53 dMdC = tensor::Tensor53({3,3,3,3,3});
    tensor::Tensor53 dMdPsi = tensor::Tensor53({3,3,3,3,3});
    tensor::Tensor63 dMdGamma = tensor::Tensor63({3,3,3,3,3,3});
    
    micro_material::get_stress( fparams,    iparams,            C, Psi, Gamma, PK2_result, SIGMA_result, M_result,
                                 dPK2dC,   dPK2dPsi,   dPK2dGamma,
                               dSIGMAdC, dSIGMAdPsi, dSIGMAdGamma,
                                   dMdC,     dMdPsi,     dMdGamma);
    
    //!Compute the numeric results for dPK2dC
    std::vector<double> vec_ten;
    vec_ten.resize(9);
    vec_ten[0] = C(0,0);
    vec_ten[1] = C(0,1);
    vec_ten[2] = C(0,2);
    vec_ten[3] = C(1,0);
    vec_ten[4] = C(1,1);
    vec_ten[5] = C(1,2);
    vec_ten[6] = C(2,0);
    vec_ten[7] = C(2,1);
    vec_ten[8] = C(2,2);
    finite_difference::FiniteDifference FD(dPK2dC_parser,2,vec_ten,1e-6);
    std::vector< std::vector< double > > temp_gradient = FD.numeric_gradient();
    
    test_results[3]=compare_vectors_matrix_transpose(temp_gradient,dPK2dC.data);
    
    //print_vector_of_vectors("dPK2dC_answer",temp_gradient);
    //std::cout << "dPK2dC_result\n" << dPK2dC.data << "\n";
    
    //!Compute the numeric results for dSIGMAdC
    FD = finite_difference::FiniteDifference(dSIGMAdC_parser,2,vec_ten,1e-6);
    temp_gradient = FD.numeric_gradient();
    
    //print_vector_of_vectors("dSIGMAdC_answer",temp_gradient);
    //std::cout << "dSIGMAdC_result\n" << dSIGMAdC.data << "\n";
    
    test_results[4] = compare_vectors_matrix_transpose(temp_gradient,dSIGMAdC.data);
    
    vec_ten[0] = Psi(0,0);
    vec_ten[1] = Psi(0,1);
    vec_ten[2] = Psi(0,2);
    vec_ten[3] = Psi(1,0);
    vec_ten[4] = Psi(1,1);
    vec_ten[5] = Psi(1,2);
    vec_ten[6] = Psi(2,0);
    vec_ten[7] = Psi(2,1);
    vec_ten[8] = Psi(2,2);
    FD = finite_difference::FiniteDifference(dPK2dPsi_parser,2,vec_ten,1e-6);
    temp_gradient = FD.numeric_gradient();
    
    //print_vector_of_vectors("dPK2dPsi_answer",temp_gradient);
    //std::cout << "dPK2dPsi_result\n" << dPK2dPsi.data << "\n";
    
    test_results[5] = compare_vectors_matrix_transpose(temp_gradient,dPK2dPsi.data);
    
    FD = finite_difference::FiniteDifference(dSIGMAdPsi_parser,2,vec_ten,1e-6);
    temp_gradient = FD.numeric_gradient();
    
    //print_vector_of_vectors("dSIGMAdPsi_answer",temp_gradient);
    //std::cout << "dSIGMAdPsi_result\n" << dSIGMAdPsi.data << "\n";
    
    test_results[6] = compare_vectors_matrix_transpose(temp_gradient,dSIGMAdPsi.data);
    
    int temp_indx=0;
    
    vec_ten.resize(27);
    for(int K=0; K<3; K++){
        for(int I=0; I<3; I++){
            for(int J=0; J<3; J++){
                vec_ten[temp_indx] = Gamma(I,J,K);
                temp_indx++;
            }
        }
    }
    
    FD = finite_difference::FiniteDifference(dPK2dGamma_parser,2,vec_ten,1e-6);
    temp_gradient = FD.numeric_gradient();
    
    //Rearrange the gradient into a form we can compare easily
    std::vector< std::vector< double > > rearranged_gradient;
    rearranged_gradient.resize(27);
    for(int I=0; I<rearranged_gradient.size(); I++){rearranged_gradient[I].resize(9);}
    
    for(int N=0; N<9; N++){
        for(int I=0; I<3; I++){
            for(int J=0; J<3; J++){
                for(int K=0; K<3; K++){
                    rearranged_gradient[K+9*I+3*J][N] = temp_gradient[9*K+3*I+J][N];
                }
            }
        }
    }
    //print_vector_of_vectors("dPK2dGamma_answer",rearranged_gradient);
    //std::cout << "dPK2dGamma_result\n" << dPK2dGamma.data << "\n";
    
    test_results[7] = compare_vectors_matrix_transpose(rearranged_gradient,dPK2dGamma.data);
    
    FD = finite_difference::FiniteDifference(dSIGMAdGamma_parser,2,vec_ten,1e-6);
    temp_gradient = FD.numeric_gradient();
    
    //Rearrange the gradient into a form we can compare easily
    
    for(int N=0; N<9; N++){
        for(int I=0; I<3; I++){
            for(int J=0; J<3; J++){
                for(int K=0; K<3; K++){
                    rearranged_gradient[K+9*I+3*J][N] = temp_gradient[9*K+3*I+J][N];
                }
            }
        }
    }
    //print_vector_of_vectors("dPK2dGamma_answer",rearranged_gradient);
    //std::cout << "dPK2dGamma_result\n" << dPK2dGamma.data << "\n";
    
    test_results[8] = compare_vectors_matrix_transpose(rearranged_gradient,dSIGMAdGamma.data);
    
    test_results[9] = dMdGamma.data.isApprox(C_stiffness.data);
    
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

