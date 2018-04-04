/*!============================================================
  |                                                           |
  |         micromorphic_linear_elasticity_voigt.cpp          |
  |                                                           |
  -------------------------------------------------------------
  | The source file for the definition of a                   |
  | micromorphic linear elasticity using voigt notation.      |
  -------------------------------------------------------------
  | Notes: Micromorphic constitutive models should be         |
  |        developed in the namespace micro_material          |
  |        and have the function get_stress. This             |
  |        function should read in the right Cauchy           |
  |        green deformation tensor, Psi, Gamma, and          |
  |        write the PK2 stress, the symmetric stress         |
  |        in the reference configuration, and the            |
  |        higher order couple stress in the reference        |
  |        configuration. (ADDITIONAL VALUES WILL BE          |
  |        ADDED OVER TIME).                                  |
  =============================================================
  | Dependencies:                                             |
  | Eigen: Open source matrix library available at            |
  |        eigen.tuxfamily.org.                               |
  =============================================================*/

#include <iostream>  
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <deformation_measures.h>
#include <micromorphic_linear_elasticity_voigt.h>

namespace micro_material{

    void get_stress(const double t,      const double dt,       const double (&params)[18],
                    const Matrix_3x3 &F, const Matrix_3x3 &chi, const Matrix_3x9 &grad_chi,
                    std::vector<double> &SDVS,    Vector_9 &PK2, Vector_9 &SIGMA, Vector_27 &M){
    
        //Extract the parameters
        double lambda  = params[ 0];
        double mu      = params[ 1];
        double eta     = params[ 2];
        double tau     = params[ 3];
        double kappa   = params[ 4];
        double nu      = params[ 5];
        double sigma   = params[ 6];
        double tau1    = params[ 7];
        double tau2    = params[ 8];
        double tau3    = params[ 9];
        double tau4    = params[10];
        double tau5    = params[11];
        double tau6    = params[12];
        double tau7    = params[13];
        double tau8    = params[14];
        double tau9    = params[15];
        double tau10   = params[16];
        double tau11   = params[17];

        //Initialize the stiffness matrices
        SpMat A( 9, 9);
        SpMat B( 9, 9);
        SpMat C(27,27);
        SpMat D( 9, 9);

        //Populate the stiffness matrices
        compute_A_voigt(lambda,mu,A);
        compute_B_voigt(eta,kappa,nu,sigma,tau,B);
        compute_C_voigt(tau1,tau2,tau3, tau4, tau5,tau6,
                        tau7,tau8,tau9,tau10,tau11,C);
        compute_D_voigt(sigma,tau,D);

        //Compute the deformation measures
        Matrix_3x3 RCG; //The right cauchy green deformation tensor
        deformation_measures::get_right_cauchy_green(F,RCG);
        Matrix_3x3 RCGinv = RCG.inverse(); //The inverse of the right cauchy green deformation tensor
        Matrix_3x3 Psi; //The second order micromorphic deformation measure
        deformation_measures::get_psi(F,chi,Psi);
        Matrix_3x9 Gamma; //The higher order micromorphic deformation measure
        deformation_measures::get_gamma(F,grad_chi,Gamma);

        //Compute the strain measures
        Matrix_3x3 E;
        Matrix_3x3 E_micro;
        deformation_measures::get_lagrange_strain(F,E);
        deformation_measures::get_micro_strain(Psi,E_micro);
        
        //std::cout << "F:\n" << F << "\n";
        //std::cout << "chi:\n" << chi << "\n";
        //std::cout << "grad chi:\n" << grad_chi << "\n";
        //std::cout << "Psi:\n" << Psi << "\n";
        //std::cout << "Gamma:\n" << Gamma << "\n";
        //std::cout << "E:\n" << E << "\n";
        //std::cout << "E_micro:\n" << E_micro << "\n";

        //Put the strain measures in voigt notation
        Vector_9  E_voigt;
        Vector_9  E_micro_voigt;
        Vector_27 Gamma_voigt;

        deformation_measures::voigt_3x3_tensor(E,       E_voigt);
        deformation_measures::voigt_3x3_tensor(E_micro, E_micro_voigt);
        deformation_measures::voigt_3x9_tensor(Gamma,   Gamma_voigt);

        //Compute the stress measures
        compute_PK2_stress(E_voigt, E_micro_voigt, Gamma_voigt, RCGinv, Psi, Gamma,
                           A,       B,             C,           D,      PK2);
                           
        compute_symmetric_stress(E_voigt, E_micro_voigt, Gamma_voigt, RCGinv, Psi, Gamma,
                           A,       B,             C,           D,      SIGMA);
                           
        compute_higher_order_stress(Gamma_voigt, C, M);

        return;
    }

    void compute_A_voigt(const double &lambda,const double &mu, SpMat &A) {
        /*!=========================
           |    compute_A_voigt    |
           =========================
           
           Compute the A stiffness matrix in voigt notation.
        */
        
        std::vector<T> tripletList;
        tripletList.reserve(21);
   
        tripletList.push_back(T(0,0,lambda + 2*mu));
        tripletList.push_back(T(0,1,lambda));
        tripletList.push_back(T(0,2,lambda));
        tripletList.push_back(T(1,0,lambda));
        tripletList.push_back(T(1,1,lambda + 2*mu));
        tripletList.push_back(T(1,2,lambda));
        tripletList.push_back(T(2,0,lambda));
        tripletList.push_back(T(2,1,lambda));
        tripletList.push_back(T(2,2,lambda + 2*mu));
        tripletList.push_back(T(3,3,mu));
        tripletList.push_back(T(3,6,mu));
        tripletList.push_back(T(4,4,mu));
        tripletList.push_back(T(4,7,mu));
        tripletList.push_back(T(5,5,mu));
        tripletList.push_back(T(5,8,mu));
        tripletList.push_back(T(6,3,mu));
        tripletList.push_back(T(6,6,mu));
        tripletList.push_back(T(7,4,mu));
        tripletList.push_back(T(7,7,mu));
        tripletList.push_back(T(8,5,mu));
        tripletList.push_back(T(8,8,mu));
    
        A.setFromTriplets(tripletList.begin(), tripletList.end());
        return;
    }
    
    void compute_B_voigt(const double &eta,   const double &kappa, const double &nu,
                         const double &sigma, const double &tau,   SpMat &B) {
        /*!=========================
           |    compute_B_voigt    |
           =========================
           
           Compute the B stiffness matrix in voigt notation.
        */
        
        std::vector<T> tripletList;
        tripletList.reserve(21);
        
        tripletList.push_back(T(0,0,eta + kappa + nu - 2*sigma - tau));
        tripletList.push_back(T(0,1,eta - tau));
        tripletList.push_back(T(0,2,eta - tau));
        tripletList.push_back(T(1,0,eta - tau));
        tripletList.push_back(T(1,1,eta + kappa + nu - 2*sigma - tau));
        tripletList.push_back(T(1,2,eta - tau));
        tripletList.push_back(T(2,0,eta - tau));
        tripletList.push_back(T(2,1,eta - tau));
        tripletList.push_back(T(2,2,eta + kappa + nu - 2*sigma - tau));
        tripletList.push_back(T(3,3,kappa - sigma));
        tripletList.push_back(T(3,6,nu - sigma));
        tripletList.push_back(T(4,4,kappa - sigma));
        tripletList.push_back(T(4,7,nu - sigma));
        tripletList.push_back(T(5,5,kappa - sigma));
        tripletList.push_back(T(5,8,nu - sigma));
        tripletList.push_back(T(6,3,nu - sigma));
        tripletList.push_back(T(6,6,kappa - sigma));
        tripletList.push_back(T(7,4,nu - sigma));
        tripletList.push_back(T(7,7,kappa - sigma));
        tripletList.push_back(T(8,5,nu - sigma));
        tripletList.push_back(T(8,8,kappa - sigma));
        
        B.setFromTriplets(tripletList.begin(), tripletList.end());
        return;
    }
    
    void compute_C_voigt(const double &tau1,  const double &tau2,  const double &tau3,
                         const double &tau4,  const double &tau5,  const double &tau6,
                         const double &tau7,  const double &tau8,  const double &tau9,
                         const double &tau10, const double &tau11, SpMat &C) {
        /*!=========================
           |    compute_C_voigt    |
           =========================
           
        Compute the C stiffness tensor in voigt 
        format.
        
        */
        std::vector<T> tripletList;
        tripletList.reserve(183);

        tripletList.push_back(T(0,0,2*tau1 + tau10 + tau11 + 2*tau2 + tau3 + tau4 + 2*tau5 + tau6 + tau7 + 2*tau8 + tau9));
        tripletList.push_back(T(0,1,tau1 + tau4 + tau5));
        tripletList.push_back(T(0,2,tau1 + tau4 + tau5));
        tripletList.push_back(T(0,14,tau2 + tau5 + tau6));
        tripletList.push_back(T(0,17,tau1 + tau2 + tau3));
        tripletList.push_back(T(0,22,tau2 + tau5 + tau6));
        tripletList.push_back(T(0,25,tau1 + tau2 + tau3));
        tripletList.push_back(T(1,0,tau1 + tau4 + tau5));
        tripletList.push_back(T(1,1,tau4 + tau7 + tau9));
        tripletList.push_back(T(1,2,tau4));
        tripletList.push_back(T(1,14,tau10 + tau5 + tau8));
        tripletList.push_back(T(1,17,tau1 + tau11 + tau8));
        tripletList.push_back(T(1,22,tau5));
        tripletList.push_back(T(1,25,tau1));
        tripletList.push_back(T(2,0,tau1 + tau4 + tau5));
        tripletList.push_back(T(2,1,tau4));
        tripletList.push_back(T(2,2,tau4 + tau7 + tau9));
        tripletList.push_back(T(2,14,tau5));
        tripletList.push_back(T(2,17,tau1));
        tripletList.push_back(T(2,22,tau10 + tau5 + tau8));
        tripletList.push_back(T(2,25,tau1 + tau11 + tau8));
        tripletList.push_back(T(3,3,tau7));
        tripletList.push_back(T(3,6,tau9));
        tripletList.push_back(T(3,13,tau10));
        tripletList.push_back(T(3,16,tau8));
        tripletList.push_back(T(3,23,tau8));
        tripletList.push_back(T(3,26,tau11));
        tripletList.push_back(T(4,4,tau10 + tau3 + tau7));
        tripletList.push_back(T(4,7,tau2 + tau8 + tau9));
        tripletList.push_back(T(4,12,tau3));
        tripletList.push_back(T(4,15,tau2));
        tripletList.push_back(T(4,18,tau1 + tau11 + tau8));
        tripletList.push_back(T(4,19,tau1));
        tripletList.push_back(T(4,20,tau1 + tau2 + tau3));
        tripletList.push_back(T(5,5,tau10 + tau3 + tau7));
        tripletList.push_back(T(5,8,tau2 + tau8 + tau9));
        tripletList.push_back(T(5,9,tau1 + tau11 + tau8));
        tripletList.push_back(T(5,10,tau1 + tau2 + tau3));
        tripletList.push_back(T(5,11,tau1));
        tripletList.push_back(T(5,21,tau2));
        tripletList.push_back(T(5,24,tau3));
        tripletList.push_back(T(6,3,tau9));
        tripletList.push_back(T(6,6,tau7));
        tripletList.push_back(T(6,13,tau8));
        tripletList.push_back(T(6,16,tau11));
        tripletList.push_back(T(6,23,tau10));
        tripletList.push_back(T(6,26,tau8));
        tripletList.push_back(T(7,4,tau2 + tau8 + tau9));
        tripletList.push_back(T(7,7,tau11 + tau6 + tau7));
        tripletList.push_back(T(7,12,tau2));
        tripletList.push_back(T(7,15,tau6));
        tripletList.push_back(T(7,18,tau10 + tau5 + tau8));
        tripletList.push_back(T(7,19,tau5));
        tripletList.push_back(T(7,20,tau2 + tau5 + tau6));
        tripletList.push_back(T(8,5,tau2 + tau8 + tau9));
        tripletList.push_back(T(8,8,tau11 + tau6 + tau7));
        tripletList.push_back(T(8,9,tau10 + tau5 + tau8));
        tripletList.push_back(T(8,10,tau2 + tau5 + tau6));
        tripletList.push_back(T(8,11,tau5));
        tripletList.push_back(T(8,21,tau6));
        tripletList.push_back(T(8,24,tau2));
        tripletList.push_back(T(9,5,tau1 + tau11 + tau8));
        tripletList.push_back(T(9,8,tau10 + tau5 + tau8));
        tripletList.push_back(T(9,9,tau4 + tau7 + tau9));
        tripletList.push_back(T(9,10,tau1 + tau4 + tau5));
        tripletList.push_back(T(9,11,tau4));
        tripletList.push_back(T(9,21,tau5));
        tripletList.push_back(T(9,24,tau1));
        tripletList.push_back(T(10,5,tau1 + tau2 + tau3));
        tripletList.push_back(T(10,8,tau2 + tau5 + tau6));
        tripletList.push_back(T(10,9,tau1 + tau4 + tau5));
        tripletList.push_back(T(10,10,2*tau1 + tau10 + tau11 + 2*tau2 + tau3 + tau4 + 2*tau5 + tau6 + tau7 + 2*tau8 + tau9));
        tripletList.push_back(T(10,11,tau1 + tau4 + tau5));
        tripletList.push_back(T(10,21,tau2 + tau5 + tau6));
        tripletList.push_back(T(10,24,tau1 + tau2 + tau3));
        tripletList.push_back(T(11,5,tau1));
        tripletList.push_back(T(11,8,tau5));
        tripletList.push_back(T(11,9,tau4));
        tripletList.push_back(T(11,10,tau1 + tau4 + tau5));
        tripletList.push_back(T(11,11,tau4 + tau7 + tau9));
        tripletList.push_back(T(11,21,tau10 + tau5 + tau8));
        tripletList.push_back(T(11,24,tau1 + tau11 + tau8));
        tripletList.push_back(T(12,4,tau3));
        tripletList.push_back(T(12,7,tau2));
        tripletList.push_back(T(12,12,tau10 + tau3 + tau7));
        tripletList.push_back(T(12,15,tau2 + tau8 + tau9));
        tripletList.push_back(T(12,18,tau1));
        tripletList.push_back(T(12,19,tau1 + tau11 + tau8));
        tripletList.push_back(T(12,20,tau1 + tau2 + tau3));
        tripletList.push_back(T(13,3,tau10));
        tripletList.push_back(T(13,6,tau8));
        tripletList.push_back(T(13,13,tau7));
        tripletList.push_back(T(13,16,tau9));
        tripletList.push_back(T(13,23,tau11));
        tripletList.push_back(T(13,26,tau8));
        tripletList.push_back(T(14,0,tau2 + tau5 + tau6));
        tripletList.push_back(T(14,1,tau10 + tau5 + tau8));
        tripletList.push_back(T(14,2,tau5));
        tripletList.push_back(T(14,14,tau11 + tau6 + tau7));
        tripletList.push_back(T(14,17,tau2 + tau8 + tau9));
        tripletList.push_back(T(14,22,tau6));
        tripletList.push_back(T(14,25,tau2));
        tripletList.push_back(T(15,4,tau2));
        tripletList.push_back(T(15,7,tau6));
        tripletList.push_back(T(15,12,tau2 + tau8 + tau9));
        tripletList.push_back(T(15,15,tau11 + tau6 + tau7));
        tripletList.push_back(T(15,18,tau5));
        tripletList.push_back(T(15,19,tau10 + tau5 + tau8));
        tripletList.push_back(T(15,20,tau2 + tau5 + tau6));
        tripletList.push_back(T(16,3,tau8));
        tripletList.push_back(T(16,6,tau11));
        tripletList.push_back(T(16,13,tau9));
        tripletList.push_back(T(16,16,tau7));
        tripletList.push_back(T(16,23,tau8));
        tripletList.push_back(T(16,26,tau10));
        tripletList.push_back(T(17,0,tau1 + tau2 + tau3));
        tripletList.push_back(T(17,1,tau1 + tau11 + tau8));
        tripletList.push_back(T(17,2,tau1));
        tripletList.push_back(T(17,14,tau2 + tau8 + tau9));
        tripletList.push_back(T(17,17,tau10 + tau3 + tau7));
        tripletList.push_back(T(17,22,tau2));
        tripletList.push_back(T(17,25,tau3));
        tripletList.push_back(T(18,4,tau1 + tau11 + tau8));
        tripletList.push_back(T(18,7,tau10 + tau5 + tau8));
        tripletList.push_back(T(18,12,tau1));
        tripletList.push_back(T(18,15,tau5));
        tripletList.push_back(T(18,18,tau4 + tau7 + tau9));
        tripletList.push_back(T(18,19,tau4));
        tripletList.push_back(T(18,20,tau1 + tau4 + tau5));
        tripletList.push_back(T(19,4,tau1));
        tripletList.push_back(T(19,7,tau5));
        tripletList.push_back(T(19,12,tau1 + tau11 + tau8));
        tripletList.push_back(T(19,15,tau10 + tau5 + tau8));
        tripletList.push_back(T(19,18,tau4));
        tripletList.push_back(T(19,19,tau4 + tau7 + tau9));
        tripletList.push_back(T(19,20,tau1 + tau4 + tau5));
        tripletList.push_back(T(20,4,tau1 + tau2 + tau3));
        tripletList.push_back(T(20,7,tau2 + tau5 + tau6));
        tripletList.push_back(T(20,12,tau1 + tau2 + tau3));
        tripletList.push_back(T(20,15,tau2 + tau5 + tau6));
        tripletList.push_back(T(20,18,tau1 + tau4 + tau5));
        tripletList.push_back(T(20,19,tau1 + tau4 + tau5));
        tripletList.push_back(T(20,20,2*tau1 + tau10 + tau11 + 2*tau2 + tau3 + tau4 + 2*tau5 + tau6 + tau7 + 2*tau8 + tau9));
        tripletList.push_back(T(21,5,tau2));
        tripletList.push_back(T(21,8,tau6));
        tripletList.push_back(T(21,9,tau5));
        tripletList.push_back(T(21,10,tau2 + tau5 + tau6));
        tripletList.push_back(T(21,11,tau10 + tau5 + tau8));
        tripletList.push_back(T(21,21,tau11 + tau6 + tau7));
        tripletList.push_back(T(21,24,tau2 + tau8 + tau9));
        tripletList.push_back(T(22,0,tau2 + tau5 + tau6));
        tripletList.push_back(T(22,1,tau5));
        tripletList.push_back(T(22,2,tau10 + tau5 + tau8));
        tripletList.push_back(T(22,14,tau6));
        tripletList.push_back(T(22,17,tau2));
        tripletList.push_back(T(22,22,tau11 + tau6 + tau7));
        tripletList.push_back(T(22,25,tau2 + tau8 + tau9));
        tripletList.push_back(T(23,3,tau8));
        tripletList.push_back(T(23,6,tau10));
        tripletList.push_back(T(23,13,tau11));
        tripletList.push_back(T(23,16,tau8));
        tripletList.push_back(T(23,23,tau7));
        tripletList.push_back(T(23,26,tau9));
        tripletList.push_back(T(24,5,tau3));
        tripletList.push_back(T(24,8,tau2));
        tripletList.push_back(T(24,9,tau1));
        tripletList.push_back(T(24,10,tau1 + tau2 + tau3));
        tripletList.push_back(T(24,11,tau1 + tau11 + tau8));
        tripletList.push_back(T(24,21,tau2 + tau8 + tau9));
        tripletList.push_back(T(24,24,tau10 + tau3 + tau7));
        tripletList.push_back(T(25,0,tau1 + tau2 + tau3));
        tripletList.push_back(T(25,1,tau1));
        tripletList.push_back(T(25,2,tau1 + tau11 + tau8));
        tripletList.push_back(T(25,14,tau2));
        tripletList.push_back(T(25,17,tau3));
        tripletList.push_back(T(25,22,tau2 + tau8 + tau9));
        tripletList.push_back(T(25,25,tau10 + tau3 + tau7));
        tripletList.push_back(T(26,3,tau11));
        tripletList.push_back(T(26,6,tau8));
        tripletList.push_back(T(26,13,tau8));
        tripletList.push_back(T(26,16,tau10));
        tripletList.push_back(T(26,23,tau9));
        tripletList.push_back(T(26,26,tau7));
        
        C.setFromTriplets(tripletList.begin(), tripletList.end());
        return;
    }
    
    void compute_D_voigt(const double &sigma, const double &tau, SpMat &D){
        /*!=========================
           |    compute_D_voigt    |
           =========================
           
        Compute the D stiffness tensor in voigt 
        format.
        
        */
        std::vector<T> tripletList;
        tripletList.reserve(21);

        tripletList.push_back(T(0,0,2*sigma + tau));
        tripletList.push_back(T(0,1,tau));
        tripletList.push_back(T(0,2,tau));
        tripletList.push_back(T(1,0,tau));
        tripletList.push_back(T(1,1,2*sigma + tau));
        tripletList.push_back(T(1,2,tau));
        tripletList.push_back(T(2,0,tau));
        tripletList.push_back(T(2,1,tau));
        tripletList.push_back(T(2,2,2*sigma + tau));
        tripletList.push_back(T(3,3,sigma));
        tripletList.push_back(T(3,6,sigma));
        tripletList.push_back(T(4,4,sigma));
        tripletList.push_back(T(4,7,sigma));
        tripletList.push_back(T(5,5,sigma));
        tripletList.push_back(T(5,8,sigma));
        tripletList.push_back(T(6,3,sigma));
        tripletList.push_back(T(6,6,sigma));
        tripletList.push_back(T(7,4,sigma));
        tripletList.push_back(T(7,7,sigma));
        tripletList.push_back(T(8,5,sigma));
        tripletList.push_back(T(8,8,sigma));       
        
        D.setFromTriplets(tripletList.begin(), tripletList.end());
        return;
    }
    
    void compute_PK2_stress(const Vector_9 &E_voigt,    const Vector_9 &E_micro_voigt, const Vector_27 &Gamma_voigt,
                            const Matrix_3x3 &RCGinv,   const Matrix_3x3 &Psi,         const Matrix_3x9 &Gamma,
                            const SpMat &A, const SpMat &B,    const SpMat &C,
                            const SpMat &D, Vector_9 &PK2){
        /*!============================
        |    compute_PK2_stress    |
        ============================
           
        Compute the second piola kirchoff stress.
           
        */
        
        PK2 = A*E_voigt;        //Compute the first terms
        PK2 += D*E_micro_voigt;
        
        //Compute the middle terms
        Matrix_3x3 Temp1;
        deformation_measures::undo_voigt_3x3_tensor(B*E_micro_voigt+D*E_voigt,Temp1);
        Vector_9 term3_4_voigt;
        deformation_measures::voigt_3x3_tensor(Temp1*(RCGinv*Psi).transpose(),term3_4_voigt);
        
        PK2 += term3_4_voigt;
        
        //Compute the end terms
        Matrix_3x9 Temp2;
        deformation_measures::undo_voigt_3x9_tensor(C*Gamma_voigt,Temp2);
        Vector_9 term5_voigt;
        deformation_measures::voigt_3x3_tensor(Temp2*(RCGinv*Gamma).transpose(),term5_voigt);
        PK2 += term5_voigt;
        return;
    }
    
    void compute_symmetric_stress(const Vector_9 &E_voigt,    const Vector_9 &E_micro_voigt, const Vector_27 &Gamma_voigt,
                            const Matrix_3x3 &RCGinv,     const Matrix_3x3 &Psi,         const Matrix_3x9 &Gamma,
                            const SpMat &A, const SpMat &B,    const SpMat &C,
                            const SpMat &D, Vector_9 &SIGMA){
        /*!=====================================
           |    compute_symmetric_stress    |
           ==================================
           
           Compute the symmetric stress in the reference configuration.
           
        */
        
        SIGMA = A*E_voigt;        //Compute the first terms
        SIGMA += D*E_micro_voigt;
        
        //Compute the middle terms
        Matrix_3x3 Temp1;
        deformation_measures::undo_voigt_3x3_tensor(B*E_micro_voigt+D*E_voigt,Temp1);
        Matrix_3x3 symmetric_part = Temp1*(RCGinv*Psi).transpose();
        
        //Compute the end terms
        Matrix_3x9 Temp2;
        deformation_measures::undo_voigt_3x9_tensor(C*Gamma_voigt,Temp2);
        symmetric_part += Temp2*(RCGinv*Gamma).transpose();
        Vector_9 vector_symm_part;
        deformation_measures::voigt_3x3_tensor((symmetric_part + symmetric_part.transpose()),vector_symm_part);
        
        SIGMA += vector_symm_part;
        
        return;
    }
    
    void compute_higher_order_stress(const Vector_27 &Gamma_voigt, const SpMat &C, Vector_27 &M){
        /*!=====================================
        |    compute_higher_order_stress    |
        =====================================
          
        Compute the higher order stress in the reference configuration.
          
        */
        
        M = C*Gamma_voigt; //Compute the stress (requires positive permutation)
        deformation_measures::perform_right_positive_cyclic_permutation(M); //Perform the permutation
    }

    void map_stresses_to_current_configuration(const Matrix_3x3 &F, const Matrix_3x3 &chi,
                                               const Vector_9 &PK2, const Vector_9 &SIGMA, const Vector_27 &M,
                                               Vector_9 &cauchy, Vector_9 &s, Vector_27 &m){
        /*!===============================================
        |    map_stresses_to_current_configuration    |
        ===============================================

        Map the stresses to the current configuration.

        */

        //Compute the jacobian of deformation
        double Jac = F.determinant();

        //Map the PK2 stress to the cauchy stress
        Matrix_3x3 PK2_mat;
        deformation_measures::undo_voigt_3x3_tensor(PK2,PK2_mat);
        deformation_measures::voigt_3x3_tensor(F*PK2_mat*F.transpose()/Jac,cauchy);

        //Map the symmetric stress to the current configuration
        Matrix_3x3 SIGMA_mat;
        deformation_measures::undo_voigt_3x3_tensor(SIGMA,SIGMA_mat);
        deformation_measures::voigt_3x3_tensor(F*SIGMA_mat*F.transpose()/Jac,s);

        //Map the higher order stress to the current configuration
        Matrix_3x9 M_mat;
        deformation_measures::undo_voigt_3x9_tensor(M,M_mat);
        deformation_measures::voigt_3x9_tensor(F*M_mat,m);                 //Map the first index
        deformation_measures::perform_left_positive_cyclic_permutation(m); //Cycle the indices
        deformation_measures::undo_voigt_3x9_tensor(m,M_mat);
        deformation_measures::voigt_3x9_tensor(F*M_mat,m);                 //Map the second index
        deformation_measures::perform_left_positive_cyclic_permutation(m); //Cycle the indices
        deformation_measures::undo_voigt_3x9_tensor(m,M_mat);
        deformation_measures::voigt_3x9_tensor(chi*M_mat,m);               //Map the third index
        deformation_measures::perform_left_positive_cyclic_permutation(m); //Cycle the indices
        m = m/Jac;

        return;
    }

    void map_fot_stress_jacobian(const Matrix_3x3 &F, const Matrix_9x9 &dSdfot, const double &J, Matrix_9x9 &term){
        /*!=================================
        |    map_fot_stress_jacobian    |
        =================================

        Map a jacobian of a second order stress 
        tensor w.r.t. a second order deformation measure.

        */
        
        //Extract the deformation gradient
        double F11 = F(0,0);
        double F12 = F(0,1);
        double F13 = F(0,2);
        double F21 = F(1,0);
        double F22 = F(1,1);
        double F23 = F(1,2);
        double F31 = F(2,0);
        double F32 = F(2,1);
        double F33 = F(2,2);

        //Extract the jacobian
        double dSdfot1111 = dSdfot(0,0);
        double dSdfot1112 = dSdfot(0,5);
        double dSdfot1113 = dSdfot(0,4);
        double dSdfot1121 = dSdfot(0,8);
        double dSdfot1122 = dSdfot(0,1);
        double dSdfot1123 = dSdfot(0,3);
        double dSdfot1131 = dSdfot(0,7);
        double dSdfot1132 = dSdfot(0,6);
        double dSdfot1133 = dSdfot(0,2);
        double dSdfot1211 = dSdfot(5,0);
        double dSdfot1212 = dSdfot(5,5);
        double dSdfot1213 = dSdfot(5,4);
        double dSdfot1221 = dSdfot(5,8);
        double dSdfot1222 = dSdfot(5,1);
        double dSdfot1223 = dSdfot(5,3);
        double dSdfot1231 = dSdfot(5,7);
        double dSdfot1232 = dSdfot(5,6);
        double dSdfot1233 = dSdfot(5,2);
        double dSdfot1311 = dSdfot(4,0);
        double dSdfot1312 = dSdfot(4,5);
        double dSdfot1313 = dSdfot(4,4);
        double dSdfot1321 = dSdfot(4,8);
        double dSdfot1322 = dSdfot(4,1);
        double dSdfot1323 = dSdfot(4,3);
        double dSdfot1331 = dSdfot(4,7);
        double dSdfot1332 = dSdfot(4,6);
        double dSdfot1333 = dSdfot(4,2);
        double dSdfot2111 = dSdfot(8,0);
        double dSdfot2112 = dSdfot(8,5);
        double dSdfot2113 = dSdfot(8,4);
        double dSdfot2121 = dSdfot(8,8);
        double dSdfot2122 = dSdfot(8,1);
        double dSdfot2123 = dSdfot(8,3);
        double dSdfot2131 = dSdfot(8,7);
        double dSdfot2132 = dSdfot(8,6);
        double dSdfot2133 = dSdfot(8,2);
        double dSdfot2211 = dSdfot(1,0);
        double dSdfot2212 = dSdfot(1,5);
        double dSdfot2213 = dSdfot(1,4);
        double dSdfot2221 = dSdfot(1,8);
        double dSdfot2222 = dSdfot(1,1);
        double dSdfot2223 = dSdfot(1,3);
        double dSdfot2231 = dSdfot(1,7);
        double dSdfot2232 = dSdfot(1,6);
        double dSdfot2233 = dSdfot(1,2);
        double dSdfot2311 = dSdfot(3,0);
        double dSdfot2312 = dSdfot(3,5);
        double dSdfot2313 = dSdfot(3,4);
        double dSdfot2321 = dSdfot(3,8);
        double dSdfot2322 = dSdfot(3,1);
        double dSdfot2323 = dSdfot(3,3);
        double dSdfot2331 = dSdfot(3,7);
        double dSdfot2332 = dSdfot(3,6);
        double dSdfot2333 = dSdfot(3,2);
        double dSdfot3111 = dSdfot(7,0);
        double dSdfot3112 = dSdfot(7,5);
        double dSdfot3113 = dSdfot(7,4);
        double dSdfot3121 = dSdfot(7,8);
        double dSdfot3122 = dSdfot(7,1);
        double dSdfot3123 = dSdfot(7,3);
        double dSdfot3131 = dSdfot(7,7);
        double dSdfot3132 = dSdfot(7,6);
        double dSdfot3133 = dSdfot(7,2);
        double dSdfot3211 = dSdfot(6,0);
        double dSdfot3212 = dSdfot(6,5);
        double dSdfot3213 = dSdfot(6,4);
        double dSdfot3221 = dSdfot(6,8);
        double dSdfot3222 = dSdfot(6,1);
        double dSdfot3223 = dSdfot(6,3);
        double dSdfot3231 = dSdfot(6,7);
        double dSdfot3232 = dSdfot(6,6);
        double dSdfot3233 = dSdfot(6,2);
        double dSdfot3311 = dSdfot(2,0);
        double dSdfot3312 = dSdfot(2,5);
        double dSdfot3313 = dSdfot(2,4);
        double dSdfot3321 = dSdfot(2,8);
        double dSdfot3322 = dSdfot(2,1);
        double dSdfot3323 = dSdfot(2,3);
        double dSdfot3331 = dSdfot(2,7);
        double dSdfot3332 = dSdfot(2,6);
        double dSdfot3333 = dSdfot(2,2);

        //Compute term
        term(0,0) = pow(F11, 2)*dSdfot1111 + F11*F12*dSdfot1211 + F11*F12*dSdfot2111 + F11*F13*dSdfot1311 + F11*F13*dSdfot3111 + pow(F12, 2)*dSdfot2211 + F12*F13*dSdfot2311 + F12*F13*dSdfot3211 + pow(F13, 2)*dSdfot3311;
        term(0,1) = pow(F11, 2)*dSdfot1122 + F11*F12*dSdfot1222 + F11*F12*dSdfot2122 + F11*F13*dSdfot1322 + F11*F13*dSdfot3122 + pow(F12, 2)*dSdfot2222 + F12*F13*dSdfot2322 + F12*F13*dSdfot3222 + pow(F13, 2)*dSdfot3322;
        term(0,2) = pow(F11, 2)*dSdfot1133 + F11*F12*dSdfot1233 + F11*F12*dSdfot2133 + F11*F13*dSdfot1333 + F11*F13*dSdfot3133 + pow(F12, 2)*dSdfot2233 + F12*F13*dSdfot2333 + F12*F13*dSdfot3233 + pow(F13, 2)*dSdfot3333;
        term(0,3) = pow(F11, 2)*dSdfot1123 + F11*F12*dSdfot1223 + F11*F12*dSdfot2123 + F11*F13*dSdfot1323 + F11*F13*dSdfot3123 + pow(F12, 2)*dSdfot2223 + F12*F13*dSdfot2323 + F12*F13*dSdfot3223 + pow(F13, 2)*dSdfot3323;
        term(0,4) = pow(F11, 2)*dSdfot1113 + F11*F12*dSdfot1213 + F11*F12*dSdfot2113 + F11*F13*dSdfot1313 + F11*F13*dSdfot3113 + pow(F12, 2)*dSdfot2213 + F12*F13*dSdfot2313 + F12*F13*dSdfot3213 + pow(F13, 2)*dSdfot3313;
        term(0,5) = pow(F11, 2)*dSdfot1112 + F11*F12*dSdfot1212 + F11*F12*dSdfot2112 + F11*F13*dSdfot1312 + F11*F13*dSdfot3112 + pow(F12, 2)*dSdfot2212 + F12*F13*dSdfot2312 + F12*F13*dSdfot3212 + pow(F13, 2)*dSdfot3312;
        term(0,6) = pow(F11, 2)*dSdfot1132 + F11*F12*dSdfot1232 + F11*F12*dSdfot2132 + F11*F13*dSdfot1332 + F11*F13*dSdfot3132 + pow(F12, 2)*dSdfot2232 + F12*F13*dSdfot2332 + F12*F13*dSdfot3232 + pow(F13, 2)*dSdfot3332;
        term(0,7) = pow(F11, 2)*dSdfot1131 + F11*F12*dSdfot1231 + F11*F12*dSdfot2131 + F11*F13*dSdfot1331 + F11*F13*dSdfot3131 + pow(F12, 2)*dSdfot2231 + F12*F13*dSdfot2331 + F12*F13*dSdfot3231 + pow(F13, 2)*dSdfot3331;
        term(0,8) = pow(F11, 2)*dSdfot1121 + F11*F12*dSdfot1221 + F11*F12*dSdfot2121 + F11*F13*dSdfot1321 + F11*F13*dSdfot3121 + pow(F12, 2)*dSdfot2221 + F12*F13*dSdfot2321 + F12*F13*dSdfot3221 + pow(F13, 2)*dSdfot3321;
        term(1,0) = pow(F21, 2)*dSdfot1111 + F21*F22*dSdfot1211 + F21*F22*dSdfot2111 + F21*F23*dSdfot1311 + F21*F23*dSdfot3111 + pow(F22, 2)*dSdfot2211 + F22*F23*dSdfot2311 + F22*F23*dSdfot3211 + pow(F23, 2)*dSdfot3311;
        term(1,1) = pow(F21, 2)*dSdfot1122 + F21*F22*dSdfot1222 + F21*F22*dSdfot2122 + F21*F23*dSdfot1322 + F21*F23*dSdfot3122 + pow(F22, 2)*dSdfot2222 + F22*F23*dSdfot2322 + F22*F23*dSdfot3222 + pow(F23, 2)*dSdfot3322;
        term(1,2) = pow(F21, 2)*dSdfot1133 + F21*F22*dSdfot1233 + F21*F22*dSdfot2133 + F21*F23*dSdfot1333 + F21*F23*dSdfot3133 + pow(F22, 2)*dSdfot2233 + F22*F23*dSdfot2333 + F22*F23*dSdfot3233 + pow(F23, 2)*dSdfot3333;
        term(1,3) = pow(F21, 2)*dSdfot1123 + F21*F22*dSdfot1223 + F21*F22*dSdfot2123 + F21*F23*dSdfot1323 + F21*F23*dSdfot3123 + pow(F22, 2)*dSdfot2223 + F22*F23*dSdfot2323 + F22*F23*dSdfot3223 + pow(F23, 2)*dSdfot3323;
        term(1,4) = pow(F21, 2)*dSdfot1113 + F21*F22*dSdfot1213 + F21*F22*dSdfot2113 + F21*F23*dSdfot1313 + F21*F23*dSdfot3113 + pow(F22, 2)*dSdfot2213 + F22*F23*dSdfot2313 + F22*F23*dSdfot3213 + pow(F23, 2)*dSdfot3313;
        term(1,5) = pow(F21, 2)*dSdfot1112 + F21*F22*dSdfot1212 + F21*F22*dSdfot2112 + F21*F23*dSdfot1312 + F21*F23*dSdfot3112 + pow(F22, 2)*dSdfot2212 + F22*F23*dSdfot2312 + F22*F23*dSdfot3212 + pow(F23, 2)*dSdfot3312;
        term(1,6) = pow(F21, 2)*dSdfot1132 + F21*F22*dSdfot1232 + F21*F22*dSdfot2132 + F21*F23*dSdfot1332 + F21*F23*dSdfot3132 + pow(F22, 2)*dSdfot2232 + F22*F23*dSdfot2332 + F22*F23*dSdfot3232 + pow(F23, 2)*dSdfot3332;
        term(1,7) = pow(F21, 2)*dSdfot1131 + F21*F22*dSdfot1231 + F21*F22*dSdfot2131 + F21*F23*dSdfot1331 + F21*F23*dSdfot3131 + pow(F22, 2)*dSdfot2231 + F22*F23*dSdfot2331 + F22*F23*dSdfot3231 + pow(F23, 2)*dSdfot3331;
        term(1,8) = pow(F21, 2)*dSdfot1121 + F21*F22*dSdfot1221 + F21*F22*dSdfot2121 + F21*F23*dSdfot1321 + F21*F23*dSdfot3121 + pow(F22, 2)*dSdfot2221 + F22*F23*dSdfot2321 + F22*F23*dSdfot3221 + pow(F23, 2)*dSdfot3321;
        term(2,0) = pow(F31, 2)*dSdfot1111 + F31*F32*dSdfot1211 + F31*F32*dSdfot2111 + F31*F33*dSdfot1311 + F31*F33*dSdfot3111 + pow(F32, 2)*dSdfot2211 + F32*F33*dSdfot2311 + F32*F33*dSdfot3211 + pow(F33, 2)*dSdfot3311;
        term(2,1) = pow(F31, 2)*dSdfot1122 + F31*F32*dSdfot1222 + F31*F32*dSdfot2122 + F31*F33*dSdfot1322 + F31*F33*dSdfot3122 + pow(F32, 2)*dSdfot2222 + F32*F33*dSdfot2322 + F32*F33*dSdfot3222 + pow(F33, 2)*dSdfot3322;
        term(2,2) = pow(F31, 2)*dSdfot1133 + F31*F32*dSdfot1233 + F31*F32*dSdfot2133 + F31*F33*dSdfot1333 + F31*F33*dSdfot3133 + pow(F32, 2)*dSdfot2233 + F32*F33*dSdfot2333 + F32*F33*dSdfot3233 + pow(F33, 2)*dSdfot3333;
        term(2,3) = pow(F31, 2)*dSdfot1123 + F31*F32*dSdfot1223 + F31*F32*dSdfot2123 + F31*F33*dSdfot1323 + F31*F33*dSdfot3123 + pow(F32, 2)*dSdfot2223 + F32*F33*dSdfot2323 + F32*F33*dSdfot3223 + pow(F33, 2)*dSdfot3323;
        term(2,4) = pow(F31, 2)*dSdfot1113 + F31*F32*dSdfot1213 + F31*F32*dSdfot2113 + F31*F33*dSdfot1313 + F31*F33*dSdfot3113 + pow(F32, 2)*dSdfot2213 + F32*F33*dSdfot2313 + F32*F33*dSdfot3213 + pow(F33, 2)*dSdfot3313;
        term(2,5) = pow(F31, 2)*dSdfot1112 + F31*F32*dSdfot1212 + F31*F32*dSdfot2112 + F31*F33*dSdfot1312 + F31*F33*dSdfot3112 + pow(F32, 2)*dSdfot2212 + F32*F33*dSdfot2312 + F32*F33*dSdfot3212 + pow(F33, 2)*dSdfot3312;
        term(2,6) = pow(F31, 2)*dSdfot1132 + F31*F32*dSdfot1232 + F31*F32*dSdfot2132 + F31*F33*dSdfot1332 + F31*F33*dSdfot3132 + pow(F32, 2)*dSdfot2232 + F32*F33*dSdfot2332 + F32*F33*dSdfot3232 + pow(F33, 2)*dSdfot3332;
        term(2,7) = pow(F31, 2)*dSdfot1131 + F31*F32*dSdfot1231 + F31*F32*dSdfot2131 + F31*F33*dSdfot1331 + F31*F33*dSdfot3131 + pow(F32, 2)*dSdfot2231 + F32*F33*dSdfot2331 + F32*F33*dSdfot3231 + pow(F33, 2)*dSdfot3331;
        term(2,8) = pow(F31, 2)*dSdfot1121 + F31*F32*dSdfot1221 + F31*F32*dSdfot2121 + F31*F33*dSdfot1321 + F31*F33*dSdfot3121 + pow(F32, 2)*dSdfot2221 + F32*F33*dSdfot2321 + F32*F33*dSdfot3221 + pow(F33, 2)*dSdfot3321;
        term(3,0) = F21*F31*dSdfot1111 + F21*F32*dSdfot1211 + F21*F33*dSdfot1311 + F22*F31*dSdfot2111 + F22*F32*dSdfot2211 + F22*F33*dSdfot2311 + F23*F31*dSdfot3111 + F23*F32*dSdfot3211 + F23*F33*dSdfot3311;
        term(3,1) = F21*F31*dSdfot1122 + F21*F32*dSdfot1222 + F21*F33*dSdfot1322 + F22*F31*dSdfot2122 + F22*F32*dSdfot2222 + F22*F33*dSdfot2322 + F23*F31*dSdfot3122 + F23*F32*dSdfot3222 + F23*F33*dSdfot3322;
        term(3,2) = F21*F31*dSdfot1133 + F21*F32*dSdfot1233 + F21*F33*dSdfot1333 + F22*F31*dSdfot2133 + F22*F32*dSdfot2233 + F22*F33*dSdfot2333 + F23*F31*dSdfot3133 + F23*F32*dSdfot3233 + F23*F33*dSdfot3333;
        term(3,3) = F21*F31*dSdfot1123 + F21*F32*dSdfot1223 + F21*F33*dSdfot1323 + F22*F31*dSdfot2123 + F22*F32*dSdfot2223 + F22*F33*dSdfot2323 + F23*F31*dSdfot3123 + F23*F32*dSdfot3223 + F23*F33*dSdfot3323;
        term(3,4) = F21*F31*dSdfot1113 + F21*F32*dSdfot1213 + F21*F33*dSdfot1313 + F22*F31*dSdfot2113 + F22*F32*dSdfot2213 + F22*F33*dSdfot2313 + F23*F31*dSdfot3113 + F23*F32*dSdfot3213 + F23*F33*dSdfot3313;
        term(3,5) = F21*F31*dSdfot1112 + F21*F32*dSdfot1212 + F21*F33*dSdfot1312 + F22*F31*dSdfot2112 + F22*F32*dSdfot2212 + F22*F33*dSdfot2312 + F23*F31*dSdfot3112 + F23*F32*dSdfot3212 + F23*F33*dSdfot3312;
        term(3,6) = F21*F31*dSdfot1132 + F21*F32*dSdfot1232 + F21*F33*dSdfot1332 + F22*F31*dSdfot2132 + F22*F32*dSdfot2232 + F22*F33*dSdfot2332 + F23*F31*dSdfot3132 + F23*F32*dSdfot3232 + F23*F33*dSdfot3332;
        term(3,7) = F21*F31*dSdfot1131 + F21*F32*dSdfot1231 + F21*F33*dSdfot1331 + F22*F31*dSdfot2131 + F22*F32*dSdfot2231 + F22*F33*dSdfot2331 + F23*F31*dSdfot3131 + F23*F32*dSdfot3231 + F23*F33*dSdfot3331;
        term(3,8) = F21*F31*dSdfot1121 + F21*F32*dSdfot1221 + F21*F33*dSdfot1321 + F22*F31*dSdfot2121 + F22*F32*dSdfot2221 + F22*F33*dSdfot2321 + F23*F31*dSdfot3121 + F23*F32*dSdfot3221 + F23*F33*dSdfot3321;
        term(4,0) = F11*F31*dSdfot1111 + F11*F32*dSdfot1211 + F11*F33*dSdfot1311 + F12*F31*dSdfot2111 + F12*F32*dSdfot2211 + F12*F33*dSdfot2311 + F13*F31*dSdfot3111 + F13*F32*dSdfot3211 + F13*F33*dSdfot3311;
        term(4,1) = F11*F31*dSdfot1122 + F11*F32*dSdfot1222 + F11*F33*dSdfot1322 + F12*F31*dSdfot2122 + F12*F32*dSdfot2222 + F12*F33*dSdfot2322 + F13*F31*dSdfot3122 + F13*F32*dSdfot3222 + F13*F33*dSdfot3322;
        term(4,2) = F11*F31*dSdfot1133 + F11*F32*dSdfot1233 + F11*F33*dSdfot1333 + F12*F31*dSdfot2133 + F12*F32*dSdfot2233 + F12*F33*dSdfot2333 + F13*F31*dSdfot3133 + F13*F32*dSdfot3233 + F13*F33*dSdfot3333;
        term(4,3) = F11*F31*dSdfot1123 + F11*F32*dSdfot1223 + F11*F33*dSdfot1323 + F12*F31*dSdfot2123 + F12*F32*dSdfot2223 + F12*F33*dSdfot2323 + F13*F31*dSdfot3123 + F13*F32*dSdfot3223 + F13*F33*dSdfot3323;
        term(4,4) = F11*F31*dSdfot1113 + F11*F32*dSdfot1213 + F11*F33*dSdfot1313 + F12*F31*dSdfot2113 + F12*F32*dSdfot2213 + F12*F33*dSdfot2313 + F13*F31*dSdfot3113 + F13*F32*dSdfot3213 + F13*F33*dSdfot3313;
        term(4,5) = F11*F31*dSdfot1112 + F11*F32*dSdfot1212 + F11*F33*dSdfot1312 + F12*F31*dSdfot2112 + F12*F32*dSdfot2212 + F12*F33*dSdfot2312 + F13*F31*dSdfot3112 + F13*F32*dSdfot3212 + F13*F33*dSdfot3312;
        term(4,6) = F11*F31*dSdfot1132 + F11*F32*dSdfot1232 + F11*F33*dSdfot1332 + F12*F31*dSdfot2132 + F12*F32*dSdfot2232 + F12*F33*dSdfot2332 + F13*F31*dSdfot3132 + F13*F32*dSdfot3232 + F13*F33*dSdfot3332;
        term(4,7) = F11*F31*dSdfot1131 + F11*F32*dSdfot1231 + F11*F33*dSdfot1331 + F12*F31*dSdfot2131 + F12*F32*dSdfot2231 + F12*F33*dSdfot2331 + F13*F31*dSdfot3131 + F13*F32*dSdfot3231 + F13*F33*dSdfot3331;
        term(4,8) = F11*F31*dSdfot1121 + F11*F32*dSdfot1221 + F11*F33*dSdfot1321 + F12*F31*dSdfot2121 + F12*F32*dSdfot2221 + F12*F33*dSdfot2321 + F13*F31*dSdfot3121 + F13*F32*dSdfot3221 + F13*F33*dSdfot3321;
        term(5,0) = F11*F21*dSdfot1111 + F11*F22*dSdfot1211 + F11*F23*dSdfot1311 + F12*F21*dSdfot2111 + F12*F22*dSdfot2211 + F12*F23*dSdfot2311 + F13*F21*dSdfot3111 + F13*F22*dSdfot3211 + F13*F23*dSdfot3311;
        term(5,1) = F11*F21*dSdfot1122 + F11*F22*dSdfot1222 + F11*F23*dSdfot1322 + F12*F21*dSdfot2122 + F12*F22*dSdfot2222 + F12*F23*dSdfot2322 + F13*F21*dSdfot3122 + F13*F22*dSdfot3222 + F13*F23*dSdfot3322;
        term(5,2) = F11*F21*dSdfot1133 + F11*F22*dSdfot1233 + F11*F23*dSdfot1333 + F12*F21*dSdfot2133 + F12*F22*dSdfot2233 + F12*F23*dSdfot2333 + F13*F21*dSdfot3133 + F13*F22*dSdfot3233 + F13*F23*dSdfot3333;
        term(5,3) = F11*F21*dSdfot1123 + F11*F22*dSdfot1223 + F11*F23*dSdfot1323 + F12*F21*dSdfot2123 + F12*F22*dSdfot2223 + F12*F23*dSdfot2323 + F13*F21*dSdfot3123 + F13*F22*dSdfot3223 + F13*F23*dSdfot3323;
        term(5,4) = F11*F21*dSdfot1113 + F11*F22*dSdfot1213 + F11*F23*dSdfot1313 + F12*F21*dSdfot2113 + F12*F22*dSdfot2213 + F12*F23*dSdfot2313 + F13*F21*dSdfot3113 + F13*F22*dSdfot3213 + F13*F23*dSdfot3313;
        term(5,5) = F11*F21*dSdfot1112 + F11*F22*dSdfot1212 + F11*F23*dSdfot1312 + F12*F21*dSdfot2112 + F12*F22*dSdfot2212 + F12*F23*dSdfot2312 + F13*F21*dSdfot3112 + F13*F22*dSdfot3212 + F13*F23*dSdfot3312;
        term(5,6) = F11*F21*dSdfot1132 + F11*F22*dSdfot1232 + F11*F23*dSdfot1332 + F12*F21*dSdfot2132 + F12*F22*dSdfot2232 + F12*F23*dSdfot2332 + F13*F21*dSdfot3132 + F13*F22*dSdfot3232 + F13*F23*dSdfot3332;
        term(5,7) = F11*F21*dSdfot1131 + F11*F22*dSdfot1231 + F11*F23*dSdfot1331 + F12*F21*dSdfot2131 + F12*F22*dSdfot2231 + F12*F23*dSdfot2331 + F13*F21*dSdfot3131 + F13*F22*dSdfot3231 + F13*F23*dSdfot3331;
        term(5,8) = F11*F21*dSdfot1121 + F11*F22*dSdfot1221 + F11*F23*dSdfot1321 + F12*F21*dSdfot2121 + F12*F22*dSdfot2221 + F12*F23*dSdfot2321 + F13*F21*dSdfot3121 + F13*F22*dSdfot3221 + F13*F23*dSdfot3321;
        term(6,0) = F21*F31*dSdfot1111 + F21*F32*dSdfot2111 + F21*F33*dSdfot3111 + F22*F31*dSdfot1211 + F22*F32*dSdfot2211 + F22*F33*dSdfot3211 + F23*F31*dSdfot1311 + F23*F32*dSdfot2311 + F23*F33*dSdfot3311;
        term(6,1) = F21*F31*dSdfot1122 + F21*F32*dSdfot2122 + F21*F33*dSdfot3122 + F22*F31*dSdfot1222 + F22*F32*dSdfot2222 + F22*F33*dSdfot3222 + F23*F31*dSdfot1322 + F23*F32*dSdfot2322 + F23*F33*dSdfot3322;
        term(6,2) = F21*F31*dSdfot1133 + F21*F32*dSdfot2133 + F21*F33*dSdfot3133 + F22*F31*dSdfot1233 + F22*F32*dSdfot2233 + F22*F33*dSdfot3233 + F23*F31*dSdfot1333 + F23*F32*dSdfot2333 + F23*F33*dSdfot3333;
        term(6,3) = F21*F31*dSdfot1123 + F21*F32*dSdfot2123 + F21*F33*dSdfot3123 + F22*F31*dSdfot1223 + F22*F32*dSdfot2223 + F22*F33*dSdfot3223 + F23*F31*dSdfot1323 + F23*F32*dSdfot2323 + F23*F33*dSdfot3323;
        term(6,4) = F21*F31*dSdfot1113 + F21*F32*dSdfot2113 + F21*F33*dSdfot3113 + F22*F31*dSdfot1213 + F22*F32*dSdfot2213 + F22*F33*dSdfot3213 + F23*F31*dSdfot1313 + F23*F32*dSdfot2313 + F23*F33*dSdfot3313;
        term(6,5) = F21*F31*dSdfot1112 + F21*F32*dSdfot2112 + F21*F33*dSdfot3112 + F22*F31*dSdfot1212 + F22*F32*dSdfot2212 + F22*F33*dSdfot3212 + F23*F31*dSdfot1312 + F23*F32*dSdfot2312 + F23*F33*dSdfot3312;
        term(6,6) = F21*F31*dSdfot1132 + F21*F32*dSdfot2132 + F21*F33*dSdfot3132 + F22*F31*dSdfot1232 + F22*F32*dSdfot2232 + F22*F33*dSdfot3232 + F23*F31*dSdfot1332 + F23*F32*dSdfot2332 + F23*F33*dSdfot3332;
        term(6,7) = F21*F31*dSdfot1131 + F21*F32*dSdfot2131 + F21*F33*dSdfot3131 + F22*F31*dSdfot1231 + F22*F32*dSdfot2231 + F22*F33*dSdfot3231 + F23*F31*dSdfot1331 + F23*F32*dSdfot2331 + F23*F33*dSdfot3331;
        term(6,8) = F21*F31*dSdfot1121 + F21*F32*dSdfot2121 + F21*F33*dSdfot3121 + F22*F31*dSdfot1221 + F22*F32*dSdfot2221 + F22*F33*dSdfot3221 + F23*F31*dSdfot1321 + F23*F32*dSdfot2321 + F23*F33*dSdfot3321;
        term(7,0) = F11*F31*dSdfot1111 + F11*F32*dSdfot2111 + F11*F33*dSdfot3111 + F12*F31*dSdfot1211 + F12*F32*dSdfot2211 + F12*F33*dSdfot3211 + F13*F31*dSdfot1311 + F13*F32*dSdfot2311 + F13*F33*dSdfot3311;
        term(7,1) = F11*F31*dSdfot1122 + F11*F32*dSdfot2122 + F11*F33*dSdfot3122 + F12*F31*dSdfot1222 + F12*F32*dSdfot2222 + F12*F33*dSdfot3222 + F13*F31*dSdfot1322 + F13*F32*dSdfot2322 + F13*F33*dSdfot3322;
        term(7,2) = F11*F31*dSdfot1133 + F11*F32*dSdfot2133 + F11*F33*dSdfot3133 + F12*F31*dSdfot1233 + F12*F32*dSdfot2233 + F12*F33*dSdfot3233 + F13*F31*dSdfot1333 + F13*F32*dSdfot2333 + F13*F33*dSdfot3333;
        term(7,3) = F11*F31*dSdfot1123 + F11*F32*dSdfot2123 + F11*F33*dSdfot3123 + F12*F31*dSdfot1223 + F12*F32*dSdfot2223 + F12*F33*dSdfot3223 + F13*F31*dSdfot1323 + F13*F32*dSdfot2323 + F13*F33*dSdfot3323;
        term(7,4) = F11*F31*dSdfot1113 + F11*F32*dSdfot2113 + F11*F33*dSdfot3113 + F12*F31*dSdfot1213 + F12*F32*dSdfot2213 + F12*F33*dSdfot3213 + F13*F31*dSdfot1313 + F13*F32*dSdfot2313 + F13*F33*dSdfot3313;
        term(7,5) = F11*F31*dSdfot1112 + F11*F32*dSdfot2112 + F11*F33*dSdfot3112 + F12*F31*dSdfot1212 + F12*F32*dSdfot2212 + F12*F33*dSdfot3212 + F13*F31*dSdfot1312 + F13*F32*dSdfot2312 + F13*F33*dSdfot3312;
        term(7,6) = F11*F31*dSdfot1132 + F11*F32*dSdfot2132 + F11*F33*dSdfot3132 + F12*F31*dSdfot1232 + F12*F32*dSdfot2232 + F12*F33*dSdfot3232 + F13*F31*dSdfot1332 + F13*F32*dSdfot2332 + F13*F33*dSdfot3332;
        term(7,7) = F11*F31*dSdfot1131 + F11*F32*dSdfot2131 + F11*F33*dSdfot3131 + F12*F31*dSdfot1231 + F12*F32*dSdfot2231 + F12*F33*dSdfot3231 + F13*F31*dSdfot1331 + F13*F32*dSdfot2331 + F13*F33*dSdfot3331;
        term(7,8) = F11*F31*dSdfot1121 + F11*F32*dSdfot2121 + F11*F33*dSdfot3121 + F12*F31*dSdfot1221 + F12*F32*dSdfot2221 + F12*F33*dSdfot3221 + F13*F31*dSdfot1321 + F13*F32*dSdfot2321 + F13*F33*dSdfot3321;
        term(8,0) = F11*F21*dSdfot1111 + F11*F22*dSdfot2111 + F11*F23*dSdfot3111 + F12*F21*dSdfot1211 + F12*F22*dSdfot2211 + F12*F23*dSdfot3211 + F13*F21*dSdfot1311 + F13*F22*dSdfot2311 + F13*F23*dSdfot3311;
        term(8,1) = F11*F21*dSdfot1122 + F11*F22*dSdfot2122 + F11*F23*dSdfot3122 + F12*F21*dSdfot1222 + F12*F22*dSdfot2222 + F12*F23*dSdfot3222 + F13*F21*dSdfot1322 + F13*F22*dSdfot2322 + F13*F23*dSdfot3322;
        term(8,2) = F11*F21*dSdfot1133 + F11*F22*dSdfot2133 + F11*F23*dSdfot3133 + F12*F21*dSdfot1233 + F12*F22*dSdfot2233 + F12*F23*dSdfot3233 + F13*F21*dSdfot1333 + F13*F22*dSdfot2333 + F13*F23*dSdfot3333;
        term(8,3) = F11*F21*dSdfot1123 + F11*F22*dSdfot2123 + F11*F23*dSdfot3123 + F12*F21*dSdfot1223 + F12*F22*dSdfot2223 + F12*F23*dSdfot3223 + F13*F21*dSdfot1323 + F13*F22*dSdfot2323 + F13*F23*dSdfot3323;
        term(8,4) = F11*F21*dSdfot1113 + F11*F22*dSdfot2113 + F11*F23*dSdfot3113 + F12*F21*dSdfot1213 + F12*F22*dSdfot2213 + F12*F23*dSdfot3213 + F13*F21*dSdfot1313 + F13*F22*dSdfot2313 + F13*F23*dSdfot3313;
        term(8,5) = F11*F21*dSdfot1112 + F11*F22*dSdfot2112 + F11*F23*dSdfot3112 + F12*F21*dSdfot1212 + F12*F22*dSdfot2212 + F12*F23*dSdfot3212 + F13*F21*dSdfot1312 + F13*F22*dSdfot2312 + F13*F23*dSdfot3312;
        term(8,6) = F11*F21*dSdfot1132 + F11*F22*dSdfot2132 + F11*F23*dSdfot3132 + F12*F21*dSdfot1232 + F12*F22*dSdfot2232 + F12*F23*dSdfot3232 + F13*F21*dSdfot1332 + F13*F22*dSdfot2332 + F13*F23*dSdfot3332;
        term(8,7) = F11*F21*dSdfot1131 + F11*F22*dSdfot2131 + F11*F23*dSdfot3131 + F12*F21*dSdfot1231 + F12*F22*dSdfot2231 + F12*F23*dSdfot3231 + F13*F21*dSdfot1331 + F13*F22*dSdfot2331 + F13*F23*dSdfot3331;
        term(8,8) = F11*F21*dSdfot1121 + F11*F22*dSdfot2121 + F11*F23*dSdfot3121 + F12*F21*dSdfot1221 + F12*F22*dSdfot2221 + F12*F23*dSdfot3221 + F13*F21*dSdfot1321 + F13*F22*dSdfot2321 + F13*F23*dSdfot3321;

        term /= J;
        return;
    }

    void map_fifthot_stress_jacobian(const Matrix_3x3 &F, const Matrix_9x27 &dSdtot, const double &J, Matrix_9x27 &term){
        /*!==================================
        |    map_fifthot_stress_jacobian    |
        =====================================

        Map a fifth order stress jacobian from the reference to the current configuration.

        */
        
        //Extract the deformation gradient
        double F11 = F(0,0);
        double F12 = F(0,1);
        double F13 = F(0,2);
        double F21 = F(1,0);
        double F22 = F(1,1);
        double F23 = F(1,2);
        double F31 = F(2,0);
        double F32 = F(2,1);
        double F33 = F(2,2);

        //Extract the gradient of the stress in the reference configuration
        double dSdtot11111 = dSdtot(0,0);
        double dSdtot11112 = dSdtot(0,5);
        double dSdtot11113 = dSdtot(0,4);
        double dSdtot11121 = dSdtot(0,8);
        double dSdtot11122 = dSdtot(0,1);
        double dSdtot11123 = dSdtot(0,3);
        double dSdtot11131 = dSdtot(0,7);
        double dSdtot11132 = dSdtot(0,6);
        double dSdtot11133 = dSdtot(0,2);
        double dSdtot11211 = dSdtot(0,9);
        double dSdtot11212 = dSdtot(0,14);
        double dSdtot11213 = dSdtot(0,13);
        double dSdtot11221 = dSdtot(0,17);
        double dSdtot11222 = dSdtot(0,10);
        double dSdtot11223 = dSdtot(0,12);
        double dSdtot11231 = dSdtot(0,16);
        double dSdtot11232 = dSdtot(0,15);
        double dSdtot11233 = dSdtot(0,11);
        double dSdtot11311 = dSdtot(0,18);
        double dSdtot11312 = dSdtot(0,23);
        double dSdtot11313 = dSdtot(0,22);
        double dSdtot11321 = dSdtot(0,26);
        double dSdtot11322 = dSdtot(0,19);
        double dSdtot11323 = dSdtot(0,21);
        double dSdtot11331 = dSdtot(0,25);
        double dSdtot11332 = dSdtot(0,24);
        double dSdtot11333 = dSdtot(0,20);
        double dSdtot12111 = dSdtot(5,0);
        double dSdtot12112 = dSdtot(5,5);
        double dSdtot12113 = dSdtot(5,4);
        double dSdtot12121 = dSdtot(5,8);
        double dSdtot12122 = dSdtot(5,1);
        double dSdtot12123 = dSdtot(5,3);
        double dSdtot12131 = dSdtot(5,7);
        double dSdtot12132 = dSdtot(5,6);
        double dSdtot12133 = dSdtot(5,2);
        double dSdtot12211 = dSdtot(5,9);
        double dSdtot12212 = dSdtot(5,14);
        double dSdtot12213 = dSdtot(5,13);
        double dSdtot12221 = dSdtot(5,17);
        double dSdtot12222 = dSdtot(5,10);
        double dSdtot12223 = dSdtot(5,12);
        double dSdtot12231 = dSdtot(5,16);
        double dSdtot12232 = dSdtot(5,15);
        double dSdtot12233 = dSdtot(5,11);
        double dSdtot12311 = dSdtot(5,18);
        double dSdtot12312 = dSdtot(5,23);
        double dSdtot12313 = dSdtot(5,22);
        double dSdtot12321 = dSdtot(5,26);
        double dSdtot12322 = dSdtot(5,19);
        double dSdtot12323 = dSdtot(5,21);
        double dSdtot12331 = dSdtot(5,25);
        double dSdtot12332 = dSdtot(5,24);
        double dSdtot12333 = dSdtot(5,20);
        double dSdtot13111 = dSdtot(4,0);
        double dSdtot13112 = dSdtot(4,5);
        double dSdtot13113 = dSdtot(4,4);
        double dSdtot13121 = dSdtot(4,8);
        double dSdtot13122 = dSdtot(4,1);
        double dSdtot13123 = dSdtot(4,3);
        double dSdtot13131 = dSdtot(4,7);
        double dSdtot13132 = dSdtot(4,6);
        double dSdtot13133 = dSdtot(4,2);
        double dSdtot13211 = dSdtot(4,9);
        double dSdtot13212 = dSdtot(4,14);
        double dSdtot13213 = dSdtot(4,13);
        double dSdtot13221 = dSdtot(4,17);
        double dSdtot13222 = dSdtot(4,10);
        double dSdtot13223 = dSdtot(4,12);
        double dSdtot13231 = dSdtot(4,16);
        double dSdtot13232 = dSdtot(4,15);
        double dSdtot13233 = dSdtot(4,11);
        double dSdtot13311 = dSdtot(4,18);
        double dSdtot13312 = dSdtot(4,23);
        double dSdtot13313 = dSdtot(4,22);
        double dSdtot13321 = dSdtot(4,26);
        double dSdtot13322 = dSdtot(4,19);
        double dSdtot13323 = dSdtot(4,21);
        double dSdtot13331 = dSdtot(4,25);
        double dSdtot13332 = dSdtot(4,24);
        double dSdtot13333 = dSdtot(4,20);
        double dSdtot21111 = dSdtot(8,0);
        double dSdtot21112 = dSdtot(8,5);
        double dSdtot21113 = dSdtot(8,4);
        double dSdtot21121 = dSdtot(8,8);
        double dSdtot21122 = dSdtot(8,1);
        double dSdtot21123 = dSdtot(8,3);
        double dSdtot21131 = dSdtot(8,7);
        double dSdtot21132 = dSdtot(8,6);
        double dSdtot21133 = dSdtot(8,2);
        double dSdtot21211 = dSdtot(8,9);
        double dSdtot21212 = dSdtot(8,14);
        double dSdtot21213 = dSdtot(8,13);
        double dSdtot21221 = dSdtot(8,17);
        double dSdtot21222 = dSdtot(8,10);
        double dSdtot21223 = dSdtot(8,12);
        double dSdtot21231 = dSdtot(8,16);
        double dSdtot21232 = dSdtot(8,15);
        double dSdtot21233 = dSdtot(8,11);
        double dSdtot21311 = dSdtot(8,18);
        double dSdtot21312 = dSdtot(8,23);
        double dSdtot21313 = dSdtot(8,22);
        double dSdtot21321 = dSdtot(8,26);
        double dSdtot21322 = dSdtot(8,19);
        double dSdtot21323 = dSdtot(8,21);
        double dSdtot21331 = dSdtot(8,25);
        double dSdtot21332 = dSdtot(8,24);
        double dSdtot21333 = dSdtot(8,20);
        double dSdtot22111 = dSdtot(1,0);
        double dSdtot22112 = dSdtot(1,5);
        double dSdtot22113 = dSdtot(1,4);
        double dSdtot22121 = dSdtot(1,8);
        double dSdtot22122 = dSdtot(1,1);
        double dSdtot22123 = dSdtot(1,3);
        double dSdtot22131 = dSdtot(1,7);
        double dSdtot22132 = dSdtot(1,6);
        double dSdtot22133 = dSdtot(1,2);
        double dSdtot22211 = dSdtot(1,9);
        double dSdtot22212 = dSdtot(1,14);
        double dSdtot22213 = dSdtot(1,13);
        double dSdtot22221 = dSdtot(1,17);
        double dSdtot22222 = dSdtot(1,10);
        double dSdtot22223 = dSdtot(1,12);
        double dSdtot22231 = dSdtot(1,16);
        double dSdtot22232 = dSdtot(1,15);
        double dSdtot22233 = dSdtot(1,11);
        double dSdtot22311 = dSdtot(1,18);
        double dSdtot22312 = dSdtot(1,23);
        double dSdtot22313 = dSdtot(1,22);
        double dSdtot22321 = dSdtot(1,26);
        double dSdtot22322 = dSdtot(1,19);
        double dSdtot22323 = dSdtot(1,21);
        double dSdtot22331 = dSdtot(1,25);
        double dSdtot22332 = dSdtot(1,24);
        double dSdtot22333 = dSdtot(1,20);
        double dSdtot23111 = dSdtot(3,0);
        double dSdtot23112 = dSdtot(3,5);
        double dSdtot23113 = dSdtot(3,4);
        double dSdtot23121 = dSdtot(3,8);
        double dSdtot23122 = dSdtot(3,1);
        double dSdtot23123 = dSdtot(3,3);
        double dSdtot23131 = dSdtot(3,7);
        double dSdtot23132 = dSdtot(3,6);
        double dSdtot23133 = dSdtot(3,2);
        double dSdtot23211 = dSdtot(3,9);
        double dSdtot23212 = dSdtot(3,14);
        double dSdtot23213 = dSdtot(3,13);
        double dSdtot23221 = dSdtot(3,17);
        double dSdtot23222 = dSdtot(3,10);
        double dSdtot23223 = dSdtot(3,12);
        double dSdtot23231 = dSdtot(3,16);
        double dSdtot23232 = dSdtot(3,15);
        double dSdtot23233 = dSdtot(3,11);
        double dSdtot23311 = dSdtot(3,18);
        double dSdtot23312 = dSdtot(3,23);
        double dSdtot23313 = dSdtot(3,22);
        double dSdtot23321 = dSdtot(3,26);
        double dSdtot23322 = dSdtot(3,19);
        double dSdtot23323 = dSdtot(3,21);
        double dSdtot23331 = dSdtot(3,25);
        double dSdtot23332 = dSdtot(3,24);
        double dSdtot23333 = dSdtot(3,20);
        double dSdtot31111 = dSdtot(7,0);
        double dSdtot31112 = dSdtot(7,5);
        double dSdtot31113 = dSdtot(7,4);
        double dSdtot31121 = dSdtot(7,8);
        double dSdtot31122 = dSdtot(7,1);
        double dSdtot31123 = dSdtot(7,3);
        double dSdtot31131 = dSdtot(7,7);
        double dSdtot31132 = dSdtot(7,6);
        double dSdtot31133 = dSdtot(7,2);
        double dSdtot31211 = dSdtot(7,9);
        double dSdtot31212 = dSdtot(7,14);
        double dSdtot31213 = dSdtot(7,13);
        double dSdtot31221 = dSdtot(7,17);
        double dSdtot31222 = dSdtot(7,10);
        double dSdtot31223 = dSdtot(7,12);
        double dSdtot31231 = dSdtot(7,16);
        double dSdtot31232 = dSdtot(7,15);
        double dSdtot31233 = dSdtot(7,11);
        double dSdtot31311 = dSdtot(7,18);
        double dSdtot31312 = dSdtot(7,23);
        double dSdtot31313 = dSdtot(7,22);
        double dSdtot31321 = dSdtot(7,26);
        double dSdtot31322 = dSdtot(7,19);
        double dSdtot31323 = dSdtot(7,21);
        double dSdtot31331 = dSdtot(7,25);
        double dSdtot31332 = dSdtot(7,24);
        double dSdtot31333 = dSdtot(7,20);
        double dSdtot32111 = dSdtot(6,0);
        double dSdtot32112 = dSdtot(6,5);
        double dSdtot32113 = dSdtot(6,4);
        double dSdtot32121 = dSdtot(6,8);
        double dSdtot32122 = dSdtot(6,1);
        double dSdtot32123 = dSdtot(6,3);
        double dSdtot32131 = dSdtot(6,7);
        double dSdtot32132 = dSdtot(6,6);
        double dSdtot32133 = dSdtot(6,2);
        double dSdtot32211 = dSdtot(6,9);
        double dSdtot32212 = dSdtot(6,14);
        double dSdtot32213 = dSdtot(6,13);
        double dSdtot32221 = dSdtot(6,17);
        double dSdtot32222 = dSdtot(6,10);
        double dSdtot32223 = dSdtot(6,12);
        double dSdtot32231 = dSdtot(6,16);
        double dSdtot32232 = dSdtot(6,15);
        double dSdtot32233 = dSdtot(6,11);
        double dSdtot32311 = dSdtot(6,18);
        double dSdtot32312 = dSdtot(6,23);
        double dSdtot32313 = dSdtot(6,22);
        double dSdtot32321 = dSdtot(6,26);
        double dSdtot32322 = dSdtot(6,19);
        double dSdtot32323 = dSdtot(6,21);
        double dSdtot32331 = dSdtot(6,25);
        double dSdtot32332 = dSdtot(6,24);
        double dSdtot32333 = dSdtot(6,20);
        double dSdtot33111 = dSdtot(2,0);
        double dSdtot33112 = dSdtot(2,5);
        double dSdtot33113 = dSdtot(2,4);
        double dSdtot33121 = dSdtot(2,8);
        double dSdtot33122 = dSdtot(2,1);
        double dSdtot33123 = dSdtot(2,3);
        double dSdtot33131 = dSdtot(2,7);
        double dSdtot33132 = dSdtot(2,6);
        double dSdtot33133 = dSdtot(2,2);
        double dSdtot33211 = dSdtot(2,9);
        double dSdtot33212 = dSdtot(2,14);
        double dSdtot33213 = dSdtot(2,13);
        double dSdtot33221 = dSdtot(2,17);
        double dSdtot33222 = dSdtot(2,10);
        double dSdtot33223 = dSdtot(2,12);
        double dSdtot33231 = dSdtot(2,16);
        double dSdtot33232 = dSdtot(2,15);
        double dSdtot33233 = dSdtot(2,11);
        double dSdtot33311 = dSdtot(2,18);
        double dSdtot33312 = dSdtot(2,23);
        double dSdtot33313 = dSdtot(2,22);
        double dSdtot33321 = dSdtot(2,26);
        double dSdtot33322 = dSdtot(2,19);
        double dSdtot33323 = dSdtot(2,21);
        double dSdtot33331 = dSdtot(2,25);
        double dSdtot33332 = dSdtot(2,24);
        double dSdtot33333 = dSdtot(2,20);

        //Map the stress gradient w.r.t. the gradient of chi
        term(0,0) = pow(F11, 2)*dSdtot11111 + F11*F12*dSdtot12111 + F11*F12*dSdtot21111 + F11*F13*dSdtot13111 + F11*F13*dSdtot31111 + pow(F12, 2)*dSdtot22111 + F12*F13*dSdtot23111 + F12*F13*dSdtot32111 + pow(F13, 2)*dSdtot33111;
        term(0,1) = pow(F11, 2)*dSdtot11122 + F11*F12*dSdtot12122 + F11*F12*dSdtot21122 + F11*F13*dSdtot13122 + F11*F13*dSdtot31122 + pow(F12, 2)*dSdtot22122 + F12*F13*dSdtot23122 + F12*F13*dSdtot32122 + pow(F13, 2)*dSdtot33122;
        term(0,2) = pow(F11, 2)*dSdtot11133 + F11*F12*dSdtot12133 + F11*F12*dSdtot21133 + F11*F13*dSdtot13133 + F11*F13*dSdtot31133 + pow(F12, 2)*dSdtot22133 + F12*F13*dSdtot23133 + F12*F13*dSdtot32133 + pow(F13, 2)*dSdtot33133;
        term(0,3) = pow(F11, 2)*dSdtot11123 + F11*F12*dSdtot12123 + F11*F12*dSdtot21123 + F11*F13*dSdtot13123 + F11*F13*dSdtot31123 + pow(F12, 2)*dSdtot22123 + F12*F13*dSdtot23123 + F12*F13*dSdtot32123 + pow(F13, 2)*dSdtot33123;
        term(0,4) = pow(F11, 2)*dSdtot11113 + F11*F12*dSdtot12113 + F11*F12*dSdtot21113 + F11*F13*dSdtot13113 + F11*F13*dSdtot31113 + pow(F12, 2)*dSdtot22113 + F12*F13*dSdtot23113 + F12*F13*dSdtot32113 + pow(F13, 2)*dSdtot33113;
        term(0,5) = pow(F11, 2)*dSdtot11112 + F11*F12*dSdtot12112 + F11*F12*dSdtot21112 + F11*F13*dSdtot13112 + F11*F13*dSdtot31112 + pow(F12, 2)*dSdtot22112 + F12*F13*dSdtot23112 + F12*F13*dSdtot32112 + pow(F13, 2)*dSdtot33112;
        term(0,6) = pow(F11, 2)*dSdtot11132 + F11*F12*dSdtot12132 + F11*F12*dSdtot21132 + F11*F13*dSdtot13132 + F11*F13*dSdtot31132 + pow(F12, 2)*dSdtot22132 + F12*F13*dSdtot23132 + F12*F13*dSdtot32132 + pow(F13, 2)*dSdtot33132;
        term(0,7) = pow(F11, 2)*dSdtot11131 + F11*F12*dSdtot12131 + F11*F12*dSdtot21131 + F11*F13*dSdtot13131 + F11*F13*dSdtot31131 + pow(F12, 2)*dSdtot22131 + F12*F13*dSdtot23131 + F12*F13*dSdtot32131 + pow(F13, 2)*dSdtot33131;
        term(0,8) = pow(F11, 2)*dSdtot11121 + F11*F12*dSdtot12121 + F11*F12*dSdtot21121 + F11*F13*dSdtot13121 + F11*F13*dSdtot31121 + pow(F12, 2)*dSdtot22121 + F12*F13*dSdtot23121 + F12*F13*dSdtot32121 + pow(F13, 2)*dSdtot33121;
        term(0,9) = pow(F11, 2)*dSdtot11211 + F11*F12*dSdtot12211 + F11*F12*dSdtot21211 + F11*F13*dSdtot13211 + F11*F13*dSdtot31211 + pow(F12, 2)*dSdtot22211 + F12*F13*dSdtot23211 + F12*F13*dSdtot32211 + pow(F13, 2)*dSdtot33211;
        term(0,10) = pow(F11, 2)*dSdtot11222 + F11*F12*dSdtot12222 + F11*F12*dSdtot21222 + F11*F13*dSdtot13222 + F11*F13*dSdtot31222 + pow(F12, 2)*dSdtot22222 + F12*F13*dSdtot23222 + F12*F13*dSdtot32222 + pow(F13, 2)*dSdtot33222;
        term(0,11) = pow(F11, 2)*dSdtot11233 + F11*F12*dSdtot12233 + F11*F12*dSdtot21233 + F11*F13*dSdtot13233 + F11*F13*dSdtot31233 + pow(F12, 2)*dSdtot22233 + F12*F13*dSdtot23233 + F12*F13*dSdtot32233 + pow(F13, 2)*dSdtot33233;
        term(0,12) = pow(F11, 2)*dSdtot11223 + F11*F12*dSdtot12223 + F11*F12*dSdtot21223 + F11*F13*dSdtot13223 + F11*F13*dSdtot31223 + pow(F12, 2)*dSdtot22223 + F12*F13*dSdtot23223 + F12*F13*dSdtot32223 + pow(F13, 2)*dSdtot33223;
        term(0,13) = pow(F11, 2)*dSdtot11213 + F11*F12*dSdtot12213 + F11*F12*dSdtot21213 + F11*F13*dSdtot13213 + F11*F13*dSdtot31213 + pow(F12, 2)*dSdtot22213 + F12*F13*dSdtot23213 + F12*F13*dSdtot32213 + pow(F13, 2)*dSdtot33213;
        term(0,14) = pow(F11, 2)*dSdtot11212 + F11*F12*dSdtot12212 + F11*F12*dSdtot21212 + F11*F13*dSdtot13212 + F11*F13*dSdtot31212 + pow(F12, 2)*dSdtot22212 + F12*F13*dSdtot23212 + F12*F13*dSdtot32212 + pow(F13, 2)*dSdtot33212;
        term(0,15) = pow(F11, 2)*dSdtot11232 + F11*F12*dSdtot12232 + F11*F12*dSdtot21232 + F11*F13*dSdtot13232 + F11*F13*dSdtot31232 + pow(F12, 2)*dSdtot22232 + F12*F13*dSdtot23232 + F12*F13*dSdtot32232 + pow(F13, 2)*dSdtot33232;
        term(0,16) = pow(F11, 2)*dSdtot11231 + F11*F12*dSdtot12231 + F11*F12*dSdtot21231 + F11*F13*dSdtot13231 + F11*F13*dSdtot31231 + pow(F12, 2)*dSdtot22231 + F12*F13*dSdtot23231 + F12*F13*dSdtot32231 + pow(F13, 2)*dSdtot33231;
        term(0,17) = pow(F11, 2)*dSdtot11221 + F11*F12*dSdtot12221 + F11*F12*dSdtot21221 + F11*F13*dSdtot13221 + F11*F13*dSdtot31221 + pow(F12, 2)*dSdtot22221 + F12*F13*dSdtot23221 + F12*F13*dSdtot32221 + pow(F13, 2)*dSdtot33221;
        term(0,18) = pow(F11, 2)*dSdtot11311 + F11*F12*dSdtot12311 + F11*F12*dSdtot21311 + F11*F13*dSdtot13311 + F11*F13*dSdtot31311 + pow(F12, 2)*dSdtot22311 + F12*F13*dSdtot23311 + F12*F13*dSdtot32311 + pow(F13, 2)*dSdtot33311;
        term(0,19) = pow(F11, 2)*dSdtot11322 + F11*F12*dSdtot12322 + F11*F12*dSdtot21322 + F11*F13*dSdtot13322 + F11*F13*dSdtot31322 + pow(F12, 2)*dSdtot22322 + F12*F13*dSdtot23322 + F12*F13*dSdtot32322 + pow(F13, 2)*dSdtot33322;
        term(0,20) = pow(F11, 2)*dSdtot11333 + F11*F12*dSdtot12333 + F11*F12*dSdtot21333 + F11*F13*dSdtot13333 + F11*F13*dSdtot31333 + pow(F12, 2)*dSdtot22333 + F12*F13*dSdtot23333 + F12*F13*dSdtot32333 + pow(F13, 2)*dSdtot33333;
        term(0,21) = pow(F11, 2)*dSdtot11323 + F11*F12*dSdtot12323 + F11*F12*dSdtot21323 + F11*F13*dSdtot13323 + F11*F13*dSdtot31323 + pow(F12, 2)*dSdtot22323 + F12*F13*dSdtot23323 + F12*F13*dSdtot32323 + pow(F13, 2)*dSdtot33323;
        term(0,22) = pow(F11, 2)*dSdtot11313 + F11*F12*dSdtot12313 + F11*F12*dSdtot21313 + F11*F13*dSdtot13313 + F11*F13*dSdtot31313 + pow(F12, 2)*dSdtot22313 + F12*F13*dSdtot23313 + F12*F13*dSdtot32313 + pow(F13, 2)*dSdtot33313;
        term(0,23) = pow(F11, 2)*dSdtot11312 + F11*F12*dSdtot12312 + F11*F12*dSdtot21312 + F11*F13*dSdtot13312 + F11*F13*dSdtot31312 + pow(F12, 2)*dSdtot22312 + F12*F13*dSdtot23312 + F12*F13*dSdtot32312 + pow(F13, 2)*dSdtot33312;
        term(0,24) = pow(F11, 2)*dSdtot11332 + F11*F12*dSdtot12332 + F11*F12*dSdtot21332 + F11*F13*dSdtot13332 + F11*F13*dSdtot31332 + pow(F12, 2)*dSdtot22332 + F12*F13*dSdtot23332 + F12*F13*dSdtot32332 + pow(F13, 2)*dSdtot33332;
        term(0,25) = pow(F11, 2)*dSdtot11331 + F11*F12*dSdtot12331 + F11*F12*dSdtot21331 + F11*F13*dSdtot13331 + F11*F13*dSdtot31331 + pow(F12, 2)*dSdtot22331 + F12*F13*dSdtot23331 + F12*F13*dSdtot32331 + pow(F13, 2)*dSdtot33331;
        term(0,26) = pow(F11, 2)*dSdtot11321 + F11*F12*dSdtot12321 + F11*F12*dSdtot21321 + F11*F13*dSdtot13321 + F11*F13*dSdtot31321 + pow(F12, 2)*dSdtot22321 + F12*F13*dSdtot23321 + F12*F13*dSdtot32321 + pow(F13, 2)*dSdtot33321;
        term(1,0) = pow(F21, 2)*dSdtot11111 + F21*F22*dSdtot12111 + F21*F22*dSdtot21111 + F21*F23*dSdtot13111 + F21*F23*dSdtot31111 + pow(F22, 2)*dSdtot22111 + F22*F23*dSdtot23111 + F22*F23*dSdtot32111 + pow(F23, 2)*dSdtot33111;
        term(1,1) = pow(F21, 2)*dSdtot11122 + F21*F22*dSdtot12122 + F21*F22*dSdtot21122 + F21*F23*dSdtot13122 + F21*F23*dSdtot31122 + pow(F22, 2)*dSdtot22122 + F22*F23*dSdtot23122 + F22*F23*dSdtot32122 + pow(F23, 2)*dSdtot33122;
        term(1,2) = pow(F21, 2)*dSdtot11133 + F21*F22*dSdtot12133 + F21*F22*dSdtot21133 + F21*F23*dSdtot13133 + F21*F23*dSdtot31133 + pow(F22, 2)*dSdtot22133 + F22*F23*dSdtot23133 + F22*F23*dSdtot32133 + pow(F23, 2)*dSdtot33133;
        term(1,3) = pow(F21, 2)*dSdtot11123 + F21*F22*dSdtot12123 + F21*F22*dSdtot21123 + F21*F23*dSdtot13123 + F21*F23*dSdtot31123 + pow(F22, 2)*dSdtot22123 + F22*F23*dSdtot23123 + F22*F23*dSdtot32123 + pow(F23, 2)*dSdtot33123;
        term(1,4) = pow(F21, 2)*dSdtot11113 + F21*F22*dSdtot12113 + F21*F22*dSdtot21113 + F21*F23*dSdtot13113 + F21*F23*dSdtot31113 + pow(F22, 2)*dSdtot22113 + F22*F23*dSdtot23113 + F22*F23*dSdtot32113 + pow(F23, 2)*dSdtot33113;
        term(1,5) = pow(F21, 2)*dSdtot11112 + F21*F22*dSdtot12112 + F21*F22*dSdtot21112 + F21*F23*dSdtot13112 + F21*F23*dSdtot31112 + pow(F22, 2)*dSdtot22112 + F22*F23*dSdtot23112 + F22*F23*dSdtot32112 + pow(F23, 2)*dSdtot33112;
        term(1,6) = pow(F21, 2)*dSdtot11132 + F21*F22*dSdtot12132 + F21*F22*dSdtot21132 + F21*F23*dSdtot13132 + F21*F23*dSdtot31132 + pow(F22, 2)*dSdtot22132 + F22*F23*dSdtot23132 + F22*F23*dSdtot32132 + pow(F23, 2)*dSdtot33132;
        term(1,7) = pow(F21, 2)*dSdtot11131 + F21*F22*dSdtot12131 + F21*F22*dSdtot21131 + F21*F23*dSdtot13131 + F21*F23*dSdtot31131 + pow(F22, 2)*dSdtot22131 + F22*F23*dSdtot23131 + F22*F23*dSdtot32131 + pow(F23, 2)*dSdtot33131;
        term(1,8) = pow(F21, 2)*dSdtot11121 + F21*F22*dSdtot12121 + F21*F22*dSdtot21121 + F21*F23*dSdtot13121 + F21*F23*dSdtot31121 + pow(F22, 2)*dSdtot22121 + F22*F23*dSdtot23121 + F22*F23*dSdtot32121 + pow(F23, 2)*dSdtot33121;
        term(1,9) = pow(F21, 2)*dSdtot11211 + F21*F22*dSdtot12211 + F21*F22*dSdtot21211 + F21*F23*dSdtot13211 + F21*F23*dSdtot31211 + pow(F22, 2)*dSdtot22211 + F22*F23*dSdtot23211 + F22*F23*dSdtot32211 + pow(F23, 2)*dSdtot33211;
        term(1,10) = pow(F21, 2)*dSdtot11222 + F21*F22*dSdtot12222 + F21*F22*dSdtot21222 + F21*F23*dSdtot13222 + F21*F23*dSdtot31222 + pow(F22, 2)*dSdtot22222 + F22*F23*dSdtot23222 + F22*F23*dSdtot32222 + pow(F23, 2)*dSdtot33222;
        term(1,11) = pow(F21, 2)*dSdtot11233 + F21*F22*dSdtot12233 + F21*F22*dSdtot21233 + F21*F23*dSdtot13233 + F21*F23*dSdtot31233 + pow(F22, 2)*dSdtot22233 + F22*F23*dSdtot23233 + F22*F23*dSdtot32233 + pow(F23, 2)*dSdtot33233;
        term(1,12) = pow(F21, 2)*dSdtot11223 + F21*F22*dSdtot12223 + F21*F22*dSdtot21223 + F21*F23*dSdtot13223 + F21*F23*dSdtot31223 + pow(F22, 2)*dSdtot22223 + F22*F23*dSdtot23223 + F22*F23*dSdtot32223 + pow(F23, 2)*dSdtot33223;
        term(1,13) = pow(F21, 2)*dSdtot11213 + F21*F22*dSdtot12213 + F21*F22*dSdtot21213 + F21*F23*dSdtot13213 + F21*F23*dSdtot31213 + pow(F22, 2)*dSdtot22213 + F22*F23*dSdtot23213 + F22*F23*dSdtot32213 + pow(F23, 2)*dSdtot33213;
        term(1,14) = pow(F21, 2)*dSdtot11212 + F21*F22*dSdtot12212 + F21*F22*dSdtot21212 + F21*F23*dSdtot13212 + F21*F23*dSdtot31212 + pow(F22, 2)*dSdtot22212 + F22*F23*dSdtot23212 + F22*F23*dSdtot32212 + pow(F23, 2)*dSdtot33212;
        term(1,15) = pow(F21, 2)*dSdtot11232 + F21*F22*dSdtot12232 + F21*F22*dSdtot21232 + F21*F23*dSdtot13232 + F21*F23*dSdtot31232 + pow(F22, 2)*dSdtot22232 + F22*F23*dSdtot23232 + F22*F23*dSdtot32232 + pow(F23, 2)*dSdtot33232;
        term(1,16) = pow(F21, 2)*dSdtot11231 + F21*F22*dSdtot12231 + F21*F22*dSdtot21231 + F21*F23*dSdtot13231 + F21*F23*dSdtot31231 + pow(F22, 2)*dSdtot22231 + F22*F23*dSdtot23231 + F22*F23*dSdtot32231 + pow(F23, 2)*dSdtot33231;
        term(1,17) = pow(F21, 2)*dSdtot11221 + F21*F22*dSdtot12221 + F21*F22*dSdtot21221 + F21*F23*dSdtot13221 + F21*F23*dSdtot31221 + pow(F22, 2)*dSdtot22221 + F22*F23*dSdtot23221 + F22*F23*dSdtot32221 + pow(F23, 2)*dSdtot33221;
        term(1,18) = pow(F21, 2)*dSdtot11311 + F21*F22*dSdtot12311 + F21*F22*dSdtot21311 + F21*F23*dSdtot13311 + F21*F23*dSdtot31311 + pow(F22, 2)*dSdtot22311 + F22*F23*dSdtot23311 + F22*F23*dSdtot32311 + pow(F23, 2)*dSdtot33311;
        term(1,19) = pow(F21, 2)*dSdtot11322 + F21*F22*dSdtot12322 + F21*F22*dSdtot21322 + F21*F23*dSdtot13322 + F21*F23*dSdtot31322 + pow(F22, 2)*dSdtot22322 + F22*F23*dSdtot23322 + F22*F23*dSdtot32322 + pow(F23, 2)*dSdtot33322;
        term(1,20) = pow(F21, 2)*dSdtot11333 + F21*F22*dSdtot12333 + F21*F22*dSdtot21333 + F21*F23*dSdtot13333 + F21*F23*dSdtot31333 + pow(F22, 2)*dSdtot22333 + F22*F23*dSdtot23333 + F22*F23*dSdtot32333 + pow(F23, 2)*dSdtot33333;
        term(1,21) = pow(F21, 2)*dSdtot11323 + F21*F22*dSdtot12323 + F21*F22*dSdtot21323 + F21*F23*dSdtot13323 + F21*F23*dSdtot31323 + pow(F22, 2)*dSdtot22323 + F22*F23*dSdtot23323 + F22*F23*dSdtot32323 + pow(F23, 2)*dSdtot33323;
        term(1,22) = pow(F21, 2)*dSdtot11313 + F21*F22*dSdtot12313 + F21*F22*dSdtot21313 + F21*F23*dSdtot13313 + F21*F23*dSdtot31313 + pow(F22, 2)*dSdtot22313 + F22*F23*dSdtot23313 + F22*F23*dSdtot32313 + pow(F23, 2)*dSdtot33313;
        term(1,23) = pow(F21, 2)*dSdtot11312 + F21*F22*dSdtot12312 + F21*F22*dSdtot21312 + F21*F23*dSdtot13312 + F21*F23*dSdtot31312 + pow(F22, 2)*dSdtot22312 + F22*F23*dSdtot23312 + F22*F23*dSdtot32312 + pow(F23, 2)*dSdtot33312;
        term(1,24) = pow(F21, 2)*dSdtot11332 + F21*F22*dSdtot12332 + F21*F22*dSdtot21332 + F21*F23*dSdtot13332 + F21*F23*dSdtot31332 + pow(F22, 2)*dSdtot22332 + F22*F23*dSdtot23332 + F22*F23*dSdtot32332 + pow(F23, 2)*dSdtot33332;
        term(1,25) = pow(F21, 2)*dSdtot11331 + F21*F22*dSdtot12331 + F21*F22*dSdtot21331 + F21*F23*dSdtot13331 + F21*F23*dSdtot31331 + pow(F22, 2)*dSdtot22331 + F22*F23*dSdtot23331 + F22*F23*dSdtot32331 + pow(F23, 2)*dSdtot33331;
        term(1,26) = pow(F21, 2)*dSdtot11321 + F21*F22*dSdtot12321 + F21*F22*dSdtot21321 + F21*F23*dSdtot13321 + F21*F23*dSdtot31321 + pow(F22, 2)*dSdtot22321 + F22*F23*dSdtot23321 + F22*F23*dSdtot32321 + pow(F23, 2)*dSdtot33321;
        term(2,0) = pow(F31, 2)*dSdtot11111 + F31*F32*dSdtot12111 + F31*F32*dSdtot21111 + F31*F33*dSdtot13111 + F31*F33*dSdtot31111 + pow(F32, 2)*dSdtot22111 + F32*F33*dSdtot23111 + F32*F33*dSdtot32111 + pow(F33, 2)*dSdtot33111;
        term(2,1) = pow(F31, 2)*dSdtot11122 + F31*F32*dSdtot12122 + F31*F32*dSdtot21122 + F31*F33*dSdtot13122 + F31*F33*dSdtot31122 + pow(F32, 2)*dSdtot22122 + F32*F33*dSdtot23122 + F32*F33*dSdtot32122 + pow(F33, 2)*dSdtot33122;
        term(2,2) = pow(F31, 2)*dSdtot11133 + F31*F32*dSdtot12133 + F31*F32*dSdtot21133 + F31*F33*dSdtot13133 + F31*F33*dSdtot31133 + pow(F32, 2)*dSdtot22133 + F32*F33*dSdtot23133 + F32*F33*dSdtot32133 + pow(F33, 2)*dSdtot33133;
        term(2,3) = pow(F31, 2)*dSdtot11123 + F31*F32*dSdtot12123 + F31*F32*dSdtot21123 + F31*F33*dSdtot13123 + F31*F33*dSdtot31123 + pow(F32, 2)*dSdtot22123 + F32*F33*dSdtot23123 + F32*F33*dSdtot32123 + pow(F33, 2)*dSdtot33123;
        term(2,4) = pow(F31, 2)*dSdtot11113 + F31*F32*dSdtot12113 + F31*F32*dSdtot21113 + F31*F33*dSdtot13113 + F31*F33*dSdtot31113 + pow(F32, 2)*dSdtot22113 + F32*F33*dSdtot23113 + F32*F33*dSdtot32113 + pow(F33, 2)*dSdtot33113;
        term(2,5) = pow(F31, 2)*dSdtot11112 + F31*F32*dSdtot12112 + F31*F32*dSdtot21112 + F31*F33*dSdtot13112 + F31*F33*dSdtot31112 + pow(F32, 2)*dSdtot22112 + F32*F33*dSdtot23112 + F32*F33*dSdtot32112 + pow(F33, 2)*dSdtot33112;
        term(2,6) = pow(F31, 2)*dSdtot11132 + F31*F32*dSdtot12132 + F31*F32*dSdtot21132 + F31*F33*dSdtot13132 + F31*F33*dSdtot31132 + pow(F32, 2)*dSdtot22132 + F32*F33*dSdtot23132 + F32*F33*dSdtot32132 + pow(F33, 2)*dSdtot33132;
        term(2,7) = pow(F31, 2)*dSdtot11131 + F31*F32*dSdtot12131 + F31*F32*dSdtot21131 + F31*F33*dSdtot13131 + F31*F33*dSdtot31131 + pow(F32, 2)*dSdtot22131 + F32*F33*dSdtot23131 + F32*F33*dSdtot32131 + pow(F33, 2)*dSdtot33131;
        term(2,8) = pow(F31, 2)*dSdtot11121 + F31*F32*dSdtot12121 + F31*F32*dSdtot21121 + F31*F33*dSdtot13121 + F31*F33*dSdtot31121 + pow(F32, 2)*dSdtot22121 + F32*F33*dSdtot23121 + F32*F33*dSdtot32121 + pow(F33, 2)*dSdtot33121;
        term(2,9) = pow(F31, 2)*dSdtot11211 + F31*F32*dSdtot12211 + F31*F32*dSdtot21211 + F31*F33*dSdtot13211 + F31*F33*dSdtot31211 + pow(F32, 2)*dSdtot22211 + F32*F33*dSdtot23211 + F32*F33*dSdtot32211 + pow(F33, 2)*dSdtot33211;
        term(2,10) = pow(F31, 2)*dSdtot11222 + F31*F32*dSdtot12222 + F31*F32*dSdtot21222 + F31*F33*dSdtot13222 + F31*F33*dSdtot31222 + pow(F32, 2)*dSdtot22222 + F32*F33*dSdtot23222 + F32*F33*dSdtot32222 + pow(F33, 2)*dSdtot33222;
        term(2,11) = pow(F31, 2)*dSdtot11233 + F31*F32*dSdtot12233 + F31*F32*dSdtot21233 + F31*F33*dSdtot13233 + F31*F33*dSdtot31233 + pow(F32, 2)*dSdtot22233 + F32*F33*dSdtot23233 + F32*F33*dSdtot32233 + pow(F33, 2)*dSdtot33233;
        term(2,12) = pow(F31, 2)*dSdtot11223 + F31*F32*dSdtot12223 + F31*F32*dSdtot21223 + F31*F33*dSdtot13223 + F31*F33*dSdtot31223 + pow(F32, 2)*dSdtot22223 + F32*F33*dSdtot23223 + F32*F33*dSdtot32223 + pow(F33, 2)*dSdtot33223;
        term(2,13) = pow(F31, 2)*dSdtot11213 + F31*F32*dSdtot12213 + F31*F32*dSdtot21213 + F31*F33*dSdtot13213 + F31*F33*dSdtot31213 + pow(F32, 2)*dSdtot22213 + F32*F33*dSdtot23213 + F32*F33*dSdtot32213 + pow(F33, 2)*dSdtot33213;
        term(2,14) = pow(F31, 2)*dSdtot11212 + F31*F32*dSdtot12212 + F31*F32*dSdtot21212 + F31*F33*dSdtot13212 + F31*F33*dSdtot31212 + pow(F32, 2)*dSdtot22212 + F32*F33*dSdtot23212 + F32*F33*dSdtot32212 + pow(F33, 2)*dSdtot33212;
        term(2,15) = pow(F31, 2)*dSdtot11232 + F31*F32*dSdtot12232 + F31*F32*dSdtot21232 + F31*F33*dSdtot13232 + F31*F33*dSdtot31232 + pow(F32, 2)*dSdtot22232 + F32*F33*dSdtot23232 + F32*F33*dSdtot32232 + pow(F33, 2)*dSdtot33232;
        term(2,16) = pow(F31, 2)*dSdtot11231 + F31*F32*dSdtot12231 + F31*F32*dSdtot21231 + F31*F33*dSdtot13231 + F31*F33*dSdtot31231 + pow(F32, 2)*dSdtot22231 + F32*F33*dSdtot23231 + F32*F33*dSdtot32231 + pow(F33, 2)*dSdtot33231;
        term(2,17) = pow(F31, 2)*dSdtot11221 + F31*F32*dSdtot12221 + F31*F32*dSdtot21221 + F31*F33*dSdtot13221 + F31*F33*dSdtot31221 + pow(F32, 2)*dSdtot22221 + F32*F33*dSdtot23221 + F32*F33*dSdtot32221 + pow(F33, 2)*dSdtot33221;
        term(2,18) = pow(F31, 2)*dSdtot11311 + F31*F32*dSdtot12311 + F31*F32*dSdtot21311 + F31*F33*dSdtot13311 + F31*F33*dSdtot31311 + pow(F32, 2)*dSdtot22311 + F32*F33*dSdtot23311 + F32*F33*dSdtot32311 + pow(F33, 2)*dSdtot33311;
        term(2,19) = pow(F31, 2)*dSdtot11322 + F31*F32*dSdtot12322 + F31*F32*dSdtot21322 + F31*F33*dSdtot13322 + F31*F33*dSdtot31322 + pow(F32, 2)*dSdtot22322 + F32*F33*dSdtot23322 + F32*F33*dSdtot32322 + pow(F33, 2)*dSdtot33322;
        term(2,20) = pow(F31, 2)*dSdtot11333 + F31*F32*dSdtot12333 + F31*F32*dSdtot21333 + F31*F33*dSdtot13333 + F31*F33*dSdtot31333 + pow(F32, 2)*dSdtot22333 + F32*F33*dSdtot23333 + F32*F33*dSdtot32333 + pow(F33, 2)*dSdtot33333;
        term(2,21) = pow(F31, 2)*dSdtot11323 + F31*F32*dSdtot12323 + F31*F32*dSdtot21323 + F31*F33*dSdtot13323 + F31*F33*dSdtot31323 + pow(F32, 2)*dSdtot22323 + F32*F33*dSdtot23323 + F32*F33*dSdtot32323 + pow(F33, 2)*dSdtot33323;
        term(2,22) = pow(F31, 2)*dSdtot11313 + F31*F32*dSdtot12313 + F31*F32*dSdtot21313 + F31*F33*dSdtot13313 + F31*F33*dSdtot31313 + pow(F32, 2)*dSdtot22313 + F32*F33*dSdtot23313 + F32*F33*dSdtot32313 + pow(F33, 2)*dSdtot33313;
        term(2,23) = pow(F31, 2)*dSdtot11312 + F31*F32*dSdtot12312 + F31*F32*dSdtot21312 + F31*F33*dSdtot13312 + F31*F33*dSdtot31312 + pow(F32, 2)*dSdtot22312 + F32*F33*dSdtot23312 + F32*F33*dSdtot32312 + pow(F33, 2)*dSdtot33312;
        term(2,24) = pow(F31, 2)*dSdtot11332 + F31*F32*dSdtot12332 + F31*F32*dSdtot21332 + F31*F33*dSdtot13332 + F31*F33*dSdtot31332 + pow(F32, 2)*dSdtot22332 + F32*F33*dSdtot23332 + F32*F33*dSdtot32332 + pow(F33, 2)*dSdtot33332;
        term(2,25) = pow(F31, 2)*dSdtot11331 + F31*F32*dSdtot12331 + F31*F32*dSdtot21331 + F31*F33*dSdtot13331 + F31*F33*dSdtot31331 + pow(F32, 2)*dSdtot22331 + F32*F33*dSdtot23331 + F32*F33*dSdtot32331 + pow(F33, 2)*dSdtot33331;
        term(2,26) = pow(F31, 2)*dSdtot11321 + F31*F32*dSdtot12321 + F31*F32*dSdtot21321 + F31*F33*dSdtot13321 + F31*F33*dSdtot31321 + pow(F32, 2)*dSdtot22321 + F32*F33*dSdtot23321 + F32*F33*dSdtot32321 + pow(F33, 2)*dSdtot33321;
        term(3,0) = F21*F31*dSdtot11111 + F21*F32*dSdtot12111 + F21*F33*dSdtot13111 + F22*F31*dSdtot21111 + F22*F32*dSdtot22111 + F22*F33*dSdtot23111 + F23*F31*dSdtot31111 + F23*F32*dSdtot32111 + F23*F33*dSdtot33111;
        term(3,1) = F21*F31*dSdtot11122 + F21*F32*dSdtot12122 + F21*F33*dSdtot13122 + F22*F31*dSdtot21122 + F22*F32*dSdtot22122 + F22*F33*dSdtot23122 + F23*F31*dSdtot31122 + F23*F32*dSdtot32122 + F23*F33*dSdtot33122;
        term(3,2) = F21*F31*dSdtot11133 + F21*F32*dSdtot12133 + F21*F33*dSdtot13133 + F22*F31*dSdtot21133 + F22*F32*dSdtot22133 + F22*F33*dSdtot23133 + F23*F31*dSdtot31133 + F23*F32*dSdtot32133 + F23*F33*dSdtot33133;
        term(3,3) = F21*F31*dSdtot11123 + F21*F32*dSdtot12123 + F21*F33*dSdtot13123 + F22*F31*dSdtot21123 + F22*F32*dSdtot22123 + F22*F33*dSdtot23123 + F23*F31*dSdtot31123 + F23*F32*dSdtot32123 + F23*F33*dSdtot33123;
        term(3,4) = F21*F31*dSdtot11113 + F21*F32*dSdtot12113 + F21*F33*dSdtot13113 + F22*F31*dSdtot21113 + F22*F32*dSdtot22113 + F22*F33*dSdtot23113 + F23*F31*dSdtot31113 + F23*F32*dSdtot32113 + F23*F33*dSdtot33113;
        term(3,5) = F21*F31*dSdtot11112 + F21*F32*dSdtot12112 + F21*F33*dSdtot13112 + F22*F31*dSdtot21112 + F22*F32*dSdtot22112 + F22*F33*dSdtot23112 + F23*F31*dSdtot31112 + F23*F32*dSdtot32112 + F23*F33*dSdtot33112;
        term(3,6) = F21*F31*dSdtot11132 + F21*F32*dSdtot12132 + F21*F33*dSdtot13132 + F22*F31*dSdtot21132 + F22*F32*dSdtot22132 + F22*F33*dSdtot23132 + F23*F31*dSdtot31132 + F23*F32*dSdtot32132 + F23*F33*dSdtot33132;
        term(3,7) = F21*F31*dSdtot11131 + F21*F32*dSdtot12131 + F21*F33*dSdtot13131 + F22*F31*dSdtot21131 + F22*F32*dSdtot22131 + F22*F33*dSdtot23131 + F23*F31*dSdtot31131 + F23*F32*dSdtot32131 + F23*F33*dSdtot33131;
        term(3,8) = F21*F31*dSdtot11121 + F21*F32*dSdtot12121 + F21*F33*dSdtot13121 + F22*F31*dSdtot21121 + F22*F32*dSdtot22121 + F22*F33*dSdtot23121 + F23*F31*dSdtot31121 + F23*F32*dSdtot32121 + F23*F33*dSdtot33121;
        term(3,9) = F21*F31*dSdtot11211 + F21*F32*dSdtot12211 + F21*F33*dSdtot13211 + F22*F31*dSdtot21211 + F22*F32*dSdtot22211 + F22*F33*dSdtot23211 + F23*F31*dSdtot31211 + F23*F32*dSdtot32211 + F23*F33*dSdtot33211;
        term(3,10) = F21*F31*dSdtot11222 + F21*F32*dSdtot12222 + F21*F33*dSdtot13222 + F22*F31*dSdtot21222 + F22*F32*dSdtot22222 + F22*F33*dSdtot23222 + F23*F31*dSdtot31222 + F23*F32*dSdtot32222 + F23*F33*dSdtot33222;
        term(3,11) = F21*F31*dSdtot11233 + F21*F32*dSdtot12233 + F21*F33*dSdtot13233 + F22*F31*dSdtot21233 + F22*F32*dSdtot22233 + F22*F33*dSdtot23233 + F23*F31*dSdtot31233 + F23*F32*dSdtot32233 + F23*F33*dSdtot33233;
        term(3,12) = F21*F31*dSdtot11223 + F21*F32*dSdtot12223 + F21*F33*dSdtot13223 + F22*F31*dSdtot21223 + F22*F32*dSdtot22223 + F22*F33*dSdtot23223 + F23*F31*dSdtot31223 + F23*F32*dSdtot32223 + F23*F33*dSdtot33223;
        term(3,13) = F21*F31*dSdtot11213 + F21*F32*dSdtot12213 + F21*F33*dSdtot13213 + F22*F31*dSdtot21213 + F22*F32*dSdtot22213 + F22*F33*dSdtot23213 + F23*F31*dSdtot31213 + F23*F32*dSdtot32213 + F23*F33*dSdtot33213;
        term(3,14) = F21*F31*dSdtot11212 + F21*F32*dSdtot12212 + F21*F33*dSdtot13212 + F22*F31*dSdtot21212 + F22*F32*dSdtot22212 + F22*F33*dSdtot23212 + F23*F31*dSdtot31212 + F23*F32*dSdtot32212 + F23*F33*dSdtot33212;
        term(3,15) = F21*F31*dSdtot11232 + F21*F32*dSdtot12232 + F21*F33*dSdtot13232 + F22*F31*dSdtot21232 + F22*F32*dSdtot22232 + F22*F33*dSdtot23232 + F23*F31*dSdtot31232 + F23*F32*dSdtot32232 + F23*F33*dSdtot33232;
        term(3,16) = F21*F31*dSdtot11231 + F21*F32*dSdtot12231 + F21*F33*dSdtot13231 + F22*F31*dSdtot21231 + F22*F32*dSdtot22231 + F22*F33*dSdtot23231 + F23*F31*dSdtot31231 + F23*F32*dSdtot32231 + F23*F33*dSdtot33231;
        term(3,17) = F21*F31*dSdtot11221 + F21*F32*dSdtot12221 + F21*F33*dSdtot13221 + F22*F31*dSdtot21221 + F22*F32*dSdtot22221 + F22*F33*dSdtot23221 + F23*F31*dSdtot31221 + F23*F32*dSdtot32221 + F23*F33*dSdtot33221;
        term(3,18) = F21*F31*dSdtot11311 + F21*F32*dSdtot12311 + F21*F33*dSdtot13311 + F22*F31*dSdtot21311 + F22*F32*dSdtot22311 + F22*F33*dSdtot23311 + F23*F31*dSdtot31311 + F23*F32*dSdtot32311 + F23*F33*dSdtot33311;
        term(3,19) = F21*F31*dSdtot11322 + F21*F32*dSdtot12322 + F21*F33*dSdtot13322 + F22*F31*dSdtot21322 + F22*F32*dSdtot22322 + F22*F33*dSdtot23322 + F23*F31*dSdtot31322 + F23*F32*dSdtot32322 + F23*F33*dSdtot33322;
        term(3,20) = F21*F31*dSdtot11333 + F21*F32*dSdtot12333 + F21*F33*dSdtot13333 + F22*F31*dSdtot21333 + F22*F32*dSdtot22333 + F22*F33*dSdtot23333 + F23*F31*dSdtot31333 + F23*F32*dSdtot32333 + F23*F33*dSdtot33333;
        term(3,21) = F21*F31*dSdtot11323 + F21*F32*dSdtot12323 + F21*F33*dSdtot13323 + F22*F31*dSdtot21323 + F22*F32*dSdtot22323 + F22*F33*dSdtot23323 + F23*F31*dSdtot31323 + F23*F32*dSdtot32323 + F23*F33*dSdtot33323;
        term(3,22) = F21*F31*dSdtot11313 + F21*F32*dSdtot12313 + F21*F33*dSdtot13313 + F22*F31*dSdtot21313 + F22*F32*dSdtot22313 + F22*F33*dSdtot23313 + F23*F31*dSdtot31313 + F23*F32*dSdtot32313 + F23*F33*dSdtot33313;
        term(3,23) = F21*F31*dSdtot11312 + F21*F32*dSdtot12312 + F21*F33*dSdtot13312 + F22*F31*dSdtot21312 + F22*F32*dSdtot22312 + F22*F33*dSdtot23312 + F23*F31*dSdtot31312 + F23*F32*dSdtot32312 + F23*F33*dSdtot33312;
        term(3,24) = F21*F31*dSdtot11332 + F21*F32*dSdtot12332 + F21*F33*dSdtot13332 + F22*F31*dSdtot21332 + F22*F32*dSdtot22332 + F22*F33*dSdtot23332 + F23*F31*dSdtot31332 + F23*F32*dSdtot32332 + F23*F33*dSdtot33332;
        term(3,25) = F21*F31*dSdtot11331 + F21*F32*dSdtot12331 + F21*F33*dSdtot13331 + F22*F31*dSdtot21331 + F22*F32*dSdtot22331 + F22*F33*dSdtot23331 + F23*F31*dSdtot31331 + F23*F32*dSdtot32331 + F23*F33*dSdtot33331;
        term(3,26) = F21*F31*dSdtot11321 + F21*F32*dSdtot12321 + F21*F33*dSdtot13321 + F22*F31*dSdtot21321 + F22*F32*dSdtot22321 + F22*F33*dSdtot23321 + F23*F31*dSdtot31321 + F23*F32*dSdtot32321 + F23*F33*dSdtot33321;
        term(4,0) = F11*F31*dSdtot11111 + F11*F32*dSdtot12111 + F11*F33*dSdtot13111 + F12*F31*dSdtot21111 + F12*F32*dSdtot22111 + F12*F33*dSdtot23111 + F13*F31*dSdtot31111 + F13*F32*dSdtot32111 + F13*F33*dSdtot33111;
        term(4,1) = F11*F31*dSdtot11122 + F11*F32*dSdtot12122 + F11*F33*dSdtot13122 + F12*F31*dSdtot21122 + F12*F32*dSdtot22122 + F12*F33*dSdtot23122 + F13*F31*dSdtot31122 + F13*F32*dSdtot32122 + F13*F33*dSdtot33122;
        term(4,2) = F11*F31*dSdtot11133 + F11*F32*dSdtot12133 + F11*F33*dSdtot13133 + F12*F31*dSdtot21133 + F12*F32*dSdtot22133 + F12*F33*dSdtot23133 + F13*F31*dSdtot31133 + F13*F32*dSdtot32133 + F13*F33*dSdtot33133;
        term(4,3) = F11*F31*dSdtot11123 + F11*F32*dSdtot12123 + F11*F33*dSdtot13123 + F12*F31*dSdtot21123 + F12*F32*dSdtot22123 + F12*F33*dSdtot23123 + F13*F31*dSdtot31123 + F13*F32*dSdtot32123 + F13*F33*dSdtot33123;
        term(4,4) = F11*F31*dSdtot11113 + F11*F32*dSdtot12113 + F11*F33*dSdtot13113 + F12*F31*dSdtot21113 + F12*F32*dSdtot22113 + F12*F33*dSdtot23113 + F13*F31*dSdtot31113 + F13*F32*dSdtot32113 + F13*F33*dSdtot33113;
        term(4,5) = F11*F31*dSdtot11112 + F11*F32*dSdtot12112 + F11*F33*dSdtot13112 + F12*F31*dSdtot21112 + F12*F32*dSdtot22112 + F12*F33*dSdtot23112 + F13*F31*dSdtot31112 + F13*F32*dSdtot32112 + F13*F33*dSdtot33112;
        term(4,6) = F11*F31*dSdtot11132 + F11*F32*dSdtot12132 + F11*F33*dSdtot13132 + F12*F31*dSdtot21132 + F12*F32*dSdtot22132 + F12*F33*dSdtot23132 + F13*F31*dSdtot31132 + F13*F32*dSdtot32132 + F13*F33*dSdtot33132;
        term(4,7) = F11*F31*dSdtot11131 + F11*F32*dSdtot12131 + F11*F33*dSdtot13131 + F12*F31*dSdtot21131 + F12*F32*dSdtot22131 + F12*F33*dSdtot23131 + F13*F31*dSdtot31131 + F13*F32*dSdtot32131 + F13*F33*dSdtot33131;
        term(4,8) = F11*F31*dSdtot11121 + F11*F32*dSdtot12121 + F11*F33*dSdtot13121 + F12*F31*dSdtot21121 + F12*F32*dSdtot22121 + F12*F33*dSdtot23121 + F13*F31*dSdtot31121 + F13*F32*dSdtot32121 + F13*F33*dSdtot33121;
        term(4,9) = F11*F31*dSdtot11211 + F11*F32*dSdtot12211 + F11*F33*dSdtot13211 + F12*F31*dSdtot21211 + F12*F32*dSdtot22211 + F12*F33*dSdtot23211 + F13*F31*dSdtot31211 + F13*F32*dSdtot32211 + F13*F33*dSdtot33211;
        term(4,10) = F11*F31*dSdtot11222 + F11*F32*dSdtot12222 + F11*F33*dSdtot13222 + F12*F31*dSdtot21222 + F12*F32*dSdtot22222 + F12*F33*dSdtot23222 + F13*F31*dSdtot31222 + F13*F32*dSdtot32222 + F13*F33*dSdtot33222;
        term(4,11) = F11*F31*dSdtot11233 + F11*F32*dSdtot12233 + F11*F33*dSdtot13233 + F12*F31*dSdtot21233 + F12*F32*dSdtot22233 + F12*F33*dSdtot23233 + F13*F31*dSdtot31233 + F13*F32*dSdtot32233 + F13*F33*dSdtot33233;
        term(4,12) = F11*F31*dSdtot11223 + F11*F32*dSdtot12223 + F11*F33*dSdtot13223 + F12*F31*dSdtot21223 + F12*F32*dSdtot22223 + F12*F33*dSdtot23223 + F13*F31*dSdtot31223 + F13*F32*dSdtot32223 + F13*F33*dSdtot33223;
        term(4,13) = F11*F31*dSdtot11213 + F11*F32*dSdtot12213 + F11*F33*dSdtot13213 + F12*F31*dSdtot21213 + F12*F32*dSdtot22213 + F12*F33*dSdtot23213 + F13*F31*dSdtot31213 + F13*F32*dSdtot32213 + F13*F33*dSdtot33213;
        term(4,14) = F11*F31*dSdtot11212 + F11*F32*dSdtot12212 + F11*F33*dSdtot13212 + F12*F31*dSdtot21212 + F12*F32*dSdtot22212 + F12*F33*dSdtot23212 + F13*F31*dSdtot31212 + F13*F32*dSdtot32212 + F13*F33*dSdtot33212;
        term(4,15) = F11*F31*dSdtot11232 + F11*F32*dSdtot12232 + F11*F33*dSdtot13232 + F12*F31*dSdtot21232 + F12*F32*dSdtot22232 + F12*F33*dSdtot23232 + F13*F31*dSdtot31232 + F13*F32*dSdtot32232 + F13*F33*dSdtot33232;
        term(4,16) = F11*F31*dSdtot11231 + F11*F32*dSdtot12231 + F11*F33*dSdtot13231 + F12*F31*dSdtot21231 + F12*F32*dSdtot22231 + F12*F33*dSdtot23231 + F13*F31*dSdtot31231 + F13*F32*dSdtot32231 + F13*F33*dSdtot33231;
        term(4,17) = F11*F31*dSdtot11221 + F11*F32*dSdtot12221 + F11*F33*dSdtot13221 + F12*F31*dSdtot21221 + F12*F32*dSdtot22221 + F12*F33*dSdtot23221 + F13*F31*dSdtot31221 + F13*F32*dSdtot32221 + F13*F33*dSdtot33221;
        term(4,18) = F11*F31*dSdtot11311 + F11*F32*dSdtot12311 + F11*F33*dSdtot13311 + F12*F31*dSdtot21311 + F12*F32*dSdtot22311 + F12*F33*dSdtot23311 + F13*F31*dSdtot31311 + F13*F32*dSdtot32311 + F13*F33*dSdtot33311;
        term(4,19) = F11*F31*dSdtot11322 + F11*F32*dSdtot12322 + F11*F33*dSdtot13322 + F12*F31*dSdtot21322 + F12*F32*dSdtot22322 + F12*F33*dSdtot23322 + F13*F31*dSdtot31322 + F13*F32*dSdtot32322 + F13*F33*dSdtot33322;
        term(4,20) = F11*F31*dSdtot11333 + F11*F32*dSdtot12333 + F11*F33*dSdtot13333 + F12*F31*dSdtot21333 + F12*F32*dSdtot22333 + F12*F33*dSdtot23333 + F13*F31*dSdtot31333 + F13*F32*dSdtot32333 + F13*F33*dSdtot33333;
        term(4,21) = F11*F31*dSdtot11323 + F11*F32*dSdtot12323 + F11*F33*dSdtot13323 + F12*F31*dSdtot21323 + F12*F32*dSdtot22323 + F12*F33*dSdtot23323 + F13*F31*dSdtot31323 + F13*F32*dSdtot32323 + F13*F33*dSdtot33323;
        term(4,22) = F11*F31*dSdtot11313 + F11*F32*dSdtot12313 + F11*F33*dSdtot13313 + F12*F31*dSdtot21313 + F12*F32*dSdtot22313 + F12*F33*dSdtot23313 + F13*F31*dSdtot31313 + F13*F32*dSdtot32313 + F13*F33*dSdtot33313;
        term(4,23) = F11*F31*dSdtot11312 + F11*F32*dSdtot12312 + F11*F33*dSdtot13312 + F12*F31*dSdtot21312 + F12*F32*dSdtot22312 + F12*F33*dSdtot23312 + F13*F31*dSdtot31312 + F13*F32*dSdtot32312 + F13*F33*dSdtot33312;
        term(4,24) = F11*F31*dSdtot11332 + F11*F32*dSdtot12332 + F11*F33*dSdtot13332 + F12*F31*dSdtot21332 + F12*F32*dSdtot22332 + F12*F33*dSdtot23332 + F13*F31*dSdtot31332 + F13*F32*dSdtot32332 + F13*F33*dSdtot33332;
        term(4,25) = F11*F31*dSdtot11331 + F11*F32*dSdtot12331 + F11*F33*dSdtot13331 + F12*F31*dSdtot21331 + F12*F32*dSdtot22331 + F12*F33*dSdtot23331 + F13*F31*dSdtot31331 + F13*F32*dSdtot32331 + F13*F33*dSdtot33331;
        term(5,26) = F11*F31*dSdtot11321 + F11*F32*dSdtot12321 + F11*F33*dSdtot13321 + F12*F31*dSdtot21321 + F12*F32*dSdtot22321 + F12*F33*dSdtot23321 + F13*F31*dSdtot31321 + F13*F32*dSdtot32321 + F13*F33*dSdtot33321;
        term(5,0) = F11*F21*dSdtot11111 + F11*F22*dSdtot12111 + F11*F23*dSdtot13111 + F12*F21*dSdtot21111 + F12*F22*dSdtot22111 + F12*F23*dSdtot23111 + F13*F21*dSdtot31111 + F13*F22*dSdtot32111 + F13*F23*dSdtot33111;
        term(5,1) = F11*F21*dSdtot11122 + F11*F22*dSdtot12122 + F11*F23*dSdtot13122 + F12*F21*dSdtot21122 + F12*F22*dSdtot22122 + F12*F23*dSdtot23122 + F13*F21*dSdtot31122 + F13*F22*dSdtot32122 + F13*F23*dSdtot33122;
        term(5,2) = F11*F21*dSdtot11133 + F11*F22*dSdtot12133 + F11*F23*dSdtot13133 + F12*F21*dSdtot21133 + F12*F22*dSdtot22133 + F12*F23*dSdtot23133 + F13*F21*dSdtot31133 + F13*F22*dSdtot32133 + F13*F23*dSdtot33133;
        term(5,3) = F11*F21*dSdtot11123 + F11*F22*dSdtot12123 + F11*F23*dSdtot13123 + F12*F21*dSdtot21123 + F12*F22*dSdtot22123 + F12*F23*dSdtot23123 + F13*F21*dSdtot31123 + F13*F22*dSdtot32123 + F13*F23*dSdtot33123;
        term(5,4) = F11*F21*dSdtot11113 + F11*F22*dSdtot12113 + F11*F23*dSdtot13113 + F12*F21*dSdtot21113 + F12*F22*dSdtot22113 + F12*F23*dSdtot23113 + F13*F21*dSdtot31113 + F13*F22*dSdtot32113 + F13*F23*dSdtot33113;
        term(5,5) = F11*F21*dSdtot11112 + F11*F22*dSdtot12112 + F11*F23*dSdtot13112 + F12*F21*dSdtot21112 + F12*F22*dSdtot22112 + F12*F23*dSdtot23112 + F13*F21*dSdtot31112 + F13*F22*dSdtot32112 + F13*F23*dSdtot33112;
        term(5,6) = F11*F21*dSdtot11132 + F11*F22*dSdtot12132 + F11*F23*dSdtot13132 + F12*F21*dSdtot21132 + F12*F22*dSdtot22132 + F12*F23*dSdtot23132 + F13*F21*dSdtot31132 + F13*F22*dSdtot32132 + F13*F23*dSdtot33132;
        term(5,7) = F11*F21*dSdtot11131 + F11*F22*dSdtot12131 + F11*F23*dSdtot13131 + F12*F21*dSdtot21131 + F12*F22*dSdtot22131 + F12*F23*dSdtot23131 + F13*F21*dSdtot31131 + F13*F22*dSdtot32131 + F13*F23*dSdtot33131;
        term(5,8) = F11*F21*dSdtot11121 + F11*F22*dSdtot12121 + F11*F23*dSdtot13121 + F12*F21*dSdtot21121 + F12*F22*dSdtot22121 + F12*F23*dSdtot23121 + F13*F21*dSdtot31121 + F13*F22*dSdtot32121 + F13*F23*dSdtot33121;
        term(5,9) = F11*F21*dSdtot11211 + F11*F22*dSdtot12211 + F11*F23*dSdtot13211 + F12*F21*dSdtot21211 + F12*F22*dSdtot22211 + F12*F23*dSdtot23211 + F13*F21*dSdtot31211 + F13*F22*dSdtot32211 + F13*F23*dSdtot33211;
        term(5,10) = F11*F21*dSdtot11222 + F11*F22*dSdtot12222 + F11*F23*dSdtot13222 + F12*F21*dSdtot21222 + F12*F22*dSdtot22222 + F12*F23*dSdtot23222 + F13*F21*dSdtot31222 + F13*F22*dSdtot32222 + F13*F23*dSdtot33222;
        term(5,11) = F11*F21*dSdtot11233 + F11*F22*dSdtot12233 + F11*F23*dSdtot13233 + F12*F21*dSdtot21233 + F12*F22*dSdtot22233 + F12*F23*dSdtot23233 + F13*F21*dSdtot31233 + F13*F22*dSdtot32233 + F13*F23*dSdtot33233;
        term(5,12) = F11*F21*dSdtot11223 + F11*F22*dSdtot12223 + F11*F23*dSdtot13223 + F12*F21*dSdtot21223 + F12*F22*dSdtot22223 + F12*F23*dSdtot23223 + F13*F21*dSdtot31223 + F13*F22*dSdtot32223 + F13*F23*dSdtot33223;
        term(5,13) = F11*F21*dSdtot11213 + F11*F22*dSdtot12213 + F11*F23*dSdtot13213 + F12*F21*dSdtot21213 + F12*F22*dSdtot22213 + F12*F23*dSdtot23213 + F13*F21*dSdtot31213 + F13*F22*dSdtot32213 + F13*F23*dSdtot33213;
        term(5,14) = F11*F21*dSdtot11212 + F11*F22*dSdtot12212 + F11*F23*dSdtot13212 + F12*F21*dSdtot21212 + F12*F22*dSdtot22212 + F12*F23*dSdtot23212 + F13*F21*dSdtot31212 + F13*F22*dSdtot32212 + F13*F23*dSdtot33212;
        term(5,15) = F11*F21*dSdtot11232 + F11*F22*dSdtot12232 + F11*F23*dSdtot13232 + F12*F21*dSdtot21232 + F12*F22*dSdtot22232 + F12*F23*dSdtot23232 + F13*F21*dSdtot31232 + F13*F22*dSdtot32232 + F13*F23*dSdtot33232;
        term(5,16) = F11*F21*dSdtot11231 + F11*F22*dSdtot12231 + F11*F23*dSdtot13231 + F12*F21*dSdtot21231 + F12*F22*dSdtot22231 + F12*F23*dSdtot23231 + F13*F21*dSdtot31231 + F13*F22*dSdtot32231 + F13*F23*dSdtot33231;
        term(5,17) = F11*F21*dSdtot11221 + F11*F22*dSdtot12221 + F11*F23*dSdtot13221 + F12*F21*dSdtot21221 + F12*F22*dSdtot22221 + F12*F23*dSdtot23221 + F13*F21*dSdtot31221 + F13*F22*dSdtot32221 + F13*F23*dSdtot33221;
        term(5,18) = F11*F21*dSdtot11311 + F11*F22*dSdtot12311 + F11*F23*dSdtot13311 + F12*F21*dSdtot21311 + F12*F22*dSdtot22311 + F12*F23*dSdtot23311 + F13*F21*dSdtot31311 + F13*F22*dSdtot32311 + F13*F23*dSdtot33311;
        term(5,19) = F11*F21*dSdtot11322 + F11*F22*dSdtot12322 + F11*F23*dSdtot13322 + F12*F21*dSdtot21322 + F12*F22*dSdtot22322 + F12*F23*dSdtot23322 + F13*F21*dSdtot31322 + F13*F22*dSdtot32322 + F13*F23*dSdtot33322;
        term(5,20) = F11*F21*dSdtot11333 + F11*F22*dSdtot12333 + F11*F23*dSdtot13333 + F12*F21*dSdtot21333 + F12*F22*dSdtot22333 + F12*F23*dSdtot23333 + F13*F21*dSdtot31333 + F13*F22*dSdtot32333 + F13*F23*dSdtot33333;
        term(5,21) = F11*F21*dSdtot11323 + F11*F22*dSdtot12323 + F11*F23*dSdtot13323 + F12*F21*dSdtot21323 + F12*F22*dSdtot22323 + F12*F23*dSdtot23323 + F13*F21*dSdtot31323 + F13*F22*dSdtot32323 + F13*F23*dSdtot33323;
        term(5,22) = F11*F21*dSdtot11313 + F11*F22*dSdtot12313 + F11*F23*dSdtot13313 + F12*F21*dSdtot21313 + F12*F22*dSdtot22313 + F12*F23*dSdtot23313 + F13*F21*dSdtot31313 + F13*F22*dSdtot32313 + F13*F23*dSdtot33313;
        term(5,23) = F11*F21*dSdtot11312 + F11*F22*dSdtot12312 + F11*F23*dSdtot13312 + F12*F21*dSdtot21312 + F12*F22*dSdtot22312 + F12*F23*dSdtot23312 + F13*F21*dSdtot31312 + F13*F22*dSdtot32312 + F13*F23*dSdtot33312;
        term(5,24) = F11*F21*dSdtot11332 + F11*F22*dSdtot12332 + F11*F23*dSdtot13332 + F12*F21*dSdtot21332 + F12*F22*dSdtot22332 + F12*F23*dSdtot23332 + F13*F21*dSdtot31332 + F13*F22*dSdtot32332 + F13*F23*dSdtot33332;
        term(5,25) = F11*F21*dSdtot11331 + F11*F22*dSdtot12331 + F11*F23*dSdtot13331 + F12*F21*dSdtot21331 + F12*F22*dSdtot22331 + F12*F23*dSdtot23331 + F13*F21*dSdtot31331 + F13*F22*dSdtot32331 + F13*F23*dSdtot33331;
        term(5,26) = F11*F21*dSdtot11321 + F11*F22*dSdtot12321 + F11*F23*dSdtot13321 + F12*F21*dSdtot21321 + F12*F22*dSdtot22321 + F12*F23*dSdtot23321 + F13*F21*dSdtot31321 + F13*F22*dSdtot32321 + F13*F23*dSdtot33321;
        term(6,0) = F21*F31*dSdtot11111 + F21*F32*dSdtot21111 + F21*F33*dSdtot31111 + F22*F31*dSdtot12111 + F22*F32*dSdtot22111 + F22*F33*dSdtot32111 + F23*F31*dSdtot13111 + F23*F32*dSdtot23111 + F23*F33*dSdtot33111;
        term(6,1) = F21*F31*dSdtot11122 + F21*F32*dSdtot21122 + F21*F33*dSdtot31122 + F22*F31*dSdtot12122 + F22*F32*dSdtot22122 + F22*F33*dSdtot32122 + F23*F31*dSdtot13122 + F23*F32*dSdtot23122 + F23*F33*dSdtot33122;
        term(6,2) = F21*F31*dSdtot11133 + F21*F32*dSdtot21133 + F21*F33*dSdtot31133 + F22*F31*dSdtot12133 + F22*F32*dSdtot22133 + F22*F33*dSdtot32133 + F23*F31*dSdtot13133 + F23*F32*dSdtot23133 + F23*F33*dSdtot33133;
        term(6,3) = F21*F31*dSdtot11123 + F21*F32*dSdtot21123 + F21*F33*dSdtot31123 + F22*F31*dSdtot12123 + F22*F32*dSdtot22123 + F22*F33*dSdtot32123 + F23*F31*dSdtot13123 + F23*F32*dSdtot23123 + F23*F33*dSdtot33123;
        term(6,4) = F21*F31*dSdtot11113 + F21*F32*dSdtot21113 + F21*F33*dSdtot31113 + F22*F31*dSdtot12113 + F22*F32*dSdtot22113 + F22*F33*dSdtot32113 + F23*F31*dSdtot13113 + F23*F32*dSdtot23113 + F23*F33*dSdtot33113;
        term(6,5) = F21*F31*dSdtot11112 + F21*F32*dSdtot21112 + F21*F33*dSdtot31112 + F22*F31*dSdtot12112 + F22*F32*dSdtot22112 + F22*F33*dSdtot32112 + F23*F31*dSdtot13112 + F23*F32*dSdtot23112 + F23*F33*dSdtot33112;
        term(6,6) = F21*F31*dSdtot11132 + F21*F32*dSdtot21132 + F21*F33*dSdtot31132 + F22*F31*dSdtot12132 + F22*F32*dSdtot22132 + F22*F33*dSdtot32132 + F23*F31*dSdtot13132 + F23*F32*dSdtot23132 + F23*F33*dSdtot33132;
        term(6,7) = F21*F31*dSdtot11131 + F21*F32*dSdtot21131 + F21*F33*dSdtot31131 + F22*F31*dSdtot12131 + F22*F32*dSdtot22131 + F22*F33*dSdtot32131 + F23*F31*dSdtot13131 + F23*F32*dSdtot23131 + F23*F33*dSdtot33131;
        term(6,8) = F21*F31*dSdtot11121 + F21*F32*dSdtot21121 + F21*F33*dSdtot31121 + F22*F31*dSdtot12121 + F22*F32*dSdtot22121 + F22*F33*dSdtot32121 + F23*F31*dSdtot13121 + F23*F32*dSdtot23121 + F23*F33*dSdtot33121;
        term(6,9) = F21*F31*dSdtot11211 + F21*F32*dSdtot21211 + F21*F33*dSdtot31211 + F22*F31*dSdtot12211 + F22*F32*dSdtot22211 + F22*F33*dSdtot32211 + F23*F31*dSdtot13211 + F23*F32*dSdtot23211 + F23*F33*dSdtot33211;
        term(6,10) = F21*F31*dSdtot11222 + F21*F32*dSdtot21222 + F21*F33*dSdtot31222 + F22*F31*dSdtot12222 + F22*F32*dSdtot22222 + F22*F33*dSdtot32222 + F23*F31*dSdtot13222 + F23*F32*dSdtot23222 + F23*F33*dSdtot33222;
        term(6,11) = F21*F31*dSdtot11233 + F21*F32*dSdtot21233 + F21*F33*dSdtot31233 + F22*F31*dSdtot12233 + F22*F32*dSdtot22233 + F22*F33*dSdtot32233 + F23*F31*dSdtot13233 + F23*F32*dSdtot23233 + F23*F33*dSdtot33233;
        term(6,12) = F21*F31*dSdtot11223 + F21*F32*dSdtot21223 + F21*F33*dSdtot31223 + F22*F31*dSdtot12223 + F22*F32*dSdtot22223 + F22*F33*dSdtot32223 + F23*F31*dSdtot13223 + F23*F32*dSdtot23223 + F23*F33*dSdtot33223;
        term(6,13) = F21*F31*dSdtot11213 + F21*F32*dSdtot21213 + F21*F33*dSdtot31213 + F22*F31*dSdtot12213 + F22*F32*dSdtot22213 + F22*F33*dSdtot32213 + F23*F31*dSdtot13213 + F23*F32*dSdtot23213 + F23*F33*dSdtot33213;
        term(6,14) = F21*F31*dSdtot11212 + F21*F32*dSdtot21212 + F21*F33*dSdtot31212 + F22*F31*dSdtot12212 + F22*F32*dSdtot22212 + F22*F33*dSdtot32212 + F23*F31*dSdtot13212 + F23*F32*dSdtot23212 + F23*F33*dSdtot33212;
        term(6,15) = F21*F31*dSdtot11232 + F21*F32*dSdtot21232 + F21*F33*dSdtot31232 + F22*F31*dSdtot12232 + F22*F32*dSdtot22232 + F22*F33*dSdtot32232 + F23*F31*dSdtot13232 + F23*F32*dSdtot23232 + F23*F33*dSdtot33232;
        term(6,16) = F21*F31*dSdtot11231 + F21*F32*dSdtot21231 + F21*F33*dSdtot31231 + F22*F31*dSdtot12231 + F22*F32*dSdtot22231 + F22*F33*dSdtot32231 + F23*F31*dSdtot13231 + F23*F32*dSdtot23231 + F23*F33*dSdtot33231;
        term(6,17) = F21*F31*dSdtot11221 + F21*F32*dSdtot21221 + F21*F33*dSdtot31221 + F22*F31*dSdtot12221 + F22*F32*dSdtot22221 + F22*F33*dSdtot32221 + F23*F31*dSdtot13221 + F23*F32*dSdtot23221 + F23*F33*dSdtot33221;
        term(6,18) = F21*F31*dSdtot11311 + F21*F32*dSdtot21311 + F21*F33*dSdtot31311 + F22*F31*dSdtot12311 + F22*F32*dSdtot22311 + F22*F33*dSdtot32311 + F23*F31*dSdtot13311 + F23*F32*dSdtot23311 + F23*F33*dSdtot33311;
        term(6,19) = F21*F31*dSdtot11322 + F21*F32*dSdtot21322 + F21*F33*dSdtot31322 + F22*F31*dSdtot12322 + F22*F32*dSdtot22322 + F22*F33*dSdtot32322 + F23*F31*dSdtot13322 + F23*F32*dSdtot23322 + F23*F33*dSdtot33322;
        term(6,20) = F21*F31*dSdtot11333 + F21*F32*dSdtot21333 + F21*F33*dSdtot31333 + F22*F31*dSdtot12333 + F22*F32*dSdtot22333 + F22*F33*dSdtot32333 + F23*F31*dSdtot13333 + F23*F32*dSdtot23333 + F23*F33*dSdtot33333;
        term(6,21) = F21*F31*dSdtot11323 + F21*F32*dSdtot21323 + F21*F33*dSdtot31323 + F22*F31*dSdtot12323 + F22*F32*dSdtot22323 + F22*F33*dSdtot32323 + F23*F31*dSdtot13323 + F23*F32*dSdtot23323 + F23*F33*dSdtot33323;
        term(6,22) = F21*F31*dSdtot11313 + F21*F32*dSdtot21313 + F21*F33*dSdtot31313 + F22*F31*dSdtot12313 + F22*F32*dSdtot22313 + F22*F33*dSdtot32313 + F23*F31*dSdtot13313 + F23*F32*dSdtot23313 + F23*F33*dSdtot33313;
        term(6,23) = F21*F31*dSdtot11312 + F21*F32*dSdtot21312 + F21*F33*dSdtot31312 + F22*F31*dSdtot12312 + F22*F32*dSdtot22312 + F22*F33*dSdtot32312 + F23*F31*dSdtot13312 + F23*F32*dSdtot23312 + F23*F33*dSdtot33312;
        term(6,24) = F21*F31*dSdtot11332 + F21*F32*dSdtot21332 + F21*F33*dSdtot31332 + F22*F31*dSdtot12332 + F22*F32*dSdtot22332 + F22*F33*dSdtot32332 + F23*F31*dSdtot13332 + F23*F32*dSdtot23332 + F23*F33*dSdtot33332;
        term(6,25) = F21*F31*dSdtot11331 + F21*F32*dSdtot21331 + F21*F33*dSdtot31331 + F22*F31*dSdtot12331 + F22*F32*dSdtot22331 + F22*F33*dSdtot32331 + F23*F31*dSdtot13331 + F23*F32*dSdtot23331 + F23*F33*dSdtot33331;
        term(6,26) = F21*F31*dSdtot11321 + F21*F32*dSdtot21321 + F21*F33*dSdtot31321 + F22*F31*dSdtot12321 + F22*F32*dSdtot22321 + F22*F33*dSdtot32321 + F23*F31*dSdtot13321 + F23*F32*dSdtot23321 + F23*F33*dSdtot33321;
        term(7,0) = F11*F31*dSdtot11111 + F11*F32*dSdtot21111 + F11*F33*dSdtot31111 + F12*F31*dSdtot12111 + F12*F32*dSdtot22111 + F12*F33*dSdtot32111 + F13*F31*dSdtot13111 + F13*F32*dSdtot23111 + F13*F33*dSdtot33111;
        term(7,1) = F11*F31*dSdtot11122 + F11*F32*dSdtot21122 + F11*F33*dSdtot31122 + F12*F31*dSdtot12122 + F12*F32*dSdtot22122 + F12*F33*dSdtot32122 + F13*F31*dSdtot13122 + F13*F32*dSdtot23122 + F13*F33*dSdtot33122;
        term(7,2) = F11*F31*dSdtot11133 + F11*F32*dSdtot21133 + F11*F33*dSdtot31133 + F12*F31*dSdtot12133 + F12*F32*dSdtot22133 + F12*F33*dSdtot32133 + F13*F31*dSdtot13133 + F13*F32*dSdtot23133 + F13*F33*dSdtot33133;
        term(7,3) = F11*F31*dSdtot11123 + F11*F32*dSdtot21123 + F11*F33*dSdtot31123 + F12*F31*dSdtot12123 + F12*F32*dSdtot22123 + F12*F33*dSdtot32123 + F13*F31*dSdtot13123 + F13*F32*dSdtot23123 + F13*F33*dSdtot33123;
        term(7,4) = F11*F31*dSdtot11113 + F11*F32*dSdtot21113 + F11*F33*dSdtot31113 + F12*F31*dSdtot12113 + F12*F32*dSdtot22113 + F12*F33*dSdtot32113 + F13*F31*dSdtot13113 + F13*F32*dSdtot23113 + F13*F33*dSdtot33113;
        term(7,5) = F11*F31*dSdtot11112 + F11*F32*dSdtot21112 + F11*F33*dSdtot31112 + F12*F31*dSdtot12112 + F12*F32*dSdtot22112 + F12*F33*dSdtot32112 + F13*F31*dSdtot13112 + F13*F32*dSdtot23112 + F13*F33*dSdtot33112;
        term(7,6) = F11*F31*dSdtot11132 + F11*F32*dSdtot21132 + F11*F33*dSdtot31132 + F12*F31*dSdtot12132 + F12*F32*dSdtot22132 + F12*F33*dSdtot32132 + F13*F31*dSdtot13132 + F13*F32*dSdtot23132 + F13*F33*dSdtot33132;
        term(7,7) = F11*F31*dSdtot11131 + F11*F32*dSdtot21131 + F11*F33*dSdtot31131 + F12*F31*dSdtot12131 + F12*F32*dSdtot22131 + F12*F33*dSdtot32131 + F13*F31*dSdtot13131 + F13*F32*dSdtot23131 + F13*F33*dSdtot33131;
        term(7,8) = F11*F31*dSdtot11121 + F11*F32*dSdtot21121 + F11*F33*dSdtot31121 + F12*F31*dSdtot12121 + F12*F32*dSdtot22121 + F12*F33*dSdtot32121 + F13*F31*dSdtot13121 + F13*F32*dSdtot23121 + F13*F33*dSdtot33121;
        term(7,9) = F11*F31*dSdtot11211 + F11*F32*dSdtot21211 + F11*F33*dSdtot31211 + F12*F31*dSdtot12211 + F12*F32*dSdtot22211 + F12*F33*dSdtot32211 + F13*F31*dSdtot13211 + F13*F32*dSdtot23211 + F13*F33*dSdtot33211;
        term(7,10) = F11*F31*dSdtot11222 + F11*F32*dSdtot21222 + F11*F33*dSdtot31222 + F12*F31*dSdtot12222 + F12*F32*dSdtot22222 + F12*F33*dSdtot32222 + F13*F31*dSdtot13222 + F13*F32*dSdtot23222 + F13*F33*dSdtot33222;
        term(7,11) = F11*F31*dSdtot11233 + F11*F32*dSdtot21233 + F11*F33*dSdtot31233 + F12*F31*dSdtot12233 + F12*F32*dSdtot22233 + F12*F33*dSdtot32233 + F13*F31*dSdtot13233 + F13*F32*dSdtot23233 + F13*F33*dSdtot33233;
        term(7,12) = F11*F31*dSdtot11223 + F11*F32*dSdtot21223 + F11*F33*dSdtot31223 + F12*F31*dSdtot12223 + F12*F32*dSdtot22223 + F12*F33*dSdtot32223 + F13*F31*dSdtot13223 + F13*F32*dSdtot23223 + F13*F33*dSdtot33223;
        term(7,13) = F11*F31*dSdtot11213 + F11*F32*dSdtot21213 + F11*F33*dSdtot31213 + F12*F31*dSdtot12213 + F12*F32*dSdtot22213 + F12*F33*dSdtot32213 + F13*F31*dSdtot13213 + F13*F32*dSdtot23213 + F13*F33*dSdtot33213;
        term(7,14) = F11*F31*dSdtot11212 + F11*F32*dSdtot21212 + F11*F33*dSdtot31212 + F12*F31*dSdtot12212 + F12*F32*dSdtot22212 + F12*F33*dSdtot32212 + F13*F31*dSdtot13212 + F13*F32*dSdtot23212 + F13*F33*dSdtot33212;
        term(7,15) = F11*F31*dSdtot11232 + F11*F32*dSdtot21232 + F11*F33*dSdtot31232 + F12*F31*dSdtot12232 + F12*F32*dSdtot22232 + F12*F33*dSdtot32232 + F13*F31*dSdtot13232 + F13*F32*dSdtot23232 + F13*F33*dSdtot33232;
        term(7,16) = F11*F31*dSdtot11231 + F11*F32*dSdtot21231 + F11*F33*dSdtot31231 + F12*F31*dSdtot12231 + F12*F32*dSdtot22231 + F12*F33*dSdtot32231 + F13*F31*dSdtot13231 + F13*F32*dSdtot23231 + F13*F33*dSdtot33231;
        term(7,17) = F11*F31*dSdtot11221 + F11*F32*dSdtot21221 + F11*F33*dSdtot31221 + F12*F31*dSdtot12221 + F12*F32*dSdtot22221 + F12*F33*dSdtot32221 + F13*F31*dSdtot13221 + F13*F32*dSdtot23221 + F13*F33*dSdtot33221;
        term(7,18) = F11*F31*dSdtot11311 + F11*F32*dSdtot21311 + F11*F33*dSdtot31311 + F12*F31*dSdtot12311 + F12*F32*dSdtot22311 + F12*F33*dSdtot32311 + F13*F31*dSdtot13311 + F13*F32*dSdtot23311 + F13*F33*dSdtot33311;
        term(7,19) = F11*F31*dSdtot11322 + F11*F32*dSdtot21322 + F11*F33*dSdtot31322 + F12*F31*dSdtot12322 + F12*F32*dSdtot22322 + F12*F33*dSdtot32322 + F13*F31*dSdtot13322 + F13*F32*dSdtot23322 + F13*F33*dSdtot33322;
        term(7,20) = F11*F31*dSdtot11333 + F11*F32*dSdtot21333 + F11*F33*dSdtot31333 + F12*F31*dSdtot12333 + F12*F32*dSdtot22333 + F12*F33*dSdtot32333 + F13*F31*dSdtot13333 + F13*F32*dSdtot23333 + F13*F33*dSdtot33333;
        term(7,21) = F11*F31*dSdtot11323 + F11*F32*dSdtot21323 + F11*F33*dSdtot31323 + F12*F31*dSdtot12323 + F12*F32*dSdtot22323 + F12*F33*dSdtot32323 + F13*F31*dSdtot13323 + F13*F32*dSdtot23323 + F13*F33*dSdtot33323;
        term(7,22) = F11*F31*dSdtot11313 + F11*F32*dSdtot21313 + F11*F33*dSdtot31313 + F12*F31*dSdtot12313 + F12*F32*dSdtot22313 + F12*F33*dSdtot32313 + F13*F31*dSdtot13313 + F13*F32*dSdtot23313 + F13*F33*dSdtot33313;
        term(7,23) = F11*F31*dSdtot11312 + F11*F32*dSdtot21312 + F11*F33*dSdtot31312 + F12*F31*dSdtot12312 + F12*F32*dSdtot22312 + F12*F33*dSdtot32312 + F13*F31*dSdtot13312 + F13*F32*dSdtot23312 + F13*F33*dSdtot33312;
        term(7,24) = F11*F31*dSdtot11332 + F11*F32*dSdtot21332 + F11*F33*dSdtot31332 + F12*F31*dSdtot12332 + F12*F32*dSdtot22332 + F12*F33*dSdtot32332 + F13*F31*dSdtot13332 + F13*F32*dSdtot23332 + F13*F33*dSdtot33332;
        term(7,25) = F11*F31*dSdtot11331 + F11*F32*dSdtot21331 + F11*F33*dSdtot31331 + F12*F31*dSdtot12331 + F12*F32*dSdtot22331 + F12*F33*dSdtot32331 + F13*F31*dSdtot13331 + F13*F32*dSdtot23331 + F13*F33*dSdtot33331;
        term(7,26) = F11*F31*dSdtot11321 + F11*F32*dSdtot21321 + F11*F33*dSdtot31321 + F12*F31*dSdtot12321 + F12*F32*dSdtot22321 + F12*F33*dSdtot32321 + F13*F31*dSdtot13321 + F13*F32*dSdtot23321 + F13*F33*dSdtot33321;
        term(8,0) = F11*F21*dSdtot11111 + F11*F22*dSdtot21111 + F11*F23*dSdtot31111 + F12*F21*dSdtot12111 + F12*F22*dSdtot22111 + F12*F23*dSdtot32111 + F13*F21*dSdtot13111 + F13*F22*dSdtot23111 + F13*F23*dSdtot33111;
        term(8,1) = F11*F21*dSdtot11122 + F11*F22*dSdtot21122 + F11*F23*dSdtot31122 + F12*F21*dSdtot12122 + F12*F22*dSdtot22122 + F12*F23*dSdtot32122 + F13*F21*dSdtot13122 + F13*F22*dSdtot23122 + F13*F23*dSdtot33122;
        term(8,2) = F11*F21*dSdtot11133 + F11*F22*dSdtot21133 + F11*F23*dSdtot31133 + F12*F21*dSdtot12133 + F12*F22*dSdtot22133 + F12*F23*dSdtot32133 + F13*F21*dSdtot13133 + F13*F22*dSdtot23133 + F13*F23*dSdtot33133;
        term(8,3) = F11*F21*dSdtot11123 + F11*F22*dSdtot21123 + F11*F23*dSdtot31123 + F12*F21*dSdtot12123 + F12*F22*dSdtot22123 + F12*F23*dSdtot32123 + F13*F21*dSdtot13123 + F13*F22*dSdtot23123 + F13*F23*dSdtot33123;
        term(8,4) = F11*F21*dSdtot11113 + F11*F22*dSdtot21113 + F11*F23*dSdtot31113 + F12*F21*dSdtot12113 + F12*F22*dSdtot22113 + F12*F23*dSdtot32113 + F13*F21*dSdtot13113 + F13*F22*dSdtot23113 + F13*F23*dSdtot33113;
        term(8,5) = F11*F21*dSdtot11112 + F11*F22*dSdtot21112 + F11*F23*dSdtot31112 + F12*F21*dSdtot12112 + F12*F22*dSdtot22112 + F12*F23*dSdtot32112 + F13*F21*dSdtot13112 + F13*F22*dSdtot23112 + F13*F23*dSdtot33112;
        term(8,6) = F11*F21*dSdtot11132 + F11*F22*dSdtot21132 + F11*F23*dSdtot31132 + F12*F21*dSdtot12132 + F12*F22*dSdtot22132 + F12*F23*dSdtot32132 + F13*F21*dSdtot13132 + F13*F22*dSdtot23132 + F13*F23*dSdtot33132;
        term(8,7) = F11*F21*dSdtot11131 + F11*F22*dSdtot21131 + F11*F23*dSdtot31131 + F12*F21*dSdtot12131 + F12*F22*dSdtot22131 + F12*F23*dSdtot32131 + F13*F21*dSdtot13131 + F13*F22*dSdtot23131 + F13*F23*dSdtot33131;
        term(8,8) = F11*F21*dSdtot11121 + F11*F22*dSdtot21121 + F11*F23*dSdtot31121 + F12*F21*dSdtot12121 + F12*F22*dSdtot22121 + F12*F23*dSdtot32121 + F13*F21*dSdtot13121 + F13*F22*dSdtot23121 + F13*F23*dSdtot33121;
        term(8,9) = F11*F21*dSdtot11211 + F11*F22*dSdtot21211 + F11*F23*dSdtot31211 + F12*F21*dSdtot12211 + F12*F22*dSdtot22211 + F12*F23*dSdtot32211 + F13*F21*dSdtot13211 + F13*F22*dSdtot23211 + F13*F23*dSdtot33211;
        term(8,10) = F11*F21*dSdtot11222 + F11*F22*dSdtot21222 + F11*F23*dSdtot31222 + F12*F21*dSdtot12222 + F12*F22*dSdtot22222 + F12*F23*dSdtot32222 + F13*F21*dSdtot13222 + F13*F22*dSdtot23222 + F13*F23*dSdtot33222;
        term(8,11) = F11*F21*dSdtot11233 + F11*F22*dSdtot21233 + F11*F23*dSdtot31233 + F12*F21*dSdtot12233 + F12*F22*dSdtot22233 + F12*F23*dSdtot32233 + F13*F21*dSdtot13233 + F13*F22*dSdtot23233 + F13*F23*dSdtot33233;
        term(8,12) = F11*F21*dSdtot11223 + F11*F22*dSdtot21223 + F11*F23*dSdtot31223 + F12*F21*dSdtot12223 + F12*F22*dSdtot22223 + F12*F23*dSdtot32223 + F13*F21*dSdtot13223 + F13*F22*dSdtot23223 + F13*F23*dSdtot33223;
        term(8,13) = F11*F21*dSdtot11213 + F11*F22*dSdtot21213 + F11*F23*dSdtot31213 + F12*F21*dSdtot12213 + F12*F22*dSdtot22213 + F12*F23*dSdtot32213 + F13*F21*dSdtot13213 + F13*F22*dSdtot23213 + F13*F23*dSdtot33213;
        term(8,14) = F11*F21*dSdtot11212 + F11*F22*dSdtot21212 + F11*F23*dSdtot31212 + F12*F21*dSdtot12212 + F12*F22*dSdtot22212 + F12*F23*dSdtot32212 + F13*F21*dSdtot13212 + F13*F22*dSdtot23212 + F13*F23*dSdtot33212;
        term(8,15) = F11*F21*dSdtot11232 + F11*F22*dSdtot21232 + F11*F23*dSdtot31232 + F12*F21*dSdtot12232 + F12*F22*dSdtot22232 + F12*F23*dSdtot32232 + F13*F21*dSdtot13232 + F13*F22*dSdtot23232 + F13*F23*dSdtot33232;
        term(8,16) = F11*F21*dSdtot11231 + F11*F22*dSdtot21231 + F11*F23*dSdtot31231 + F12*F21*dSdtot12231 + F12*F22*dSdtot22231 + F12*F23*dSdtot32231 + F13*F21*dSdtot13231 + F13*F22*dSdtot23231 + F13*F23*dSdtot33231;
        term(8,17) = F11*F21*dSdtot11221 + F11*F22*dSdtot21221 + F11*F23*dSdtot31221 + F12*F21*dSdtot12221 + F12*F22*dSdtot22221 + F12*F23*dSdtot32221 + F13*F21*dSdtot13221 + F13*F22*dSdtot23221 + F13*F23*dSdtot33221;
        term(8,18) = F11*F21*dSdtot11311 + F11*F22*dSdtot21311 + F11*F23*dSdtot31311 + F12*F21*dSdtot12311 + F12*F22*dSdtot22311 + F12*F23*dSdtot32311 + F13*F21*dSdtot13311 + F13*F22*dSdtot23311 + F13*F23*dSdtot33311;
        term(8,19) = F11*F21*dSdtot11322 + F11*F22*dSdtot21322 + F11*F23*dSdtot31322 + F12*F21*dSdtot12322 + F12*F22*dSdtot22322 + F12*F23*dSdtot32322 + F13*F21*dSdtot13322 + F13*F22*dSdtot23322 + F13*F23*dSdtot33322;
        term(8,20) = F11*F21*dSdtot11333 + F11*F22*dSdtot21333 + F11*F23*dSdtot31333 + F12*F21*dSdtot12333 + F12*F22*dSdtot22333 + F12*F23*dSdtot32333 + F13*F21*dSdtot13333 + F13*F22*dSdtot23333 + F13*F23*dSdtot33333;
        term(8,21) = F11*F21*dSdtot11323 + F11*F22*dSdtot21323 + F11*F23*dSdtot31323 + F12*F21*dSdtot12323 + F12*F22*dSdtot22323 + F12*F23*dSdtot32323 + F13*F21*dSdtot13323 + F13*F22*dSdtot23323 + F13*F23*dSdtot33323;
        term(8,22) = F11*F21*dSdtot11313 + F11*F22*dSdtot21313 + F11*F23*dSdtot31313 + F12*F21*dSdtot12313 + F12*F22*dSdtot22313 + F12*F23*dSdtot32313 + F13*F21*dSdtot13313 + F13*F22*dSdtot23313 + F13*F23*dSdtot33313;
        term(8,23) = F11*F21*dSdtot11312 + F11*F22*dSdtot21312 + F11*F23*dSdtot31312 + F12*F21*dSdtot12312 + F12*F22*dSdtot22312 + F12*F23*dSdtot32312 + F13*F21*dSdtot13312 + F13*F22*dSdtot23312 + F13*F23*dSdtot33312;
        term(8,24) = F11*F21*dSdtot11332 + F11*F22*dSdtot21332 + F11*F23*dSdtot31332 + F12*F21*dSdtot12332 + F12*F22*dSdtot22332 + F12*F23*dSdtot32332 + F13*F21*dSdtot13332 + F13*F22*dSdtot23332 + F13*F23*dSdtot33332;
        term(8,25) = F11*F21*dSdtot11331 + F11*F22*dSdtot21331 + F11*F23*dSdtot31331 + F12*F21*dSdtot12331 + F12*F22*dSdtot22331 + F12*F23*dSdtot32331 + F13*F21*dSdtot13331 + F13*F22*dSdtot23331 + F13*F23*dSdtot33331;
        term(8,26) = F11*F21*dSdtot11321 + F11*F22*dSdtot21321 + F11*F23*dSdtot31321 + F12*F21*dSdtot12321 + F12*F22*dSdtot22321 + F12*F23*dSdtot32321 + F13*F21*dSdtot13321 + F13*F22*dSdtot23321 + F13*F23*dSdtot33321;
        
        term /= J;
        return;
    }

    void compute_dcauchydF(const Matrix_3x3 &F, const Vector_9 &cauchy, const Vector_9 &PK2, const Matrix_9x9 &dPK2dF, Matrix_9x9 &dcauchydF){
        /*!===========================
        |    compute_dcauchydF    |
        ===========================

        Compute the derivative of the cauchy stress with respect to the deformation gradient.

        */

        //Initialize temporary matrices
        Matrix_3x3 T1;

        //Compute the jacobian of deformation
        double J = F.determinant();

        //Convert the deformation gradient to voigt notation
        Vector_9 FinvT_voigt;
        deformation_measures::voigt_3x3_tensor(F.inverse().transpose(),FinvT_voigt);

        Matrix_3x3 PK2_mat;
        deformation_measures::undo_voigt_3x3_tensor(PK2,PK2_mat);

        //Populate term 1 of the jacobian
        Matrix_9x9 term1;
        for (int i=0; i<9; i++){
            for (int j=0; j<9; j++){
                term1(i,j) = -cauchy(i)*FinvT_voigt(j);
            }
        }

        //Populate term2 of the jacobian
        Matrix_9x9 term2 = Matrix_9x9::Zero();
        
        T1 = PK2_mat*F.transpose()/J;
        term2(0,0) = T1(0,0);
        term2(0,4) = T1(2,0);
        term2(0,5) = T1(1,0);
        term2(1,1) = T1(1,1);
        term2(1,3) = T1(2,1);
        term2(1,8) = T1(0,1);
        term2(2,2) = T1(2,2);
        term2(2,6) = T1(1,2);
        term2(2,7) = T1(0,2);
        term2(3,1) = T1(1,2);
        term2(3,3) = T1(2,2);
        term2(3,8) = T1(0,2);
        term2(4,0) = T1(0,2);
        term2(4,4) = T1(2,2);
        term2(4,5) = T1(1,2);
        term2(5,0) = T1(0,1);
        term2(5,4) = T1(2,1);
        term2(5,5) = T1(1,1);
        term2(6,2) = T1(2,1);
        term2(6,6) = T1(1,1);
        term2(6,7) = T1(0,1);
        term2(7,2) = T1(2,0);
        term2(7,6) = T1(1,0);
        term2(7,7) = T1(0,0);
        term2(8,1) = T1(1,0);
        term2(8,3) = T1(2,0);
        term2(8,8) = T1(0,0);

        //Populate term3 of the jacobian
        Matrix_9x9 term3 = Matrix_9x9::Zero();

        map_fot_stress_jacobian(F, dPK2dF, J, term3);

        //Populate term4 of the jacobian
        Matrix_9x9 term4 = Matrix_9x9::Zero();

        T1 = F*PK2_mat/J;
        term4(0,0) = T1(0,0);
        term4(0,4) = T1(0,2);
        term4(0,5) = T1(0,1);
        term4(1,1) = T1(1,1);
        term4(1,3) = T1(1,2);
        term4(1,8) = T1(1,0);
        term4(2,2) = T1(2,2);
        term4(2,6) = T1(2,1);
        term4(2,7) = T1(2,0);
        term4(3,2) = T1(1,2);
        term4(3,6) = T1(1,1);
        term4(3,7) = T1(1,0);
        term4(4,2) = T1(0,2);
        term4(4,6) = T1(0,1);
        term4(4,7) = T1(0,0);
        term4(5,1) = T1(0,1);
        term4(5,3) = T1(0,2);
        term4(5,8) = T1(0,0);
        term4(6,1) = T1(2,1);
        term4(6,3) = T1(2,2);
        term4(6,8) = T1(2,0);
        term4(7,0) = T1(2,0);
        term4(7,4) = T1(2,2);
        term4(7,5) = T1(2,1);
        term4(8,0) = T1(1,0);
        term4(8,4) = T1(1,2);
        term4(8,5) = T1(1,1);

        dcauchydF = term1 + term2 + term3 + term4;

        return;
    }

    void compute_dcauchydchi(const Matrix_3x3 &F, const Matrix_9x9 &dPK2dchi, Matrix_9x9 &dcauchydchi){
        /*!==========================
        |    compute_dcauchydchi    |
        =============================

        Compute the derivative of the cauchy stress w.r.t. chi.

        */

         double J = F.determinant();
         map_fot_stress_jacobian(F, dPK2dchi, J, dcauchydchi);

        return;
    }

    void compute_dcauchydgrad_chi(const Matrix_3x3 &F, const Matrix_9x27 &dPK2dgrad_chi, Matrix_9x27 &dcauchydgrad_chi){
        /*!==================================
        |    compute_dcauchydgrad_chi    |
        ==================================

        Compute the derivative of the cauchy stress w.r.t. the gradient of chi.

        */

        double J = F.determinant();
        map_fifthot_stress_jacobian(F, dPK2dgrad_chi, J, dcauchydgrad_chi);

        return;
    }


    void compute_dsdF(const Matrix_3x3 &F, const Vector_9 &s, const Vector_9 &SIGMA, const Matrix_9x9 &dSIGMAdF, Matrix_9x9 &dsdF){
        /*!======================
        |    compute_dsdF    |
        ======================

        Compute the derivative of the symmetric stress with respect to the deformation gradient.

        */

        compute_dcauchydF(F,s,SIGMA,dSIGMAdF,dsdF); //Note: The derivatives take (essentially) the same form as the cauchy stress.

        return;
    }

    void compute_dsdchi(const Matrix_3x3 &F, const Matrix_9x9 &dSIGMAdchi, Matrix_9x9 &dsdchi){
        /*!========================
        |    compute_dsdchi    |
        ========================

        Compute the derivative of the symmetric stress w.r.t. chi.

        */

        compute_dcauchydchi(F,dSIGMAdchi,dsdchi);

        return;
    }

    void compute_dsdgrad_chi(const Matrix_3x3 &F, const Matrix_9x27 &dSIGMAdgrad_chi, Matrix_9x27 &dsdgrad_chi){
        /*!=============================
        |    compute_dsdgrad_chi    |
        =============================

        Compute the derivative of the symmetric stress w.r.t. the gradient of chi.

        */

        compute_dcauchydgrad_chi(F,dsdgrad_chi,dsdgrad_chi);

        return;
    }
}
