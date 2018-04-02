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
  
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <deformation_measures.h>
#include <micromorphic_linear_elasticity_voigt.h>

namespace micro_material{

    void get_stress(double (&params)[18]){
    
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

        //Compute the stress measures

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
        tripletList.push_back(T(3,3,2*mu));
        tripletList.push_back(T(3,6,2*mu));
        tripletList.push_back(T(4,4,2*mu));
        tripletList.push_back(T(4,7,2*mu));
        tripletList.push_back(T(5,5,2*mu));
        tripletList.push_back(T(5,8,2*mu));
        tripletList.push_back(T(6,3,2*mu));
        tripletList.push_back(T(6,6,2*mu));
        tripletList.push_back(T(7,4,2*mu));
        tripletList.push_back(T(7,7,2*mu));
        tripletList.push_back(T(8,5,2*mu));
        tripletList.push_back(T(8,8,2*mu)); 
    
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
        tripletList.push_back(T(3,3,2*sigma));
        tripletList.push_back(T(3,6,2*sigma));
        tripletList.push_back(T(4,4,2*sigma));
        tripletList.push_back(T(4,7,2*sigma));
        tripletList.push_back(T(5,5,2*sigma));
        tripletList.push_back(T(5,8,2*sigma));
        tripletList.push_back(T(6,3,2*sigma));
        tripletList.push_back(T(6,6,2*sigma));
        tripletList.push_back(T(7,4,2*sigma));
        tripletList.push_back(T(7,7,2*sigma));
        tripletList.push_back(T(8,5,2*sigma));
        tripletList.push_back(T(8,8,2*sigma));        
        
        D.setFromTriplets(tripletList.begin(), tripletList.end());
        return;
    }
    
//    void compute_PK2_stress(const Vector_9 &E_voigt,    const Vector_9 &E_micro_voigt, const Vector_27 &Gamma_voigt,
//                            const Matrix_3x3 &Cinv,     const Matrix_3x3 &Psi,         const Matrix_3x9 &Gamma,
//                            const SpMat &A, const SpMat &B,    const SpMat &C,
//                            const SpMat &D, Vector_9 &PK2){
//        /*!============================
//           |    compute_PK2_stress    |
//           ============================
//           
//           Compute the second piola kirchoff stress.
//           
//        */
//        
//        PK2 = A*E_voigt;        //Compute the first terms
//        PK2 += B*E_micro_voigt;
//        
//        //Compute the middle terms
//        Matrix_3x3 Temp1;
//        deformation_measures::undo_voigt_3x3_tensor(B*E_micro_voigt+D*E_voigt,Temp1);
//        PK2 += Temp1*(Cinv*Psi).transpose();
//        
//        //Compute the end terms
//        Matrix_3x9 Temp2;
//        deformation_measures::undo_voigt_3x9_tensor(C*Gamma_voigt,Temp2);
//        PK2 += Temp2*(Cinv*Gamma).transpose();
//        
//        return;
//    }
//    
//    void compute_symmetric_stress(const Vector_9 &E_voigt,    const Vector_9 &E_micro_voigt, const Vector_27 &Gamma_voigt,
//                            const Matrix_3x3 &Cinv,     const Matrix_3x3 &Psi,         const Matrix_3x9 &Gamma,
//                            const SpMat &A, const SpMat &B,    const SpMat &C,
//                            const SpMat &D, Vector_9 &SIGMA){
//        /*!============================
//           |    compute_symmetric_stress    |
//           ============================
//           
//           Compute the symmetric stress.
//           
//        */
//        
//        SIGMA = A*E_voigt;        //Compute the first terms
//        SIGMA += B*E_micro_voigt;
//        
//        //Compute the middle terms
//        Matrix_3x3 Temp1;
//        deformation_measures::undo_voigt_3x3_tensor(B*E_micro_voigt+D*E_voigt,Temp1);
//        Matrix_3x3 symmetric_part = Temp1*(Cinv*Psi).transpose();
//        
//        //Compute the end terms
//        Matrix_3x9 Temp2;
//        deformation_measures::undo_voigt_3x9_tensor(C*Gamma_voigt,Temp2);
//        symmetric_part += Temp2*(Cinv*Gamma).transpose();
//        
//        SIGMA += (symmetric_part + symmetric_part.transpose());
//        
//        return;
//    }
//    
//    void compute_higher_order_stress(const Vector_27 &Gamma_voigt, const SpMat &C, Vector_27 &M){
//        /*!=====================================
//           |    compute_higher_order_stress    |
//           =====================================
//          
//           Compute the higher order stress.
//          
//        */
//        
//        M = C*Gamma_voigt; //Compute the stress (requires positive permutation
//        Matrix_3x9 swap;                                     //Declare swap space
//        deformation_measures::undo_voigt_3x9_tensor(M,swap); //Populate the swap space
//        deformation_measures::perform_positive_cyclic_permutation(swap); //Perform the permutation
//        deformation_measures::voigt_3x9_tensor(swap,M); //Put the tensor back in voigt notation
//    }
}
