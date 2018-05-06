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

    void LinearElasticity::evaluate_model(const std::vector<double> &time,        const std::vector<double> (&fparams),
                                          const double (&grad_u)[3][3],           const double (&phi)[9],
                                          const double (&grad_phi)[9][3],         std::vector<double> &SDVS,
                                          const std::vector<double> &ADD_DOF,     const std::vector<std::vector<double>> &ADD_grad_DOF,
                                          Vector_9 &cauchy, Vector_9 &s, Vector_27 &m, std::vector<Eigen::VectorXd> &ADD_TERMS){
        /*!
        =======================
        |    evaluate_model   |
        =======================
        
        Evaluate the constitutive model 
        from the general incoming values.
        
        Only returns stresses and additional 
        terms.
        
        */
        
        //Extract the time
        double t  = time[0];
        double dt = time[1];

        //Extract the parameters        
        double params[18];
        if(fparams.size() == 18){
            for(int i=0; i<18; i++){
                params[i] = fparams[i];
            }
        }
        else{std::cout << "Error: Material parameters incorrectly specified\n";}
        
        //Compute the required deformation measures
        Matrix_3x3 F;
        Matrix_3x3 chi;
        Matrix_3x9 grad_chi;
        get_deformation_measures(grad_u, phi, grad_phi, F, chi, grad_chi);
        
        //Compute the stresses
        Vector_9  PK2;
        Vector_9  SIGMA;
        Vector_27 M;
        
        get_stress(t, dt, params, F, chi, grad_chi, SDVS, PK2, SIGMA, M);
        deformation_measures::map_stresses_to_current_configuration(F, chi, PK2, SIGMA, M, cauchy, s, m);
        
        return;
    }
                                
    void LinearElasticity::evaluate_model(const std::vector<double> &time,        const std::vector<double> (&fparams),
                                          const double (&grad_u)[3][3],           const double (&phi)[9],
                                          const double (&grad_phi)[9][3],         std::vector<double> &SDVS,
                                          const std::vector<double> &ADD_DOF,     const std::vector<std::vector<double>> &ADD_grad_DOF,
                                          Vector_9    &cauchy,    Vector_9    &s,           Vector_27    &m,
                                          Matrix_9x9  &DcauchyDgrad_u, Matrix_9x9  &DcauchyDphi, Matrix_9x27  &DcauchyDgrad_phi,
                                          Matrix_9x9  &DsDgrad_u,      Matrix_9x9  &DsDphi,      Matrix_9x27  &DsDgrad_phi,
                                          Matrix_27x9 &DmDgrad_u,      Matrix_27x9 &DmDphi,      Matrix_27x27 &DmDgrad_phi,
                                          std::vector<Eigen::VectorXd> &ADD_TERMS,               std::vector<Eigen::MatrixXd> &ADD_JACOBIANS){
        /*!
        ========================
        |    evaluate_model    |
        ========================
        
        Evaluate the constitutive model 
        from the general incoming values.
        
        Returns stresses, additional 
        terms, and their jacobians.
        
        */

        //Extract the time
        double t  = time[0];
        double dt = time[1];

        //Extract the parameters        
        double params[18];
        if(fparams.size() == 18){
            for(int i=0; i<18; i++){
                params[i] = fparams[i];
            }
        }
        else{std::cout << "Error: Material parameters incorrectly specified\n";assert(-21==-20);}

        //Compute the required deformation measures
        Matrix_3x3 F;
        Matrix_3x3 chi;
        Matrix_3x9 grad_chi;
        get_deformation_measures(grad_u, phi, grad_phi, F, chi, grad_chi);
        
        //Compute the stresses and jacobians
        Vector_9  PK2;
        Vector_9  SIGMA;
        Vector_27 M;
        
        Matrix_9x9   dPK2dF;
        Matrix_9x9   dPK2dchi;
        Matrix_9x27  dPK2dgrad_chi;
        Matrix_9x9   dSIGMAdF;
        Matrix_9x9   dSIGMAdchi;
        Matrix_9x27  dSIGMAdgrad_chi;
        Matrix_27x9  dMdF;
        Matrix_27x9  dMdchi;
        Matrix_27x27 dMdgrad_chi;
        
        //Note: D(x)Dphi = d(x)dchi
        Matrix_9x9   dcauchydF;
        Matrix_9x27  dcauchydgrad_chi;
        Matrix_9x9   dsdF;
        Matrix_9x27  dsdgrad_chi;
        Matrix_27x9  dmdF;
        Matrix_27x27 dmdgrad_chi;

        //assert(13==14);

        get_stress(t, dt, params, F, chi, grad_chi, SDVS, PK2, SIGMA, M,
                   dPK2dF,   dPK2dchi,   dPK2dgrad_chi,
                   dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi,
                   dMdF,     dMdchi,     dMdgrad_chi);
                   
        //assert(14==15);
        deformation_measures::map_stresses_to_current_configuration(F, chi, PK2, SIGMA, M, cauchy, s, m);

        //assert(15==16);
        
        deformation_measures::map_jacobians_to_current_configuration(F,      chi,      PK2,           SIGMA,     M,           cauchy, s,   m,
                                                                     dPK2dF, dPK2dchi, dPK2dgrad_chi, dSIGMAdF,  dSIGMAdchi,  dSIGMAdgrad_chi,
                                                                     dMdF,   dMdchi,   dMdgrad_chi,   dcauchydF, DcauchyDphi, dcauchydgrad_chi,
                                                                     dsdF,   DsDphi,   dsdgrad_chi,   dmdF,      DmDphi,      dmdgrad_chi);
        //assert(16==17);                                                                    
        Matrix_3x9 _grad_phi;
        Vector_27  _grad_phi_v;
        Matrix_3x3 eye = Matrix_3x3::Identity();
        deformation_measures::assemble_grad_chi(grad_phi, eye, _grad_phi); //Put grad_phi into an Eigen Matrix
        deformation_measures::voigt_3x9_tensor(_grad_phi,_grad_phi_v);     //Put grad_phi into voigt notation
        deformation_measures::compute_total_derivatives(F, _grad_phi_v,
                                                        dcauchydF,      dcauchydgrad_chi, dsdF,      dsdgrad_chi, dmdF,      dmdgrad_chi,
                                                        DcauchyDgrad_u, DcauchyDgrad_phi, DsDgrad_u, DsDgrad_phi, DmDgrad_u, DmDgrad_phi);
        //assert(17==18);
        return;
    }
    
    void LinearElasticity::get_deformation_measures(const double (&grad_u)[3][3], const double (&phi)[9], const double (&grad_phi)[9][3],
                                                    Matrix_3x3 &F,                Matrix_3x3 &chi,        Matrix_3x9 &grad_chi){
        /*!
        ==================================
        |    get_deformation_measures    |
        ==================================
        
        Compute the deformation measures from the degrees of freedom and their 
        gradients.
        
        */
        
        deformation_measures::get_deformation_gradient(grad_u, F);

        deformation_measures::assemble_chi(phi, chi);

        deformation_measures::assemble_grad_chi(grad_phi, F, grad_chi);
        
        return;
    }

    void get_stress(const double &t,     const double &dt,      const double (&params)[18],
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
        //SpMat A( 9, 9);
        //SpMat B( 9, 9);
        //SpMat C(27,27);
        //SpMat D( 9, 9);

        Matrix_9x9   A;
        Matrix_9x9   B;
        Matrix_27x27 C;
        Matrix_9x9   D;

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
    
    void get_stress(const double &t,     const double &dt,      const double (&params)[18],
                    const Matrix_3x3 &F, const Matrix_3x3 &chi, const Matrix_3x9 &grad_chi,
                    std::vector<double> &SDVS,    Vector_9 &PK2, Vector_9 &SIGMA, Vector_27 &M,
                    Matrix_9x9  &dPK2dF,   Matrix_9x9  &dPK2dchi,   Matrix_9x27  &dPK2dgrad_chi,
                    Matrix_9x9  &dSIGMAdF, Matrix_9x9  &dSIGMAdchi, Matrix_9x27  &dSIGMAdgrad_chi,
                    Matrix_27x9 &dMdF,     Matrix_27x9 &dMdchi,     Matrix_27x27 &dMdgrad_chi){
        /*!=================
        |    get_stress    |
        ====================
        
        Computes the stress measures and their jacobians.
        
        */
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
        //SpMat A( 9, 9);
        //SpMat B( 9, 9);
        //SpMat C(27,27);
        //SpMat D( 9, 9);

        Matrix_9x9   A;
        Matrix_9x9   B;
        Matrix_27x27 C;
        Matrix_9x9   D;

        //assert(101==102);
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
        
        //Compute the jacobians w.r.t. the derived deformation measures
        Matrix_9x9  dPK2dRCG;
        Matrix_9x9  dPK2dPsi;
        Matrix_9x27 dPK2dGamma;
        
        Matrix_9x9  dSIGMAdRCG;
        Matrix_9x9  dSIGMAdPsi;
        Matrix_9x27 dSIGMAdGamma;
        
        Matrix_27x27 dMdGamma; //Note: other gradients are zero for M in this form.
        
        //Terms to speed up computation.
        Matrix_9x9  dPK2dRCGterms[2];
        Matrix_9x9  dPK2dPsiterms[2];
        Matrix_9x27 dPK2dGammaterms[2];
        
        //Gradients of the PK2 stress
        compute_dPK2dRCG(RCG, RCGinv, Gamma, Gamma_voigt, E,             E_micro, E_voigt, E_micro_voigt,
                         A,   B,      C,     D,           dPK2dRCGterms, dPK2dRCG);
                        
        compute_dPK2dPsi(RCGinv, E_micro, E_voigt,       E_micro_voigt,
                          B,      D,       dPK2dPsiterms, dPK2dPsi);
                          
        compute_dPK2dGamma(RCGinv, Gamma,           Gamma_voigt,
                           C,      dPK2dGammaterms, dPK2dGamma);
        
        //Gradients of the symmetric stress (reference configuration)
        compute_dSIGMAdRCG(dPK2dRCGterms, dSIGMAdRCG);
    
        compute_dSIGMAdPsi(dPK2dPsiterms, dSIGMAdPsi);
    
        compute_dSIGMAdGamma(dPK2dGammaterms, dSIGMAdGamma);
        
        compute_dMdGamma(Matrix_27x27(C), dMdGamma);
        
        //Gradients of the derived measures
        //Note: Replaced sparse matricies with dense matrices
        //      This is not the most efficient but it seems required
        //      for use in MOOSE.
        Matrix_9x9   dRCGdF;
        
        Matrix_9x9   dPsidF;
        Matrix_9x9   dPsidchi;
        
        Matrix_27x9  dGammadF;
        Matrix_27x27 dGammadgrad_chi;

        deformation_measures::compute_dRCGdF(F,dRCGdF);

        deformation_measures::compute_dPsidF(chi,dPsidF);
        deformation_measures::compute_dPsidchi(F,dPsidchi);
        
        Vector_27 grad_chi_voigt;
        deformation_measures::voigt_3x9_tensor(grad_chi,grad_chi_voigt);
        
        deformation_measures::compute_dGammadF(grad_chi_voigt,dGammadF);
        deformation_measures::compute_dGammadgrad_chi(F, dGammadgrad_chi);

        //Compute the jacobians of the stresses w.r.t. the fundamental deformation measures.
        dPK2dF   = dPK2dRCG*dRCGdF   + dPK2dPsi*dPsidF   + dPK2dGamma*dGammadF;
        dSIGMAdF = dSIGMAdRCG*dRCGdF + dSIGMAdPsi*dPsidF + dSIGMAdGamma*dGammadF;
        dMdF     = dMdGamma*dGammadF; //Note: all other derivatives are zero.
        
        dPK2dchi   = dPK2dPsi*dPsidchi;
        dSIGMAdchi = dSIGMAdPsi*dPsidchi;
        dMdchi     = Matrix_27x9::Zero(); //Note: M is independent of the magnitude of chi (not so for grad_chi)
        
        dPK2dgrad_chi   = dPK2dGamma*dGammadgrad_chi;
        dSIGMAdgrad_chi = dSIGMAdGamma*dGammadgrad_chi;
        dMdgrad_chi     = dMdGamma*dGammadgrad_chi;
        
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

    void compute_A_voigt(const double &lambda,const double &mu, Matrix_9x9 &A) {
        /*!=========================
           |    compute_A_voigt    |
           =========================
           
           Compute the A stiffness matrix in voigt notation.
        */

        A = Matrix_9x9::Zero();

        A(0,0) = lambda + 2*mu;
        A(0,1) = lambda;
        A(0,2) = lambda;
        A(1,0) = lambda;
        A(1,1) = lambda + 2*mu;
        A(1,2) = lambda;
        A(2,0) = lambda;
        A(2,1) = lambda;
        A(2,2) = lambda + 2*mu;
        A(3,3) = mu;
        A(3,6) = mu;
        A(4,4) = mu;
        A(4,7) = mu;
        A(5,5) = mu;
        A(5,8) = mu;
        A(6,3) = mu;
        A(6,6) = mu;
        A(7,4) = mu;
        A(7,7) = mu;
        A(8,5) = mu;
        A(8,8) = mu;

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
    
    void compute_B_voigt(const double &eta,   const double &kappa, const double &nu,
                         const double &sigma, const double &tau,   Matrix_9x9 &B) {
        /*!=========================
           |    compute_B_voigt    |
           =========================
           
           Compute the B stiffness matrix in voigt notation.
        */

        B = Matrix_9x9::Zero();
        
        B(0,0) = eta + kappa + nu - 2*sigma - tau;
        B(0,1) = eta - tau;
        B(0,2) = eta - tau;
        B(1,0) = eta - tau;
        B(1,1) = eta + kappa + nu - 2*sigma - tau;
        B(1,2) = eta - tau;
        B(2,0) = eta - tau;
        B(2,1) = eta - tau;
        B(2,2) = eta + kappa + nu - 2*sigma - tau;
        B(3,3) = kappa - sigma;
        B(3,6) = nu - sigma;
        B(4,4) = kappa - sigma;
        B(4,7) = nu - sigma;
        B(5,5) = kappa - sigma;
        B(5,8) = nu - sigma;
        B(6,3) = nu - sigma;
        B(6,6) = kappa - sigma;
        B(7,4) = nu - sigma;
        B(7,7) = kappa - sigma;
        B(8,5) = nu - sigma;
        B(8,8) = kappa - sigma;
        
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
    
    void compute_C_voigt(const double &tau1,  const double &tau2,  const double &tau3,
                         const double &tau4,  const double &tau5,  const double &tau6,
                         const double &tau7,  const double &tau8,  const double &tau9,
                         const double &tau10, const double &tau11, Matrix_27x27 &C) {
        /*!=========================
           |    compute_C_voigt    |
           =========================
           
        Compute the C stiffness tensor in voigt 
        format.
        
        */
        C = Matrix_27x27::Zero();

        C( 0, 0) = 2*tau1 + tau10 + tau11 + 2*tau2 + tau3 + tau4 + 2*tau5 + tau6 + tau7 + 2*tau8 + tau9;
        C( 0, 1) = tau1 + tau4 + tau5;
        C( 0, 2) = tau1 + tau4 + tau5;
        C( 0,14) = tau2 + tau5 + tau6;
        C( 0,17) = tau1 + tau2 + tau3;
        C( 0,22) = tau2 + tau5 + tau6;
        C( 0,25) = tau1 + tau2 + tau3;
        C( 1, 0) = tau1 + tau4 + tau5;
        C( 1, 1) = tau4 + tau7 + tau9;
        C( 1, 2) = tau4;
        C( 1,14) = tau10 + tau5 + tau8;
        C( 1,17) = tau1 + tau11 + tau8;
        C( 1,22) = tau5;
        C( 1,25) = tau1;
        C( 2, 0) = tau1 + tau4 + tau5;
        C( 2, 1) = tau4;
        C( 2, 2) = tau4 + tau7 + tau9;
        C( 2,14) = tau5;
        C( 2,17) = tau1;
        C( 2,22) = tau10 + tau5 + tau8;
        C( 2,25) = tau1 + tau11 + tau8;
        C( 3, 3) = tau7;
        C( 3, 6) = tau9;
        C( 3,13) = tau10;
        C( 3,16) = tau8;
        C( 3,23) = tau8;
        C( 3,26) = tau11;
        C( 4, 4) = tau10 + tau3 + tau7;
        C( 4, 7) = tau2 + tau8 + tau9;
        C( 4,12) = tau3;
        C( 4,15) = tau2;
        C( 4,18) = tau1 + tau11 + tau8;
        C( 4,19) = tau1;
        C( 4,20) = tau1 + tau2 + tau3;
        C( 5, 5) = tau10 + tau3 + tau7;
        C( 5, 8) = tau2 + tau8 + tau9;
        C( 5, 9) = tau1 + tau11 + tau8;
        C( 5,10) = tau1 + tau2 + tau3;
        C( 5,11) = tau1;
        C( 5,21) = tau2;
        C( 5,24) = tau3;
        C( 6, 3) = tau9;
        C( 6, 6) = tau7;
        C( 6,13) = tau8;
        C( 6,16) = tau11;
        C( 6,23) = tau10;
        C( 6,26) = tau8;
        C( 7, 4) = tau2 + tau8 + tau9;
        C( 7, 7) = tau11 + tau6 + tau7;
        C( 7,12) = tau2;
        C( 7,15) = tau6;
        C( 7,18) = tau10 + tau5 + tau8;
        C( 7,19) = tau5;
        C( 7,20) = tau2 + tau5 + tau6;
        C( 8, 5) = tau2 + tau8 + tau9;
        C( 8, 8) = tau11 + tau6 + tau7;
        C( 8, 9) = tau10 + tau5 + tau8;
        C( 8,10) = tau2 + tau5 + tau6;
        C( 8,11) = tau5;
        C( 8,21) = tau6;
        C( 8,24) = tau2;
        C( 9, 5) = tau1 + tau11 + tau8;
        C( 9, 8) = tau10 + tau5 + tau8;
        C( 9, 9) = tau4 + tau7 + tau9;
        C( 9,10) = tau1 + tau4 + tau5;
        C( 9,11) = tau4;
        C( 9,21) = tau5;
        C( 9,24) = tau1;
        C(10, 5) = tau1 + tau2 + tau3;
        C(10, 8) = tau2 + tau5 + tau6;
        C(10, 9) = tau1 + tau4 + tau5;
        C(10,10) = 2*tau1 + tau10 + tau11 + 2*tau2 + tau3 + tau4 + 2*tau5 + tau6 + tau7 + 2*tau8 + tau9;
        C(10,11) = tau1 + tau4 + tau5;
        C(10,21) = tau2 + tau5 + tau6;
        C(10,24) = tau1 + tau2 + tau3;
        C(11, 5) = tau1;
        C(11, 8) = tau5;
        C(11, 9) = tau4;
        C(11,10) = tau1 + tau4 + tau5;
        C(11,11) = tau4 + tau7 + tau9;
        C(11,21) = tau10 + tau5 + tau8;
        C(11,24) = tau1 + tau11 + tau8;
        C(12, 4) = tau3;
        C(12, 7) = tau2;
        C(12,12) = tau10 + tau3 + tau7;
        C(12,15) = tau2 + tau8 + tau9;
        C(12,18) = tau1;
        C(12,19) = tau1 + tau11 + tau8;
        C(12,20) = tau1 + tau2 + tau3;
        C(13, 3) = tau10;
        C(13, 6) = tau8;
        C(13,13) = tau7;
        C(13,16) = tau9;
        C(13,23) = tau11;
        C(13,26) = tau8;
        C(14, 0) = tau2 + tau5 + tau6;
        C(14, 1) = tau10 + tau5 + tau8;
        C(14, 2) = tau5;
        C(14,14) = tau11 + tau6 + tau7;
        C(14,17) = tau2 + tau8 + tau9;
        C(14,22) = tau6;
        C(14,25) = tau2;
        C(15, 4) = tau2;
        C(15, 7) = tau6;
        C(15,12) = tau2 + tau8 + tau9;
        C(15,15) = tau11 + tau6 + tau7;
        C(15,18) = tau5;
        C(15,19) = tau10 + tau5 + tau8;
        C(15,20) = tau2 + tau5 + tau6;
        C(16, 3) = tau8;
        C(16, 6) = tau11;
        C(16,13) = tau9;
        C(16,16) = tau7;
        C(16,23) = tau8;
        C(16,26) = tau10;
        C(17, 0) = tau1 + tau2 + tau3;
        C(17, 1) = tau1 + tau11 + tau8;
        C(17, 2) = tau1;
        C(17,14) = tau2 + tau8 + tau9;
        C(17,17) = tau10 + tau3 + tau7;
        C(17,22) = tau2;
        C(17,25) = tau3;
        C(18, 4) = tau1 + tau11 + tau8;
        C(18, 7) = tau10 + tau5 + tau8;
        C(18,12) = tau1;
        C(18,15) = tau5;
        C(18,18) = tau4 + tau7 + tau9;
        C(18,19) = tau4;
        C(18,20) = tau1 + tau4 + tau5;
        C(19, 4) = tau1;
        C(19, 7) = tau5;
        C(19,12) = tau1 + tau11 + tau8;
        C(19,15) = tau10 + tau5 + tau8;
        C(19,18) = tau4;
        C(19,19) = tau4 + tau7 + tau9;
        C(19,20) = tau1 + tau4 + tau5;
        C(20, 4) = tau1 + tau2 + tau3;
        C(20, 7) = tau2 + tau5 + tau6;
        C(20,12) = tau1 + tau2 + tau3;
        C(20,15) = tau2 + tau5 + tau6;
        C(20,18) = tau1 + tau4 + tau5;
        C(20,19) = tau1 + tau4 + tau5;
        C(20,20) = 2*tau1 + tau10 + tau11 + 2*tau2 + tau3 + tau4 + 2*tau5 + tau6 + tau7 + 2*tau8 + tau9;
        C(21, 5) = tau2;
        C(21, 8) = tau6;
        C(21, 9) = tau5;
        C(21,10) = tau2 + tau5 + tau6;
        C(21,11) = tau10 + tau5 + tau8;
        C(21,21) = tau11 + tau6 + tau7;
        C(21,24) = tau2 + tau8 + tau9;
        C(22, 0) = tau2 + tau5 + tau6;
        C(22, 1) = tau5;
        C(22, 2) = tau10 + tau5 + tau8;
        C(22,14) = tau6;
        C(22,17) = tau2;
        C(22,22) = tau11 + tau6 + tau7;
        C(22,25) = tau2 + tau8 + tau9;
        C(23, 3) = tau8;
        C(23, 6) = tau10;
        C(23,13) = tau11;
        C(23,16) = tau8;
        C(23,23) = tau7;
        C(23,26) = tau9;
        C(24, 5) = tau3;
        C(24, 8) = tau2;
        C(24, 9) = tau1;
        C(24,10) = tau1 + tau2 + tau3;
        C(24,11) = tau1 + tau11 + tau8;
        C(24,21) = tau2 + tau8 + tau9;
        C(24,24) = tau10 + tau3 + tau7;
        C(25, 0) = tau1 + tau2 + tau3;
        C(25, 1) = tau1;
        C(25, 2) = tau1 + tau11 + tau8;
        C(25,14) = tau2;
        C(25,17) = tau3;
        C(25,22) = tau2 + tau8 + tau9;
        C(25,25) = tau10 + tau3 + tau7;
        C(26, 3) = tau11;
        C(26, 6) = tau8;
        C(26,13) = tau8;
        C(26,16) = tau10;
        C(26,23) = tau9;
        C(26,26) = tau7;

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


    void compute_D_voigt(const double &sigma, const double &tau, Matrix_9x9 &D){
        /*!=========================
           |    compute_D_voigt    |
           =========================
           
        Compute the D stiffness tensor in voigt 
        format.
        
        */
        D = Matrix_9x9::Zero();

        D( 0, 0) = 2*sigma + tau;
        D( 0, 1) = tau;
        D( 0, 2) = tau;
        D( 1, 0) = tau;
        D( 1, 1) = 2*sigma + tau;
        D( 1, 2) = tau;
        D( 2, 0) = tau;
        D( 2, 1) = tau;
        D( 2, 2) = 2*sigma + tau;
        D( 3, 3) = sigma;
        D( 3, 6) = sigma;
        D( 4, 4) = sigma;
        D( 4, 7) = sigma;
        D( 5, 5) = sigma;
        D( 5, 8) = sigma;
        D( 6, 3) = sigma;
        D( 6, 6) = sigma;
        D( 7, 4) = sigma;
        D( 7, 7) = sigma;
        D( 8, 5) = sigma;
        D( 8, 8) = sigma;

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
    
    void compute_PK2_stress(const Vector_9 &E_voigt,    const Vector_9 &E_micro_voigt, const Vector_27 &Gamma_voigt,
                            const Matrix_3x3 &RCGinv,   const Matrix_3x3 &Psi,         const Matrix_3x9 &Gamma,
                            const Matrix_9x9 &A,        const Matrix_9x9 &B,           const Matrix_27x27 &C,
                            const Matrix_9x9 &D,        Vector_9 &PK2){
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
    
    void compute_symmetric_stress(const Vector_9 &E_voigt,    const Vector_9 &E_micro_voigt, const Vector_27 &Gamma_voigt,
                                  const Matrix_3x3 &RCGinv,   const Matrix_3x3 &Psi,         const Matrix_3x9 &Gamma,
                                  const Matrix_9x9 &A,        const Matrix_9x9 &B,           const Matrix_27x27 &C,
                                  const Matrix_9x9 &D,        Vector_9 &SIGMA){
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
    
    void compute_higher_order_stress(const Vector_27 &Gamma_voigt, const Matrix_27x27 &C, Vector_27 &M){
        /*!=====================================
        |    compute_higher_order_stress    |
        =====================================
          
        Compute the higher order stress in the reference configuration.
          
        */
        
        M = C*Gamma_voigt; //Compute the stress (requires positive permutation)
        deformation_measures::perform_right_positive_cyclic_permutation(M); //Perform the permutation
    }

    void compute_dPK2dRCG(const Matrix_3x3 &RCG, const Matrix_3x3 &RCGinv, const Matrix_3x9 &Gamma, const Vector_27 &Gamma_voigt,
                          const Matrix_3x3 &E, const Matrix_3x3 &E_micro, const Vector_9 &E_voigt, const Vector_9 &E_micro_voigt,
                          const SpMat &A,      const SpMat &B,        const SpMat &C,  const SpMat &D, Matrix_9x9 &dPK2dRCG){
        /*!==========================
        |    compute_dPK2dRCG    |
        ==========================
        
        Compute the derivative of the PK2 stress w.r.t. 
        the deformation gradient.
        
        */
        
        //Initialize temporary terms
        Vector_9   V1;
        Vector_27  V2;
        Matrix_3x3 T1;
        Matrix_9x9 T2;
        Matrix_3x9 T3;
        
        Matrix_9x9 dRCGinvdRCG;
        deformation_measures::compute_dAinvdA(RCGinv,dRCGinvdRCG);
        
        //Compute term1
        Matrix_9x9 term1;
        term1 = 0.5*A;
        
        //Compute term2
        Matrix_9x9 term2;
        T1 = RCGinv*(E_micro+Matrix_3x3::Identity());
        deformation_measures::dot_2ot_4ot(1,1,T1,0.5*D,term2);
        
        //Compute term3
        Matrix_9x9 term3;
        V1 = (B*E_micro_voigt+D*E_voigt);
        deformation_measures::undo_voigt_3x3_tensor(V1,T1);
        deformation_measures::dot_2ot_4ot(1,0,T1*(E_micro+Matrix_3x3::Identity()).transpose(),dRCGinvdRCG,term3);
        
        //Compute term4
        Matrix_9x9 term4;
        deformation_measures::undo_voigt_3x9_tensor(C*Gamma_voigt,T3);
        T1 = T3*Gamma.transpose();
        deformation_measures::dot_2ot_4ot(1,0,T1,dRCGinvdRCG,term4);
        
        //Assemble the derivative
        dPK2dRCG = (term1+term2+term3+term4);
        
        return;
    }
    
    void compute_dPK2dRCG(const Matrix_3x3 &RCG, const Matrix_3x3 &RCGinv,  const Matrix_3x9 &Gamma, const Vector_27 &Gamma_voigt,
                          const Matrix_3x3 &E,   const Matrix_3x3 &E_micro, const Vector_9 &E_voigt, const Vector_9 &E_micro_voigt,
                          const Matrix_9x9 &A,   const Matrix_9x9 &B,       const Matrix_27x27 &C,   const Matrix_9x9 &D, Matrix_9x9 &dPK2dRCG){
        /*!==========================
        |    compute_dPK2dRCG    |
        ==========================
        
        Compute the derivative of the PK2 stress w.r.t. 
        the deformation gradient.
        
        */
        
        //Initialize temporary terms
        Vector_9   V1;
        Vector_27  V2;
        Matrix_3x3 T1;
        Matrix_9x9 T2;
        Matrix_3x9 T3;
        
        Matrix_9x9 dRCGinvdRCG;
        deformation_measures::compute_dAinvdA(RCGinv,dRCGinvdRCG);
        
        //Compute term1
        Matrix_9x9 term1;
        term1 = 0.5*A;
        
        //Compute term2
        Matrix_9x9 term2;
        T1 = RCGinv*(E_micro+Matrix_3x3::Identity());
        deformation_measures::dot_2ot_4ot(1,1,T1,0.5*D,term2);
        
        //Compute term3
        Matrix_9x9 term3;
        V1 = (B*E_micro_voigt+D*E_voigt);
        deformation_measures::undo_voigt_3x3_tensor(V1,T1);
        deformation_measures::dot_2ot_4ot(1,0,T1*(E_micro+Matrix_3x3::Identity()).transpose(),dRCGinvdRCG,term3);
        
        //Compute term4
        Matrix_9x9 term4;
        deformation_measures::undo_voigt_3x9_tensor(C*Gamma_voigt,T3);
        T1 = T3*Gamma.transpose();
        deformation_measures::dot_2ot_4ot(1,0,T1,dRCGinvdRCG,term4);
        
        //Assemble the derivative
        dPK2dRCG = (term1+term2+term3+term4);
        
        return;
    }

    void compute_dPK2dRCG(const Matrix_3x3 &RCG, const Matrix_3x3 &RCGinv, const Matrix_3x9 &Gamma, const Vector_27 &Gamma_voigt,
                          const Matrix_3x3 &E, const Matrix_3x3 &E_micro, const Vector_9 &E_voigt, const Vector_9 &E_micro_voigt,
                          const SpMat &A,      const SpMat &B,        const SpMat &C,  const SpMat &D, Matrix_9x9 (&terms)[2], Matrix_9x9 &dPK2dRCG){
        /*!==========================
        |    compute_dPK2dRCG    |
        ==========================
        
        Compute the derivative of the PK2 stress w.r.t. 
        the deformation gradient.
        
        */
        
        //Initialize temporary terms
        Vector_9   V1;
        Vector_27  V2;
        Matrix_3x3 T1;
        Matrix_9x9 T2;
        Matrix_3x9 T3;
        Matrix_9x9 temp;
        
        Matrix_9x9 dRCGinvdRCG;
        deformation_measures::compute_dAinvdA(RCGinv,dRCGinvdRCG);
        
        //Compute term1
        terms[0] = 0.5*A;
        
        //Compute term2
        T1 = RCGinv*(E_micro+Matrix_3x3::Identity());
        deformation_measures::dot_2ot_4ot(1,1,T1,0.5*D,terms[1]);
        
        //Compute term3
        V1 = (B*E_micro_voigt+D*E_voigt);
        deformation_measures::undo_voigt_3x3_tensor(V1,T1);
        deformation_measures::dot_2ot_4ot(1,0,T1*(E_micro+Matrix_3x3::Identity()).transpose(),dRCGinvdRCG,temp);
        terms[1] += temp;        

        //Compute term4
        deformation_measures::undo_voigt_3x9_tensor(C*Gamma_voigt,T3);
        T1 = T3*Gamma.transpose();
        deformation_measures::dot_2ot_4ot(1,0,T1,dRCGinvdRCG,temp);
        terms[1] += temp;        

        //Assemble the derivative
        dPK2dRCG = (terms[0] + terms[1]);
        
        return;
    }
    
    void compute_dPK2dRCG(const Matrix_3x3 &RCG, const Matrix_3x3 &RCGinv,  const Matrix_3x9 &Gamma, const Vector_27 &Gamma_voigt,
                          const Matrix_3x3 &E,   const Matrix_3x3 &E_micro, const Vector_9 &E_voigt, const Vector_9 &E_micro_voigt,
                          const Matrix_9x9 &A,   const Matrix_9x9 &B,       const Matrix_27x27 &C,   const Matrix_9x9 &D, Matrix_9x9 (&terms)[2], Matrix_9x9 &dPK2dRCG){
        /*!==========================
        |    compute_dPK2dRCG    |
        ==========================
        
        Compute the derivative of the PK2 stress w.r.t. 
        the deformation gradient.
        
        */
        
        //Initialize temporary terms
        Vector_9   V1;
        Vector_27  V2;
        Matrix_3x3 T1;
        Matrix_9x9 T2;
        Matrix_3x9 T3;
        Matrix_9x9 temp;
        
        Matrix_9x9 dRCGinvdRCG;
        deformation_measures::compute_dAinvdA(RCGinv,dRCGinvdRCG);
        
        //Compute term1
        terms[0] = 0.5*A;
        
        //Compute term2
        T1 = RCGinv*(E_micro+Matrix_3x3::Identity());
        deformation_measures::dot_2ot_4ot(1,1,T1,0.5*D,terms[1]);
        
        //Compute term3
        V1 = (B*E_micro_voigt+D*E_voigt);
        deformation_measures::undo_voigt_3x3_tensor(V1,T1);
        deformation_measures::dot_2ot_4ot(1,0,T1*(E_micro+Matrix_3x3::Identity()).transpose(),dRCGinvdRCG,temp);
        terms[1] += temp;
        
        //Compute term4
        deformation_measures::undo_voigt_3x9_tensor(C*Gamma_voigt,T3);
        T1 = T3*Gamma.transpose();
        deformation_measures::dot_2ot_4ot(1,0,T1,dRCGinvdRCG,temp);
        terms[1] += temp;
        
        //Assemble the derivative
        dPK2dRCG = (terms[0] + terms[1]);
        
        return;
    }

    void compute_dPK2dPsi(const Matrix_3x3 &RCGinv, const Matrix_3x3 &E_micro, const Vector_9 &E_voigt, const Vector_9 &E_micro_voigt,
                          const SpMat &B, const SpMat &D, Matrix_9x9 &dPK2dPsi){
        /*!==========================
        |    compute_dPK2dPsi    |
        ==========================
        
        Compute the derivative of the second piola kirchoff 
        stress w.r.t. the deformation measure Psi.
        
        */
        
        //Initialize temporary terms
        Vector_9   V1;
        Vector_27  V2;
        Matrix_3x3 T1;
        Matrix_9x9 T2;
        Matrix_3x9 T3;
        
        //Add term1
        dPK2dPsi = D;
        
        //Add term2
        deformation_measures::dot_2ot_4ot(1, 1, RCGinv*(E_micro + Matrix_3x3::Identity()), B, T2);
        dPK2dPsi += T2;
        
        //Add term3
        V1 = B*E_micro_voigt + D*E_voigt;
        deformation_measures::undo_voigt_3x3_tensor(V1,T1);
        deformation_measures::two_sot_to_fot(2,T1,RCGinv,T2);
        dPK2dPsi += T2;
        
        return;
    }
    
    void compute_dPK2dPsi(const Matrix_3x3 &RCGinv, const Matrix_3x3 &E_micro, const Vector_9 &E_voigt, const Vector_9 &E_micro_voigt,
                          const Matrix_9x9 &B,      const Matrix_9x9 &D,       Matrix_9x9 &dPK2dPsi){
        /*!==========================
        |    compute_dPK2dPsi    |
        ==========================
        
        Compute the derivative of the second piola kirchoff 
        stress w.r.t. the deformation measure Psi.
        
        */
        
        //Initialize temporary terms
        Vector_9   V1;
        Vector_27  V2;
        Matrix_3x3 T1;
        Matrix_9x9 T2;
        Matrix_3x9 T3;
        
        //Add term1
        dPK2dPsi = D;
        
        //Add term2
        deformation_measures::dot_2ot_4ot(1, 1, RCGinv*(E_micro + Matrix_3x3::Identity()), B, T2);
        dPK2dPsi += T2;
        
        //Add term3
        V1 = B*E_micro_voigt + D*E_voigt;
        deformation_measures::undo_voigt_3x3_tensor(V1,T1);
        deformation_measures::two_sot_to_fot(2,T1,RCGinv,T2);
        dPK2dPsi += T2;
        
        return;
    }

    void compute_dPK2dPsi(const Matrix_3x3 &RCGinv, const Matrix_3x3 &E_micro, const Vector_9 &E_voigt, const Vector_9 &E_micro_voigt,
                          const SpMat &B, const SpMat &D, Matrix_9x9 (&terms)[2], Matrix_9x9 &dPK2dPsi){
        /*!==========================
        |    compute_dPK2dPsi    |
        ==========================
        
        Compute the derivative of the second piola kirchoff 
        stress w.r.t. the deformation measure Psi.
        
        */
        
        //Initialize temporary terms
        Vector_9   V1;
        Vector_27  V2;
        Matrix_3x3 T1;
        Matrix_9x9 T2;
        Matrix_3x9 T3;
        
        //Add term1
        terms[0] = D;
        
        //Add term2
        deformation_measures::dot_2ot_4ot(1, 1, RCGinv*(E_micro + Matrix_3x3::Identity()), B, T2);
        terms[1] = T2;
        
        //Add term3
        V1 = B*E_micro_voigt + D*E_voigt;
        deformation_measures::undo_voigt_3x3_tensor(V1,T1);
        deformation_measures::two_sot_to_fot(2,T1,RCGinv,T2);
        terms[1] += T2;
        
        dPK2dPsi = terms[0] + terms[1];
        
        return;
    }
    
    void compute_dPK2dPsi(const Matrix_3x3 &RCGinv, const Matrix_3x3 &E_micro, const Vector_9 &E_voigt, const Vector_9 &E_micro_voigt,
                          const Matrix_9x9 &B,      const Matrix_9x9 &D,       Matrix_9x9 (&terms)[2],  Matrix_9x9 &dPK2dPsi){
        /*!==========================
        |    compute_dPK2dPsi    |
        ==========================
        
        Compute the derivative of the second piola kirchoff 
        stress w.r.t. the deformation measure Psi.
        
        */
        
        //Initialize temporary terms
        Vector_9   V1;
        Vector_27  V2;
        Matrix_3x3 T1;
        Matrix_9x9 T2;
        Matrix_3x9 T3;
        
        //Add term1
        terms[0] = D;
        
        //Add term2
        deformation_measures::dot_2ot_4ot(1, 1, RCGinv*(E_micro + Matrix_3x3::Identity()), B, T2);
        terms[1] = T2;
        
        //Add term3
        V1 = B*E_micro_voigt + D*E_voigt;
        deformation_measures::undo_voigt_3x3_tensor(V1,T1);
        deformation_measures::two_sot_to_fot(2,T1,RCGinv,T2);
        terms[1] += T2;
        
        dPK2dPsi = terms[0] + terms[1];
        
        return;
    }

    void compute_dPK2dGamma_term2(const Vector_27 &term1, const Matrix_27x27 &C, Matrix_9x27 &term2){
        /*!==================================
        |    compute_dPK2dGamma_term2    |
        ==================================
        
        Compute term2 of dPK2dGamma. Note that this is identical as 
        term2 for dSIGMAdGamma if symmetrized.
        
        */

        term2 = Matrix_9x27::Zero();

        int sot_to_voigt_map[3][3] = {{0,5,4},
                                      {8,1,3},
                                      {7,6,2}};
        
        int tot_to_voigt_map[3][3][3];
        deformation_measures::get_tot_to_voigt_map(tot_to_voigt_map);

        int Ihat;
        int Jhat;
        int Khat;
        int Lhat;

        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++){
                Ihat = sot_to_voigt_map[i][j];
                for (int t=0; t<3; t++){
                    for (int u=0; u<3; u++){
                        for (int v=0; v<3; v++){
                            Jhat = tot_to_voigt_map[t][u][v];

                            for (int q=0; q<3; q++){
                                for (int r=0; r<3; r++){
                                    Khat = tot_to_voigt_map[i][q][r];
                                    Lhat = tot_to_voigt_map[j][q][r];
                                    term2(Ihat,Jhat) += C(Khat,Jhat)*term1(Lhat);
                                }
                            }
                        }
                    }
                }
            }
        }
        return;
    }
    
    void compute_dPK2dGamma_term3(const Vector_27 &term1, const Matrix_3x3 &RCGinv, Matrix_9x27 &term3){
        /*!==================================
        |    compute_dPK2dGamma_term3    |
        ==================================
        
        Compute term3 of dPK2dGamma. Note, also the same as dSIGMAdGamma if symmetrized.
        
        */
        
        int sot_to_voigt_map[3][3] = {{0,5,4},
                                      {8,1,3},
                                      {7,6,2}};
        
        int tot_to_voigt_map[3][3][3];
        deformation_measures::get_tot_to_voigt_map(tot_to_voigt_map);

        int Ihat;
        int Jhat;
        int Khat;
        int Lhat;

        double tmp1;
        
        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++){
                Ihat = sot_to_voigt_map[i][j];
                for (int t=0; t<3; t++){
                    tmp1 = RCGinv(j,t);
                    for (int u=0; u<3; u++){
                        for (int v=0; v<3; v++){
                            Jhat = tot_to_voigt_map[t][u][v];
                            Khat = tot_to_voigt_map[i][u][v];
                            term3(Ihat,Jhat) = term1(Khat)*tmp1;
                        }
                    }
                }
            }
        }
        return;
    }
    
    void compute_dPK2dGamma(const Matrix_3x3 &RCGinv, const Matrix_3x9 &Gamma, const Vector_27 &Gamma_voigt,
                            SpMat &C, Matrix_9x27 &dPK2dGamma){
        /*!============================
        |    compute_dPK2dGamma    |
        ============================
        
        Compute the derivative of the second piola kirchoff stress 
        w.r.t. the deformation measure Gamma.
        
        */
        
        //Compute term1
        Vector_27 term1;
        deformation_measures::voigt_3x9_tensor(RCGinv*Gamma,term1);
        
        //Compute term2
        Matrix_9x27 term2;
        Matrix_27x27 _C = C; //Copy the sparse matrix to a dense matrix.
        compute_dPK2dGamma_term2(term1,_C,term2);
        
        //Compute term3
        Matrix_9x27 term3;
        compute_dPK2dGamma_term3(C*Gamma_voigt,RCGinv,term3);
        
        dPK2dGamma = term2 + term3;
    }
    
    void compute_dPK2dGamma(const Matrix_3x3 &RCGinv, const Matrix_3x9 &Gamma, const Vector_27 &Gamma_voigt,
                            Matrix_27x27 &C, Matrix_9x27 &dPK2dGamma){
        /*!============================
        |    compute_dPK2dGamma    |
        ============================
        
        Compute the derivative of the second piola kirchoff stress 
        w.r.t. the deformation measure Gamma.
        
        */
        
        //Compute term1
        Vector_27 term1;
        deformation_measures::voigt_3x9_tensor(RCGinv*Gamma,term1);
        
        //Compute term2
        Matrix_9x27 term2;
        compute_dPK2dGamma_term2(term1, C,term2);
        
        //Compute term3
        Matrix_9x27 term3;
        compute_dPK2dGamma_term3(C*Gamma_voigt,RCGinv,term3);
        
        dPK2dGamma = term2 + term3;
    }

    void compute_dPK2dGamma(const Matrix_3x3 &RCGinv, const Matrix_3x9 &Gamma, const Vector_27 &Gamma_voigt,
                            SpMat &C, Matrix_9x27 (&terms)[2], Matrix_9x27 &dPK2dGamma){
        /*!============================
        |    compute_dPK2dGamma    |
        ============================
        
        Compute the derivative of the second piola kirchoff stress 
        w.r.t. the deformation measure Gamma.
        
        */
        
        //Compute term1
        Vector_27 term1;
        deformation_measures::voigt_3x9_tensor(RCGinv*Gamma,term1);
        
        //Compute term2
        Matrix_27x27 _C = C; //Copy the sparse matrix to a dense matrix.
        compute_dPK2dGamma_term2(term1,_C,terms[0]);
        
        //Compute term3
        compute_dPK2dGamma_term3(C*Gamma_voigt,RCGinv,terms[1]);
        
        dPK2dGamma = terms[0] + terms[1];
    }
    
    void compute_dPK2dGamma(const Matrix_3x3 &RCGinv, const Matrix_3x9 &Gamma, const Vector_27 &Gamma_voigt,
                            Matrix_27x27 &C, Matrix_9x27 (&terms)[2], Matrix_9x27 &dPK2dGamma){
        /*!============================
        |    compute_dPK2dGamma    |
        ============================
        
        Compute the derivative of the second piola kirchoff stress 
        w.r.t. the deformation measure Gamma.
        
        */
        
        //Compute term1
        Vector_27 term1;
        deformation_measures::voigt_3x9_tensor(RCGinv*Gamma,term1);
        
        //Compute term2
        compute_dPK2dGamma_term2(term1,C,terms[0]);
        
        //Compute term3
        compute_dPK2dGamma_term3(C*Gamma_voigt,RCGinv,terms[1]);
        
        dPK2dGamma = terms[0] + terms[1];
    }

    void compute_dSIGMAdRCG(Matrix_9x9 (&terms)[2], Matrix_9x9 &dSIGMAdRCG){
        /*!=========================
        |    compute_dSIGMAdRCG    |
        ============================
        
        Compute the derivative of the symmetric stress in the 
        reference configuration using the terms from the 
        computation of dPK2dRCG.
        
        This is the preferred method since the computation has 
        already been done for the PK2 derivative.
        
        */
        
        dSIGMAdRCG = terms[0];
        
        dSIGMAdRCG        += terms[1];
        dSIGMAdRCG.row(0) += terms[1].row(0);
        dSIGMAdRCG.row(1) += terms[1].row(1);
        dSIGMAdRCG.row(2) += terms[1].row(2);
        
        dSIGMAdRCG.row(3) += terms[1].row(6);
        dSIGMAdRCG.row(4) += terms[1].row(7);
        dSIGMAdRCG.row(5) += terms[1].row(8);
        dSIGMAdRCG.row(6) += terms[1].row(3);
        dSIGMAdRCG.row(7) += terms[1].row(4);
        dSIGMAdRCG.row(8) += terms[1].row(5);
 
//        for (int i=1; i<4; i++){
//            temp = terms[i];
//            
//            temp.row(0) *= 2;
//            temp.row(1) *= 2;
//            temp.row(2) *= 2;
//            temp.row(3) += terms[i].row(6);
//            temp.row(4) += terms[i].row(7);
//            temp.row(5) += terms[i].row(8);
//            temp.row(6) += terms[i].row(3);
//            temp.row(7) += terms[i].row(4);
//            temp.row(8) += terms[i].row(5);
//            
//            dSIGMAdRCG += temp;
//        }
        return;
    }
    
    void compute_dSIGMAdPsi(Matrix_9x9 (&terms)[2], Matrix_9x9 &dSIGMAdPsi){
        /*!============================
        |    compute_dSIGMAdPsi    |
        ============================
        
        Compute the derivative of the symmetric stress in the 
        reference configuration using the terms from the 
        computation of dPK2dPsi.
        
        This is the preferred method since the computation has already 
        been done for the PK2 derivative.
        
        */
        
        dSIGMAdPsi = terms[0]+terms[1];
        
        dSIGMAdPsi.row(0) += terms[1].row(0);
        dSIGMAdPsi.row(1) += terms[1].row(1);
        dSIGMAdPsi.row(2) += terms[1].row(2);
        
        dSIGMAdPsi.row(3) += terms[1].row(6);
        dSIGMAdPsi.row(4) += terms[1].row(7);
        dSIGMAdPsi.row(5) += terms[1].row(8);
        dSIGMAdPsi.row(6) += terms[1].row(3);
        dSIGMAdPsi.row(7) += terms[1].row(4);
        dSIGMAdPsi.row(8) += terms[1].row(5);
        
//        for (int i=1; i<3; i++){
//            temp = terms[i];
//            
//            temp.row(0) *= 2;
//            temp.row(1) *= 2;
//            temp.row(2) *= 2;
//            
//            temp.row(3) += terms[i].row(6);
//            temp.row(4) += terms[i].row(7);
//            temp.row(5) += terms[i].row(8);
//            temp.row(6) += terms[i].row(3);
//            temp.row(7) += terms[i].row(4);
//            temp.row(8) += terms[i].row(5);
//            
//            dSIGMAdPsi += temp;
//        }
        return;
    }
    
    void compute_dSIGMAdGamma(Matrix_9x27 (&terms)[2], Matrix_9x27 &dSIGMAdGamma){
        /*!===========================
        |    compute_dSIGMAdGamma    |
        ==============================
        
        Compute the derivative of the symmetric stress in the 
        reference configuration using the terms from the computation 
        of dPK2dGamma.
        
        This is the preferred method since the computation has already
        been done for the PK2 derivative.
        
        */
        
        dSIGMAdGamma = Matrix_9x27::Zero();
        
        Matrix_9x27 temp;
        
        for (int i=0; i<2; i++){
            temp = terms[i];
            
            temp.row(0) *= 2.;
            temp.row(1) *= 2.;
            temp.row(2) *= 2.;
            temp.row(3) += terms[i].row(6);
            temp.row(4) += terms[i].row(7);
            temp.row(5) += terms[i].row(8);
            temp.row(6) += terms[i].row(3);
            temp.row(7) += terms[i].row(4);
            temp.row(8) += terms[i].row(5);
            
            dSIGMAdGamma += temp;
        }
        return;
    }
    
    void compute_dMdGamma(const Matrix_27x27 &C, Matrix_27x27 &dMdGamma){
        /*!==========================
        |    compute_dMdGamma    |
        ==========================
        
        Compute the derivative of the higher order stress tensor w.r.t. Gamma.
        
        */

        int tot_to_voigt_map[3][3][3];
        deformation_measures::get_tot_to_voigt_map(tot_to_voigt_map);

        int Ihat;
        int Jhat;
        int Khat;

        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++){
                for (int k=0; k<3; k++){
                    Ihat = tot_to_voigt_map[i][j][k];
                    Khat = tot_to_voigt_map[j][k][i];
                    for (int l=0; l<3; l++){
                        for (int m=0; m<3; m++){
                            for (int n=0; n<3; n++){
                                Jhat = tot_to_voigt_map[l][m][n];
                                dMdGamma(Ihat,Jhat) = C(Khat,Jhat);
                            }
                        }
                    }
                }
            }
        }
        return;
    }
}
