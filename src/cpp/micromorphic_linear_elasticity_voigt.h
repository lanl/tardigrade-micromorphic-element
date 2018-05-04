/*!============================================================
  |                                                           |
  |         micromorphic_linear_elasticity_voigt.h            |
  |                                                           |
  -------------------------------------------------------------
  | The header file for the definition of a                   |
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
  
#ifndef MICROMORPHIC_LINEAR_ELASTICITY_VOIGT_H
#define MICROMORPHIC_LINEAR_ELASTICITY_VOIGT_H
  
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <micromorphic_material_library.h>

//Dense matrix type definitions
typedef Eigen::Matrix<double,3,3> Matrix_3x3;
typedef Eigen::Matrix<double,3,9> Matrix_3x9;

//Sparse matrix type definitions
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

namespace micro_material{
    
    class LinearElasticity: public micromorphic_material_library::IMaterial{
        /*!
        ==========================
        |    LinearElasticity    |
        ==========================
        
        The class which is called when evaluating a 
        linear elastic micromorphic constitutive model.
        */
        
        public:
            void evaluate_model(const std::vector<double> &time,        const std::vector<double> (&fparams),
                                const double (&grad_u)[3][3],           const double (&phi)[9],
                                const double (&grad_phi)[9][3],         std::vector<double> &SDVS,
                                const std::vector<double> &ADD_DOF,     const std::vector<std::vector<double>> &ADD_grad_DOF,
                                Vector_9 &cauchy, Vector_9 &s, Vector_27 &m, std::vector<Eigen::VectorXd> &ADD_TERMS);

            void evaluate_model(const std::vector<double> &time,        const std::vector<double> (&fparams),
                                const double (&grad_u)[3][3],           const double (&phi)[9],
                                const double (&grad_phi)[9][3],         std::vector<double> &SDVS,
                                const std::vector<double> &ADD_DOF,     const std::vector<std::vector<double>> &ADD_grad_DOF,
                                Vector_9    &cauchy,    Vector_9    &s,           Vector_27    &m,
                                Matrix_9x9  &DcauchyDgrad_u, Matrix_9x9  &DcauchyDphi, Matrix_9x27  &DcauchyDgrad_phi,
                                Matrix_9x9  &DsDgrad_u,      Matrix_9x9  &DsDphi,      Matrix_9x27  &DsDgrad_phi,
                                Matrix_27x9 &DmDgrad_u,      Matrix_27x9 &DmDphi,      Matrix_27x27 &DmDgrad_phi,
                                std::vector<Eigen::VectorXd> &ADD_TERMS,               std::vector<Eigen::MatrixXd> &ADD_JACOBIANS);
                                
            void get_deformation_measures(const double (&grad_u)[3][3], const double (&phi)[9], const double (&grad_phi)[9][3],
                                          Matrix_3x3 &F,          Matrix_3x3 &chi,  Matrix_3x9 &grad_chi);
    };

    void get_stress(const double &t,     const double &dt,      const double (&params)[18], 
                    const Matrix_3x3 &F, const Matrix_3x3 &chi, const Matrix_3x9 &grad_chi,
                    std::vector<double> &SDVS, Vector_9 &PK2, Vector_9 &SIGMA, Vector_27 &M);
                    
    void get_stress(const double &t,     const double &dt,      const double (&params)[18], 
                    const Matrix_3x3 &F, const Matrix_3x3 &chi, const Matrix_3x9 &grad_chi,
                    std::vector<double> &SDVS, Vector_9 &PK2, Vector_9 &SIGMA, Vector_27 &M,
                    Matrix_9x9  &dPK2dF,   Matrix_9x9  &dPK2dchi,   Matrix_9x27  &dPK2dgrad_chi,
                    Matrix_9x9  &dSIGMAdF, Matrix_9x9  &dSIGMAdchi, Matrix_9x27  &dSIGMAdgrad_chi,
                    Matrix_27x9 &dMdF,     Matrix_27x9 &dMdchi,     Matrix_27x27 &dMdgrad_chi);

    void compute_A_voigt(const double &lambda,const double &mu, SpMat &A);

    void compute_A_voigt(const double &lambda,const double &mu, Matrix_9x9 &A);

    void compute_B_voigt(const double &eta, const double &kappa,
                         const double &nu,  const double &sigma,
                         const double &tau, SpMat &B);

    void compute_B_voigt(const double &eta, const double &kappa,
                         const double &nu,  const double &sigma,
                         const double &tau, Matrix_9x9 &B);

    void compute_C_voigt(const double &tau1,  const double &tau2,  const double &tau3,
                         const double &tau4,  const double &tau5,  const double &tau6,
                         const double &tau7,  const double &tau8,  const double &tau9,
                         const double &tau10, const double &tau11, SpMat &C);

    void compute_C_voigt(const double &tau1,  const double &tau2,  const double &tau3,
                         const double &tau4,  const double &tau5,  const double &tau6,
                         const double &tau7,  const double &tau8,  const double &tau9,
                         const double &tau10, const double &tau11, Matrix_27x27 &C);

    void compute_D_voigt(const double &sigma, const double &tau, SpMat &D);

    void compute_D_voigt(const double &sigma, const double &tau, Matrix_9x9 &D);

    void compute_PK2_stress(const Vector_9 &E_voigt,      const Vector_9 &E_micro_voigt,
                            const Vector_27 &Gamma_voigt, const Matrix_3x3 &RCGinv,
                            const Matrix_3x3 &Psi,        const Matrix_3x9 &Gamma,
                            const SpMat &A,               const SpMat &B,
                            const SpMat &C,               const SpMat &D,
                            Vector_9 &PK2);
                            
    void compute_PK2_stress(const Vector_9 &E_voigt,      const Vector_9 &E_micro_voigt,
                            const Vector_27 &Gamma_voigt, const Matrix_3x3 &RCGinv,
                            const Matrix_3x3 &Psi,        const Matrix_3x9 &Gamma,
                            const Matrix_9x9 &A,          const Matrix_9x9 &B,
                            const Matrix_27x27 &C,        const Matrix_9x9 &D,
                            Vector_9 &PK2);

    void compute_symmetric_stress(const Vector_9 &E_voigt,      const Vector_9 &E_micro_voigt,
                                  const Vector_27 &Gamma_voigt, const Matrix_3x3 &RCGinv,
                                  const Matrix_3x3 &Psi,        const Matrix_3x9 &Gamma,
                                  const SpMat &A,               const SpMat &B,
                                  const SpMat &C,               const SpMat &D,
                                  Vector_9 &SIGMA);

    void compute_symmetric_stress(const Vector_9 &E_voigt,      const Vector_9 &E_micro_voigt,
                                  const Vector_27 &Gamma_voigt, const Matrix_3x3 &RCGinv,
                                  const Matrix_3x3 &Psi,        const Matrix_3x9 &Gamma,
                                  const Matrix_9x9 &A,          const Matrix_9x9 &B,
                                  const Matrix_27x27 &C,        const Matrix_9x9 &D,
                                  Vector_9 &SIGMA);
                                  
    void compute_higher_order_stress(const Vector_27 &Gamma_voigt, const SpMat &C, Vector_27 &M);

    void compute_higher_order_stress(const Vector_27 &Gamma_voigt, const Matrix_27x27 &C, Vector_27 &M);

    //Use internal terms in computation of tangents
    void compute_dPK2dRCG(const Matrix_3x3 &RCG, const Matrix_3x3 &RCGinv, const Matrix_3x9 &Gamma, const Vector_27 &Gamma_voigt,
                        const Matrix_3x3 &E, const Matrix_3x3 &E_micro, const Vector_9 &E_voigt, const Vector_9 &E_micro_voigt,
                        const SpMat &A,      const SpMat &B,        const SpMat &C,  const SpMat &D, Matrix_9x9 &dPK2dRCG);
                        
    void compute_dPK2dRCG(const Matrix_3x3 &RCG, const Matrix_3x3 &RCGinv,    const Matrix_3x9 &Gamma, const Vector_27 &Gamma_voigt,
                          const Matrix_3x3 &E,     const Matrix_3x3 &E_micro, const Vector_9 &E_voigt, const Vector_9 &E_micro_voigt,
                          const Matrix_9x9 &A,     const Matrix_9x9 &B,       const Matrix_27x27 &C,   const Matrix_9x9 &D, Matrix_9x9 &dPK2dRCG);

    void compute_dPK2dPsi(const Matrix_3x3 &RCGinv, const Matrix_3x3 &E_micro, const Vector_9 &E_voigt, const Vector_9 &E_micro_voigt,
                          const SpMat &B, const SpMat &D, Matrix_9x9 &dPK2dPsi);
                          
    void compute_dPK2dPsi(const Matrix_3x3 &RCGinv, const Matrix_3x3 &E_micro, const Vector_9 &E_voigt, const Vector_9 &E_micro_voigt,
                          const Matrix_9x9 &B,      const Matrix_9x9 &D, Matrix_9x9 &dPK2dPsi);

    void compute_dPK2dGamma(const Matrix_3x3 &RCGinv, const Matrix_3x9 &Gamma, const Vector_27 &Gamma_voigt,
                            SpMat &C, Matrix_9x27 &dPK2dGamma);
    
    void compute_dPK2dGamma(const Matrix_3x3 &RCGinv, const Matrix_3x9 &Gamma, const Vector_27 &Gamma_voigt,
                            Matrix_27x27 &C, Matrix_9x27 &dPK2dGamma);

    //Use external terms in computation of tangents (preferred)
    
    //Gradients of the PK2 stress
    void compute_dPK2dRCG(const Matrix_3x3 &RCG, const Matrix_3x3 &RCGinv, const Matrix_3x9 &Gamma, const Vector_27 &Gamma_voigt,
                          const Matrix_3x3 &E, const Matrix_3x3 &E_micro, const Vector_9 &E_voigt, const Vector_9 &E_micro_voigt,
                          const SpMat &A,      const SpMat &B,        const SpMat &C,  const SpMat &D, Matrix_9x9 (&terms)[2], Matrix_9x9 &dPK2dRCG);
                        
    void compute_dPK2dRCG(const Matrix_3x3 &RCG, const Matrix_3x3 &RCGinv,  const Matrix_3x9 &Gamma, const Vector_27 &Gamma_voigt,
                          const Matrix_3x3 &E,   const Matrix_3x3 &E_micro, const Vector_9 &E_voigt, const Vector_9 &E_micro_voigt,
                          const Matrix_9x9 &A,   const Matrix_9x9 &B,       const Matrix_27x27 &C,   const Matrix_9x9 &D, Matrix_9x9 (&terms)[2], Matrix_9x9 &dPK2dRCG);

    void compute_dPK2dPsi(const Matrix_3x3 &RCGinv, const Matrix_3x3 &E_micro, const Vector_9 &E_voigt, const Vector_9 &E_micro_voigt,
                          const SpMat &B, const SpMat &D, Matrix_9x9 (&terms)[2], Matrix_9x9 &dPK2dPsi);
                          
    void compute_dPK2dPsi(const Matrix_3x3 &RCGinv, const Matrix_3x3 &E_micro, const Vector_9 &E_voigt, const Vector_9 &E_micro_voigt,
                          const Matrix_9x9 &B,      const Matrix_9x9 &D,       Matrix_9x9 (&terms)[2],  Matrix_9x9 &dPK2dPsi);

    void compute_dPK2dGamma(const Matrix_3x3 &RCGinv, const Matrix_3x9 &Gamma, const Vector_27 &Gamma_voigt,
                            SpMat &C, Matrix_9x27 (&terms)[2], Matrix_9x27 &dPK2dGamma);
                            
    void compute_dPK2dGamma(const Matrix_3x3 &RCGinv, const Matrix_3x9 &Gamma, const Vector_27 &Gamma_voigt,
                            Matrix_27x27 &C,          Matrix_9x27 (&terms)[2], Matrix_9x27 &dPK2dGamma);

    //Gradients of the symmetric stress (reference configuration)
    void compute_dSIGMAdRCG(Matrix_9x9 (&terms)[2], Matrix_9x9 &dSIGMAdRCG);
    
    void compute_dSIGMAdPsi(Matrix_9x9 (&terms)[2], Matrix_9x9 &dSIGMAdPsi);
    
    void compute_dSIGMAdGamma(Matrix_9x27 (&terms)[2], Matrix_9x27 &dSIGMAdGamma);
    
    void compute_dMdGamma(const Matrix_27x27 &C, Matrix_27x27 &dMdGamma);
                          
    void compute_dPK2dGamma_term2(const Vector_27 &term1, const Matrix_27x27 &C, Matrix_9x27 &term2);
    
    void compute_dPK2dGamma_term3(const Vector_27 &term1, const Matrix_3x3 &RCGinv, Matrix_9x27 &term3);

    REGISTER_MATERIAL(LinearElasticity)
}

#endif
