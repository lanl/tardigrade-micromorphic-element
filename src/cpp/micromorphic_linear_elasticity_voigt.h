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

//Dense matrix type definitions
typedef Eigen::Matrix<double,3,3> Matrix_3x3;
typedef Eigen::Matrix<double,3,9> Matrix_3x9;

//Sparse matrix type definitions
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

namespace micro_material{

    void get_stress(const double t,      const double dt,       const double (&params)[18], 
                    const Matrix_3x3 &F, const Matrix_3x3 &chi, const Matrix_3x9 &grad_chi,
                    std::vector<double> &SDVS, Vector_9 &PK2, Vector_9 &SIGMA, Vector_27 &M);

    void compute_A_voigt(const double &lambda,const double &mu, SpMat &A);

    void compute_B_voigt(const double &eta, const double &kappa,
                         const double &nu,  const double &sigma,
                         const double &tau, SpMat &B);

    void compute_C_voigt(const double &tau1,  const double &tau2,  const double &tau3,
                         const double &tau4,  const double &tau5,  const double &tau6,
                         const double &tau7,  const double &tau8,  const double &tau9,
                         const double &tau10, const double &tau11, SpMat &C);

    void compute_D_voigt(const double &sigma, const double &tau, SpMat &D);

    void compute_PK2_stress(const Vector_9 &E_voigt,      const Vector_9 &E_micro_voigt,
                            const Vector_27 &Gamma_voigt, const Matrix_3x3 &Cinv,
                            const Matrix_3x3 &Psi,        const Matrix_3x9 &Gamma,
                            const SpMat &A,               const SpMat &B,
                            const SpMat &C,               const SpMat &D,
                            Vector_9 &PK2);
                            
    void compute_symmetric_stress(const Vector_9 &E_voigt,      const Vector_9 &E_micro_voigt,
                                  const Vector_27 &Gamma_voigt, const Matrix_3x3 &Cinv,
                                  const Matrix_3x3 &Psi,        const Matrix_3x9 &Gamma,
                                  const SpMat &A,               const SpMat &B,
                                  const SpMat &C,               const SpMat &D,
                                  Vector_9 &SIGMA);
                                  
    void compute_higher_order_stress(const Vector_27 &Gamma_voigt, const SpMat &C, Vector_27 &M);
}

#endif