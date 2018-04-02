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

    void get_stress(double (&params)[18]);

    void compute_A_voigt(const double &lambda,const double &mu, SpMat &A);

    void compute_B_voigt(const double &eta, const double &kappa,
                         const double &nu,  const double &sigma,
                         const double &tau, SpMat &B);

    void compute_C_voigt(const double &tau1,  const double &tau2,  const double &tau3,
                         const double &tau4,  const double &tau5,  const double &tau6,
                         const double &tau7,  const double &tau8,  const double &tau9,
                         const double &tau10, const double &tau11, SpMat &C);

    void compute_D_voigt(const double &sigma, const double &tau, SpMat &D);
}
