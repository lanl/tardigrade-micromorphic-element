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
typedef Eigen::SparseMatrix<double,9,6> Sparse_Matrix_9x6;
typedef Eigen::Triplet<double> T;

namespace micro_material{

    void compute_A_voigt(const double &lambda,const double &mu, Sparse_Matrix_9x6 &A);

}