/*!============================================================
  |                                                           |
  |                    balance_equations.h                    |
  |                                                           |
  -------------------------------------------------------------
  | The header file for the definition of the balance         |
  | equations used for standard as well as micromorphic       |
  | continuum mechanics.                                      |
  -------------------------------------------------------------
  | Note: N is the test function and eta is the interpolation |
  |       function. We use this notation since these symbols  |
  |       are not used by the micromorphic formulation.       |
  =============================================================
  | Dependencies:                                             |
  | Eigen: Open source matrix library available at            |
  |        eigen.tuxfamily.org.                               |
  =============================================================*/

#ifndef BALANCE_EQUATIONS_H
#define BALANCE_EQUATIONS_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <deformation_measures.h>

typedef Eigen::Matrix<double, 3, 12> Matrix_3x12;
typedef Eigen::Matrix<double, 9, 12> Matrix_9x12;
typedef Eigen::Matrix<double,27, 12> Matrix_27x12;

namespace balance_equations{

    //Forces from the balance of linear momentum
    void compute_internal_force(const double (&dNdx)[3], const Vector_9 &cauchy, double (&fint)[3]);

    void compute_internal_force(const double (&dNdx)[3], const std::vector<double> &cauchy, double (&fint)[3]);
    
    void compute_internal_force(const int &i, const double (&dNdx)[3], const Vector_9 &cauchy, double &fint_i);
    
    void compute_internal_force(const int &i, const double (&dNdx)[3], const std::vector<double> &cauchy, double &fint_i);

    void compute_body_force(const double &N, const double &density, const double (&b)[3], double (&fb)[3]);
    
    void compute_body_force(const int &i, const double &N, const double &density, const double (&b)[3], double &fb);

    void compute_kinematic_force(const double &N, const double &density, const double (&a)[3], double (&fkin)[3]);
    
    void compute_kinematic_force(const int &i, const double &N, const double &density, const double (&a)[3], double &fkin);

    //Stresses from the balance of first moment of momentum
    void compute_internal_couple(const double &N, const double (&dNdx)[3], const Vector_9 &cauchy, const Vector_9 &s, const Vector_27 &m, double (&cint)[9]);
    
    void compute_internal_couple(const int &i, const int &j, const double &N, const double (&dNdx)[3], const Vector_9 &cauchy, const Vector_9 &s, const Vector_27 &m, double &cint_ij);

    void compute_internal_couple(const double &N, const double (&dNdx)[3], const std::vector<double> &cauchy, const std::vector<double> &s, const std::vector<double> &m, double (&cint)[9]);
    
    void compute_internal_couple(const int &i, const int &j, const double &N, const double (&dNdx)[3], const std::vector<double> &cauchy, const std::vector<double> &s, const std::vector<double> &m, double &cint_ij);

    void compute_body_couple(const double &N, const double &density, const double (&l)[9], double (&cb)[9]);
    
    void compute_body_couple(const int &i, const int &j, const double &N, const double &density, const double (&l)[9], double &cb_ij);
    
    void compute_kinematic_couple(const double &N, const double &density, const double (&omega)[9], double (&ckin)[9]);
    
    void compute_kinematic_couple(const int &i, const int &j, const double &N, const double &density, const double (&omega)[9], double &ckin_ij);
    
    //Compute the jacobians of the balance of linear momentum w.r.t. the indicated displacement dof
    //Note: These seem to be incorrect but are retained to be the foundations for a Total-Lagrangian implementation
    void compute_internal_force_jacobian(const double &N, const double(&dNdx)[3], const double &eta, const double(&detadx)[3], const Matrix_9x9 &DcauchyDgrad_u, const Matrix_9x9 &DcauchyDphi, const Matrix_9x27 &DcauchyDgrad_phi, Matrix_3x12 &DfintDU);
    
    void compute_internal_force_jacobian(const int &i, const int &dof_num, const double &N, const double(&dNdx)[3], const double &eta, const double (&detadx)[3], const Matrix_9x9 &DcauchyDgrad_u, const Matrix_9x9 &DcauchyDphi, const Matrix_9x27 &DcauchyDgrad_phi, double &DfintDU_iA);
    
    void compute_internal_force_jacobian(const double &N, const double(&dNdx)[3], const double &eta, const double(&detadx)[3], const std::vector<std::vector<double>> &DcauchyDgrad_u, const std::vector<std::vector<double>> &DcauchyDphi, const std::vector<std::vector<double>> &DcauchyDgrad_phi, std::vector<std::vector<double>> &DfintDU);
    
    void compute_internal_force_jacobian(const int &i, const int &dof_num, const double &N, const double(&dNdx)[3], const double &eta, const double (&detadx)[3], const std::vector<std::vector<double>> &DcauchyDgrad_u, const std::vector<std::vector<double>> &DcauchyDphi, const std::vector<std::vector<double>> &DcauchyDgrad_phi, double &DfintDU_iA);

    //Compute the jacobians of the balance of first moment of momentum
    void compute_internal_couple_jacobian(const double &N, const double (&dNdx)[3], const double &eta, const double (&detadx)[3],
                                          const Matrix_9x9  &DcauchyDgrad_u, const Matrix_9x9  &DcauchyDphi, const Matrix_9x27 &DcauchyDgrad_phi,
                                          const Matrix_9x9  &DsDgrad_u,      const Matrix_9x9  &DsDphi,      const Matrix_9x27 &DsDgrad_phi,
                                          const Matrix_27x9 &DmDgrad_u,      const Matrix_27x9 &DmDphi,      const Matrix_27x27 &DmDgrad_phi,
                                          Matrix_9x12 &DcintDU);
                                          
    void compute_internal_couple_jacobian(const int &i, const int &j, const int &dof_num, const double &N, const double (&dNdx)[3], const double &eta, const double (&detadx)[3],
                                          const Matrix_9x9  &DcauchyDgrad_u, const Matrix_9x9  &DcauchyDphi, const Matrix_9x27 &DcauchyDgrad_phi,
                                          const Matrix_9x9  &DsDgrad_u,      const Matrix_9x9  &DsDphi,      const Matrix_9x27 &DsDgrad_phi,
                                          const Matrix_27x9 &DmDgrad_u,      const Matrix_27x9 &DmDphi,      const Matrix_27x27 &DmDgrad_phi,
                                          double &DcintDU_ijA);

    void compute_internal_couple_jacobian(const double &N, const double (&dNdx)[3], const double &eta, const double (&detadx)[3],
                                          const std::vector<std::vector<double>> &DcauchyDgrad_u, const std::vector<std::vector<double>> &DcauchyDphi, const std::vector<std::vector<double>> &DcauchyDgrad_phi,
                                          const std::vector<std::vector<double>> &DsDgrad_u,      const std::vector<std::vector<double>> &DsDphi,      const std::vector<std::vector<double>> &DsDgrad_phi,
                                          const std::vector<std::vector<double>> &DmDgrad_u,      const std::vector<std::vector<double>> &DmDphi,      const std::vector<std::vector<double>> &DmDgrad_phi,
                                          std::vector<std::vector<double>> &DcintDU);
    
    void compute_internal_couple_jacobian(const int &i, const int &j, const int &dof_num, const double &N, const double (&dNdx)[3], const double &eta, const double (&detadx)[3],
                                          const std::vector<std::vector<double>> &DcauchyDgrad_u, const std::vector<std::vector<double>> &DcauchyDphi, const std::vector<std::vector<double>> &DcauchyDgrad_phi,
                                          const std::vector<std::vector<double>> &DsDgrad_u,      const std::vector<std::vector<double>> &DsDphi,      const std::vector<std::vector<double>> &DsDgrad_phi,
                                          const std::vector<std::vector<double>> &DmDgrad_u,      const std::vector<std::vector<double>> &DmDphi,      const std::vector<std::vector<double>> &DmDgrad_phi,
                                          double &DcintDU_ijA);

    //Compute the jacobians of the balance of linear momentum (Current configuration)
    void compute_internal_force_jacobian(const double &N,        const double(&dNdx)[3],           const double &eta,             const double(&detadx)[3], 
                                         const double (&grad_u)[3][3], const double (&grad_phi)[9][3],
                                         const Vector_9 &cauchy, const Matrix_9x9 &DcauchyDgrad_u, const Matrix_9x9 &DcauchyDphi, const Matrix_9x27 &DcauchyDgrad_phi,
                                         Matrix_3x12 &DfintDU);

    void compute_internal_force_jacobian(const int &component,   const int &dof_num,
                                         const double &N,        const double(&dNdx)[3],           const double &eta,             const double(&detadx)[3], 
                                         const double (&grad_u)[3][3], const double (&grad_phi)[9][3],
                                         const Vector_9 &cauchy, const Matrix_9x9 &DcauchyDgrad_u, const Matrix_9x9 &DcauchyDphi, const Matrix_9x27 &DcauchyDgrad_phi,
                                         double &DfintDU_iA);

    void compute_internal_force_jacobian(const double &N, const double(&dNdx)[3], const double &eta, const double(&detadx)[3], 
                                         const double (&grad_u)[3][3], const double (&grad_phi)[9][3],
                                         const std::vector<double> &cauchy, const std::vector<std::vector<double>> &DcauchyDgrad_u, const std::vector<std::vector<double>> &DcauchyDphi,
                                         const std::vector<std::vector<double>> &DcauchyDgrad_phi,
                                         std::vector<std::vector<double>> &DfintDU);

    void compute_internal_force_jacobian(const int &component, const int &dof_num,
                                         const double &N, const double(&dNdx)[3], const double &eta, const double(&detadx)[3], 
                                         const double (&grad_u)[3][3], const double (&grad_phi)[9][3],
                                         const std::vector<double> &cauchy, const std::vector<std::vector<double>> &DcauchyDgrad_u, const std::vector<std::vector<double>> &DcauchyDphi,
                                         const std::vector<std::vector<double>> &DcauchyDgrad_phi,
                                         double &DfintDU_iA);

    //Compute the jacobians of the balance of the first moment of momentum (Current configuration)
    void compute_internal_couple_jacobian(const double &N,  const double (&dNdx)[3], const double &eta, const double (&detadx)[3],
                                          const double (&grad_u)[3][3], const double (&grad_phi)[9][3],
                                          const Matrix_9x9      &DcauchyDgrad_u,              const Matrix_9x9  &DcauchyDphi, const Matrix_9x27 &DcauchyDgrad_phi,
                                          const Matrix_9x9      &DsDgrad_u,                   const Matrix_9x9  &DsDphi,      const Matrix_9x27 &DsDgrad_phi,
                                          const Matrix_27x9     &DmDgrad_u,                   const Matrix_27x9 &DmDphi,      const Matrix_27x27 &DmDgrad_phi,
                                          Matrix_9x12 &DcintDU);
                                          
    void compute_internal_couple_jacobian(const int &component_i, const int &component_j, const int &dof_num,
                                          const double &N,  const double (&dNdx)[3], const double &eta, const double (&detadx)[3],
                                          const double (&grad_u)[3][3], const double (&grad_phi)[9][3],
                                          const Matrix_9x9      &DcauchyDgrad_u,              const Matrix_9x9  &DcauchyDphi, const Matrix_9x27 &DcauchyDgrad_phi,
                                          const Matrix_9x9      &DsDgrad_u,                   const Matrix_9x9  &DsDphi,      const Matrix_9x27 &DsDgrad_phi,
                                          const Matrix_27x9     &DmDgrad_u,                   const Matrix_27x9 &DmDphi,      const Matrix_27x27 &DmDgrad_phi,
                                          double &DcintDU_ijA);
                                          
    void compute_internal_couple_jacobian(const double &N, const double (&dNdx)[3], const double &eta, const double (&detadx)[3],
                                          const double (&grad_u)[3][3], const double (&grad_phi)[9][3],
                                          const std::vector<std::vector<double>> &DcauchyDgrad_u, const std::vector<std::vector<double>> &DcauchyDphi, const std::vector<std::vector<double>> &DcauchyDgrad_phi,
                                          const std::vector<std::vector<double>> &DsDgrad_u,      const std::vector<std::vector<double>> &DsDphi,      const std::vector<std::vector<double>> &DsDgrad_phi,
                                          const std::vector<std::vector<double>> &DmDgrad_u,      const std::vector<std::vector<double>> &DmDphi,      const std::vector<std::vector<double>> &DmDgrad_phi,
                                          std::vector<std::vector<double>> &DcintDU);
    
    void compute_internal_couple_jacobian(const int &i, const int &j, const int &dof_num,
                                          const double &N, const double (&dNdx)[3], const double &eta, const double (&detadx)[3],
                                          const double (&grad_u)[3][3], const double (&grad_phi)[9][3],
                                          const std::vector<std::vector<double>> &DcauchyDgrad_u, const std::vector<std::vector<double>> &DcauchyDphi, const std::vector<std::vector<double>> &DcauchyDgrad_phi,
                                          const std::vector<std::vector<double>> &DsDgrad_u,      const std::vector<std::vector<double>> &DsDphi,      const std::vector<std::vector<double>> &DsDgrad_phi,
                                          const std::vector<std::vector<double>> &DmDgrad_u,      const std::vector<std::vector<double>> &DmDphi,      const std::vector<std::vector<double>> &DmDgrad_phi,
                                          double &DcintDU_ijA);

    //The jacobians of u and phi w.r.t. the DOF vector (Total-Lagrangian)
    void construct_dgrad_udU(const double (&detadx)[3], SpMat &dgrad_udU);    

    void construct_dgrad_udU(const double (&detadx)[3], Matrix_9x12 &dgrad_udU);
    
    void construct_dphidU(const double &eta, SpMat &dphidU);

    void construct_dphidU(const double &eta, Matrix_9x12 &dphidU);
    
    void construct_dgrad_phidU(const double (&detadx)[3], SpMat &dgrad_phidU);

    void construct_dgrad_phidU(const double (&detadx)[3], Matrix_27x12 &dgrad_phidU);
    
    //The jacobians of u and phi w.r.t. the DOF vector (Current Configuration)
    void construct_dgrad_udU(const Matrix_3x3 &Finv, const double (&detadx)[3], Matrix_9x12 &dgrad_udU);
    void construct_dgrad_phidU(const double (&grad_phi)[9][3], const Matrix_3x3 &Finv, const double (&detadx)[3], Matrix_27x12 &dgrad_phidU);

    void map_eigen_to_vector(const Vector_9  &V,       std::vector<double> &v);
    void map_eigen_to_vector(const Vector_27 &V,       std::vector<double> &v);
    void map_eigen_to_vector(const Eigen::VectorXd &V, std::vector<double> &v);
    void map_eigen_to_vector(const Matrix_3x3   &M,    std::vector<std::vector<double>> &v);
    void map_eigen_to_vector(const Matrix_9x9   &M,    std::vector<std::vector<double>> &v);
    void map_eigen_to_vector(const Matrix_9x27  &M,    std::vector<std::vector<double>> &v);
    void map_eigen_to_vector(const Matrix_27x9  &M,    std::vector<std::vector<double>> &v);
    void map_eigen_to_vector(const Matrix_27x27 &M,    std::vector<std::vector<double>> &v);
    void map_eigen_to_vector(const Matrix_3x12  &M,    std::vector<std::vector<double>> &v);
    void map_eigen_to_vector(const Matrix_9x12  &M,    std::vector<std::vector<double>> &v);
    void map_eigen_to_vector(const Eigen::MatrixXd &M, std::vector<std::vector<double>> &v);

    void map_vector_to_eigen(const std::vector<double> &v, Vector_9  &V      );
    void map_vector_to_eigen(const std::vector<double> &v, Vector_27 &V      );
    void map_vector_to_eigen(const std::vector<double> &v, Eigen::VectorXd &V);
    void map_vector_to_eigen(const std::vector<std::vector<double>> &v, Matrix_3x3      &M);
    void map_vector_to_eigen(const std::vector<std::vector<double>> &v, Matrix_9x9      &M);
    void map_vector_to_eigen(const std::vector<std::vector<double>> &v, Matrix_9x27     &M);
    void map_vector_to_eigen(const std::vector<std::vector<double>> &v, Matrix_27x9     &M);
    void map_vector_to_eigen(const std::vector<std::vector<double>> &v, Matrix_27x27    &M);
    void map_vector_to_eigen(const std::vector<std::vector<double>> &v, Matrix_3x12     &M);
    void map_vector_to_eigen(const std::vector<std::vector<double>> &v, Matrix_9x12     &M);
    void map_vector_to_eigen(const std::vector<std::vector<double>> &v, Eigen::MatrixXd &M);
}

#endif
