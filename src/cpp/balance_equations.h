/*!============================================================
  |                                                           |
  |                    balance_equations.h                    |
  |                                                           |
  -------------------------------------------------------------
  | The header file for the definition of the balance         |
  | equations used for standard as well as micromorphic       |
  | continuum mechanics.                                      |
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
    
    void compute_internal_force(const int &i, const double (&dNdx)[3], const Vector_9 &cauchy, double &fint_i);

    void compute_body_force(const double &N, const double &density, const double (&b)[3], double (&fb)[3]);
    
    void compute_body_force(const int &i, const double &N, const double &density, const double (&b)[3], double &fb);

    void compute_kinematic_force(const double &N, const double &density, const double (&a)[3], double (&fkin)[3]);
    
    void compute_kinematic_force(const int &i, const double &N, const double &density, const double (&a)[3], double &fkin);

    //Stresses from the balance of first moment of momentum
    void compute_internal_couple(const double &N, const double (&dNdx)[3], const Vector_9 &cauchy, const Vector_9 &s, const Vector_27 &m, double (&cint)[9]);
    
    void compute_internal_couple(const int &i, const int &j, const double &N, const double (&dNdx)[3], const Vector_9 &cauchy, const Vector_9 &s, const Vector_27 &m, double &cint_ij);

    void compute_body_couple(const double &N, const double &density, const double (&l)[9], double (&cb)[9]);
    
    void compute_body_couple(const int &i, const int &j, const double &N, const double &density, const double (&l)[9], double &cb_ij);
    
    void compute_kinematic_couple(const double &N, const double &density, const double (&omega)[9], double (&ckin)[9]);
    
    void compute_kinematic_couple(const int &i, const int &j, const double &N, const double &density, const double (&omega)[9], double &ckin_ij);
    
    //Compute the jacobians of the balance of linear momentum w.r.t. the indicated displacement dof
    void compute_internal_force_jacobian(const double &N, const double(&dNdx)[3], const Matrix_9x9 &DcauchyDgrad_u, const Matrix_9x9 &DcauchyDphi, const Matrix_9x27 &DcauchyDgrad_phi, Matrix_3x12 &DfintDU);
    
    void compute_internal_force_jacobian(const int &i, const int &dof_num, const double &N, const double(&dNdx)[3], const Matrix_9x9 &DcauchyDgrad_u, const Matrix_9x9 &DcauchyDphi, const Matrix_9x27 &DcauchyDgrad_phi, double &DfintDU_iA);
    
    //Compute the jacobians of the balance of first moment of momentum
    void compute_internal_couple_jacobian(const double &N, const double (&dNdx)[3],
                                          const Matrix_9x9  &DcauchyDgrad_u, const Matrix_9x9  &DcauchyDphi, const Matrix_9x27 &DcauchyDgrad_phi,
                                          const Matrix_9x9  &DsDgrad_u,      const Matrix_9x9  &DsDphi,      const Matrix_9x27 &DsDgrad_phi,
                                          const Matrix_27x9 &DmDgrad_u,      const Matrix_27x9 &DmDphi,      const Matrix_27x27 &DmDgrad_phi,
                                          Matrix_9x12 &DcintDU);
    
    //The jacobians of u and phi w.r.t. the DOF vector
    void construct_dgrad_udU(const double (&dNdx)[3], SpMat &dgrad_udU);    
    
    void construct_dphidU(const double &N, SpMat &dphidU);
    
    void construct_dgrad_phidU(const double (&dNdx)[3], SpMat &dgrad_phidU);
}

#endif
